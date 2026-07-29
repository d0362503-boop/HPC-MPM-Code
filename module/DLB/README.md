# `module/DLB/` — Dynamic Load Balancing

> Runtime repartitioning of the background mesh across MPI ranks, driven by the
> material-point distribution. Currently wired into the standalone fluid MPM solver
> (`src_fluid/MPM`); solid and FSI keep `do_dlb = false`.

---

## 1. What It Does

MPM particle load follows the physics (e.g. a dam-break front sweeping across the
domain), so a static uniform partition drifts out of balance over time. DLB periodically

1. samples the local particle coordinates on every rank,
2. recomputes a per-rank background-element region that equalizes particle counts,
3. migrates particles (and their scheme-specific state) to the new owners,
4. rebuilds all partition-dependent data: overlap tables, mesh/connectivity, nodal
   volumes, boundary conditions, and the PETSc matrix layout.

The **rank topology is fixed**: the Cartesian rank grid `nxyr` never changes; only the
element cuts between neighboring ranks move. This keeps neighbor discovery local and
the change incremental.

## 2. Files and Entry Points

| File                     | Role                                                                                                |
| ------------------------ | --------------------------------------------------------------------------------------------------- |
| `module/DLB/mpm_dlb.h`   | `Region` struct + the six public functions                                                          |
| `module/DLB/mpm_dlb.cpp` | implementation; owns the single `g_current_regions` instance                                        |
| `module/mpi_data.cpp`    | `MaterialPoint::{DetermineParticleRank, DetermineDLBParticleRank, RebalanceDLBParticles, ApplyDLB}` |
| `module/bc.cpp`          | `BoundaryCondition::{CaptureGlobal*, RebuildLocal*}` — BC persistence across repartitions           |

Public API (all collective unless noted):

```cpp
void CollectCurrentRegions();                       // Allgather current aelemmin/aelemmax
const std::vector<Region> &CurrentRegions();        // read-only access (not collective)
int ComputeSampleSkip(std::size_t local_npts);      // strided-sampling interval
std::vector<std::array<double,3>> SelectSamples(coord);  // picks skip internally (collective)
std::vector<Region> ComputeDLBRegions(local_samples);  // gather → cut on root → Bcast
void UpdateDLBRegions(const std::vector<Region> &);    // apply + rebuild overlap metadata
```

## 3. Activation and Call Sites

- `do_dlb` is the second token of the first parameter line of `input.txt`
  (`<solswitch, do_dlb, rstflag, nlstep>`), read in
  `work/src_fluid/MPM/datain_para.cpp`. Solid/FSI hardcode `do_dlb = false` in their
  constructors.
- The main loop (`work/src_fluid/MPM/main.cpp`) replaces `MoveParticle()` with
  `ApplyDLB()` when `do_dlb && istep % iout == 0` (and `nprocs > 1`):

```cpp
const bool is_dlb_step = wp.do_dlb && (istep % iout == 0);
if (nprocs != 1) is_dlb_step ? wp.ApplyDLB() : wp.MoveParticle();
```

- `StabilizedMPM::ApplyDLB()` calls `MaterialPoint::ApplyDLB()` and then rebuilds the
  Navier–Stokes matrix structure (`NS_.BuildCrsMat(16)`); the implicit solid override
  does the same with `SM_.BuildCrsMat(9)`.

## 4. The `ApplyDLB()` Pipeline

`MaterialPoint::ApplyDLB()` (`module/mpi_data.cpp`) runs six stages in order:

```text
MoveParticle()                       // 1. ordinary neighbor migration + inflow,
                                     //    so sampling sees current positions
SelectSamples(coord)                 // 2. strided local samples; the interval
                                     //    is chosen internally by ComputeSampleSkip
ComputeDLBRegions(local_samples)     // 3. root cuts new regions, MPI_Bcast
RebalanceDLBParticles(regions)       // 4. global particle migration to new owners
UpdateDLBRegions(regions)            // 5. new mesh sizes + overlap tables + weights
BuildMesh(); BuildControlPoint();
MakNodalVol(); RebuildBC();   // 6. partition-dependent data
```

The physics class then rebuilds its solver system (stage 8, in the overrides).

## 5. Region Model and Invariants

```cpp
struct Region { std::array<int,3> elem_min, elem_max; };  // inclusive, global element indices
```

- **Tiling invariant:** the `nprocs` regions must be pairwise disjoint and cover the
  whole `xyelemw` mesh exactly. It holds *by construction* (§7) and is re-validated at
  runtime: `UpdateDLBRegions` bounds-checks every region, and the migration builders
  `MPI_Abort` on any particle that is covered by zero or two regions.
- **Single instance:** `g_current_regions` is defined only in `mpm_dlb.cpp` and is not
  declared in the header — outside code can only read it through the const-reference
  accessor `CurrentRegions()` and can only change it through `CollectCurrentRegions()`
  (contains the `MPI_Allgather`) and `UpdateDLBRegions()` (input comes from the
  `MPI_Bcast` in `ComputeDLBRegions`). This makes accidental rank-local edits a
  compile error instead of a silent desync.
- **Cross-rank consistency** of the content is therefore guaranteed by collectives,
  not by C++ linkage — every rank holds the same table because it was gathered or
  broadcast, not because memory is shared.
- `BuildMesh()` calls `CollectCurrentRegions()`, so the table is valid from start-up
  (initial uniform partition from `griddata`) through every DLB update.

## 6. Sampling

Goal: ~`kSamplesPerRank` (128) coordinate samples per rank, proportional to local load.

```cpp
// ComputeSampleSkip
sample_rate   = max(kMinSampleRate /*1e-4*/, 128 * nprocs / global_count)  // MPI_Allreduce
target        = ceil(local_count * sample_rate)
skip          = max(1, local_count / target)
```

`SelectSamples` then takes every `skip`-th particle — deterministic, O(num), and
density-faithful because every particle is equally likely to be sampled regardless of
rank. `GatherSamplesToRoot` packs `(x,y,z)` into a flat double array and uses
`MPI_Gather`/`MPI_Gatherv`; only the root keeps the global sample set.

## 7. Repartitioning: Recursive Coordinate Bisection

`ComputeDLBRegions` (root only) cuts the fixed rank grid `nxyr = (nx, ny, nz)` in z,
then y, then x:

```text
z_cuts = CutsFromSamples(all samples, dir=z, 0 .. xyelemw[2]-1, nz)
for each z-slice:
    y_cuts = CutsFromSamples(samples inside the z-slice, dir=y, ...)
    for each y-slice:
        x_cuts = CutsFromSamples(samples inside the z,y-slice, dir=x, ...)
        region(ix,iy,iz) = [x_cuts[ix] .. x_cuts[ix+1]-1] × ... 
        rank = ix + nx*iy + nx*ny*iz
```

`CutsFromSamples` per direction:

- sorts the sample coordinates and picks the quantile boundary between child `part-1`
  and `part` — the midpoint of the two samples straddling `size*part/parts`, converted
  to an element index with `round((x - xyminw)/dxy)`;
- **clamps every cut** to `[prev_cut+1, hi+1-(parts-part)]`, so each child gets at
  least one element and the children tile `[lo, hi]` without gaps or overlaps;
- falls back to the uniform cut `lo + (hi+1-lo)*part/parts` when the parent box has no
  samples (e.g. an empty region) — balance degrades gracefully, correctness never does;
- `FilterSamples` restricts the sample set to a child box by *element index*
  (`LocateGlobalElement`), so every sample belongs to exactly one child.

Result is broadcast (`MPI_Bcast`) so all ranks hold the identical region table.

## 8. Particle Redistribution

Two different migration builders exist; DLB requires the global one:

|                | `DetermineParticleRank`         | `DetermineDLBParticleRank`             |
| -------------- | ------------------------------- | -------------------------------------- |
| use            | per-step drift (`MoveParticle`) | repartition (`RebalanceDLBParticles`)  |
| search         | neighbor regions `naid` only    | **all** ranks                          |
| aborts if      | particle crossed a non-neighbor | regions overlap / don't cover the mesh |
| count exchange | pair-wise Isend/Irecv           | `MPI_Alltoall`                         |

Because a repartition can move a particle across several old regions at once, the
neighbor restriction would lose particles — hence the global O(num × nprocs) owner
search at DLB time. A peer is registered when `send>0 || recv>0`, so sender and
receiver peer lists are always symmetric, and `PointVarSendrecv` counts match
pair-wise. Out-of-mesh particles are dropped via `nrmp` as usual.

`RebalanceDLBParticles` then updates `num += nrps - nmps - nrmp` and calls the virtual
`MigrateParticleData()`, which reuses the same `PointVarComm` machinery as ordinary
migration — scalars, vectors, and the TPIC/APIC matrix state all travel.

Note the ordering: particles are redistributed **while the old mesh data is still in
place**, because destinations are computed from global element indices
(`LocateGlobalElement` + region bounds), not from local mesh structures.

## 9. Region and Metadata Update

`UpdateDLBRegions` applies the new partition and rebuilds everything that depends on
it, in this order:

1. mesh extents: `aelemmin/aelemmax`, `xyelem/xynode/xynodec`, `xymin/xymax`,
   `nelem/node/nodec`, and the DOF offsets `nu..npc`;
2. overlap tables `naid/nsbc/nsbl/nsubc/nsubl`: recomputed geometrically — a peer is
   any rank whose region intersects this rank's **extended** control-point range
   `[aelemmin, aelemmax+idimc]` (for `nsubc`) or node range `[aelemmin, aelemmax+1]`
   (for `nsubl`);
3. ghost weights `dbc/dbl = 1/share_count`, matching the `dbn` semantics of the
   griddata path, so `NodeVarComm`'s sum-then-weight averaging is unchanged.

The caller then rebuilds derived data: `BuildMesh()` (also re-collects
`g_current_regions`), `BuildControlPoint()` (`ncc`, `xyc`), `MakNodalVol()` (`nvol`,
with the new overlap tables), `RebuildBC()` (§10), and the physics
class rebuilds its `CrsMat` — with `use_petsc`, `BuildCrsMat` first calls
`ResetPetscSolver()` to destroy the old `Mat/Vec/KSP/lgmap/scatter` (all handles are
`nullptr`-initialized, so the first call is a safe no-op).

Nodal *state* is intentionally not preserved: the fluid MPM rebuilds all nodal fields
from particles at the next `Particle2Node()` anyway (arrays are re-`VectorAssign`ed
per step with the current `nodec`).

## 10. Boundary-Condition Persistence

Local BC IDs are meaningless after a repartition, so BCs are cached as **global IDs**
at input time (`InputBCData`, when `do_dlb = 1`):

- `CaptureGlobalControlPointBC()` converts each local control-point BC to its global
  ID (using `xynodecw`); `CaptureGlobalElementBC()` does the same for inflow
  element BCs (using `xyelemw`);
- `CacheGlobalEntries` `MPI_Allgatherv`s every rank's list and **de-duplicates**
  (overlap control points are stored by every owning rank). Duplicated entries must
  carry the same prescribed value within `1e-12`, otherwise the run aborts — this
  catches inconsistent partitioned input instead of silently picking a side;
- after a repartition, `RebuildLocalControlPointBC()` / `RebuildLocalElementBC()`
  filter the (rank-identical) global cache against the new region — no communication
  needed at rebuild time.

`StabilizedMPM::RebuildBC()` additionally resets the inflow counters
(`InitializeInflowBC`, `infbc_isfilled = false`), so the next inflow check
re-evaluates boundary-cell fill under the new partition.

## 11. Other Integration Points

- `StabilizedMPM::MoveParticle()` and the solid equivalent pass
  `mpm_dlb::CurrentRegions()` to `DetermineParticleRank`, so per-step migration always
  uses the *current* partition without touching the update machinery.
- `PairwiseRepulsiveParticleShifting` builds its ghost set by intersecting each
  particle's support ball with the peers' physical region bounds read from
  `CurrentRegions()` — correct for non-uniform regions, unlike the old Cartesian rank
  arithmetic. The regions reference is used transiently inside the call (converted to
  local bounds immediately), so it is safe across DLB updates.
- `CrsMat` ownership (`BuildOwnershipMetadata`) gathers `aelemmin/aelemmax` at matrix
  build time, so after `UpdateDLBRegions` + `BuildCrsMat` the PETSc ownership/lgmaps
  are consistent with the new regions automatically.

## 12. Failure Checks (all `MPI_Abort`)

| site                       | condition                                                                             |
| -------------------------- | ------------------------------------------------------------------------------------- |
| `CutsFromSamples`          | invalid dir/range/cell size, or fewer elements than splits                            |
| `ComputeDLBRegions`        | `nxyr` product ≠ `nprocs`; mesh smaller than rank grid                                |
| `UpdateDLBRegions`         | region count ≠ `nprocs`; region outside global mesh                                   |
| `DetermineParticleRank`    | bad region count/peer rank; overlapping peers; particle crossed a non-neighbor region |
| `DetermineDLBParticleRank` | region overlap; regions don't cover the mesh                                          |
| `CacheGlobalEntries`       | conflicting prescribed values on one global BC node                                   |
| `FilterSamples`            | sample outside the global mesh                                                        |

## 13. Known Limitations

- **Restart:** `*_re.txt` files hold the DLB-adapted per-rank distribution, but a
  restart always boots into the uniform `griddata` partition. The first
  `MoveParticle()` (neighbor-only search) then aborts on particles whose uniform
  owner is not a neighbor. Workaround until this is handled in code: after
  `RestartInput()` + `BuildMesh()`, call
  `RebalanceDLBParticles(mpm_dlb::CurrentRegions())` once (its global search handles
  any source distribution).
- **Cost:** owner search is O(num × nprocs) per DLB; acceptable at `iout` cadence. If
  it ever shows up in profiles, build an element→owner lookup table once per
  repartition for O(1) per particle.
- **Balance quality** is limited by the sample count (≈128/rank) and the
  element-granular cuts; expect approximate, not exact, load balance.
- **Fixed topology:** `nxyr` cannot change at runtime, and `nprocs` must always equal
  the partition topology stored in `griddata`.
