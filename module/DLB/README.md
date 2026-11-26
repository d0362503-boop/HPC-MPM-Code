# `module/DLB/` — Dynamic Load Balancing

> Runtime repartitioning of the background mesh across MPI ranks, driven by the
> material-point distribution. The calling material owner selects whether and
> when to enable repartitioning through `do_dlb` and `ApplyDLB()`.

## Purpose

DLB adjusts background-element cuts to follow particle load, then migrates
particles and rebuilds partition-dependent data. The Cartesian rank grid
`nxyr` stays fixed; only the cuts move.

## File Inventory and Entry Points

| File                     | Role                                                                                                |
| ------------------------ | --------------------------------------------------------------------------------------------------- |
| `module/DLB/mpm_dlb.h`   | `Region`, sampling, region update, and diagnostic interfaces                                         |
| `module/DLB/mpm_dlb.cpp` | implementation; reads and writes the regions stored in `Partition::regions`                         |
| `module/mpi_data.cpp`    | `MaterialPoint::{DetermineParticleRank, DetermineDLBParticleRank, RebalanceDLBParticles, ApplyDLB}` |
| `module/bc.cpp`          | `BoundaryCondition::{CaptureGlobal*, RebuildLocal*}` — BC persistence across repartitions           |

Public API (all collective unless noted):

```cpp
void CollectCurrentRegions(const Mesh &, Partition &);   // Allgather aelemmin/aelemmax into part.regions
int ComputeSampleSkip(std::size_t local_npts);           // strided-sampling interval
std::vector<std::array<double,3>> SelectSamples(coord);  // picks skip internally (collective)
std::vector<Region> ComputeDLBRegions(const Grid &, local_samples);  // gather → cut on root → Bcast
void UpdateDLBRegions(const std::vector<Region> &, Grid &);          // apply + rebuild overlap metadata
```

## Activation and Call Sites

- The caller sets `do_dlb` and selects a collective repartitioning cadence.
- One possible caller policy replaces ordinary migration at output steps
  when DLB is enabled and more than one rank participates:

```cpp
const bool is_dlb_step = wp.do_dlb && (ts.istep % params.iout == 0);
if (nprocs != 1) is_dlb_step ? wp.ApplyDLB() : wp.MoveParticle();
```

- `StabilizedMPM::ApplyDLB()` calls `MaterialPoint::ApplyDLB()` and then rebuilds the
  Navier–Stokes matrix structure (`NS_.BuildCrsMat(16)`); the implicit solid override
  does the same with `SM_.BuildCrsMat(9)`.

## The `ApplyDLB()` Pipeline

`MaterialPoint::ApplyDLB()` (`module/mpi_data.cpp`) runs six stages in order:

```text
MoveParticle()                       // 1. ordinary neighbor migration + inflow,
                                     //    so sampling sees current positions
SelectSamples(coord)                 // 2. strided local samples; the interval
                                     //    is chosen internally by ComputeSampleSkip
ComputeDLBRegions(grid, local_samples)  // 3. root cuts new regions, MPI_Bcast
RebalanceDLBParticles(regions)          // 4. global particle migration to new owners
UpdateDLBRegions(regions, grid)         // 5. new mesh sizes + overlap tables + weights
BuildMesh(); BuildControlPoint();
ComputeNodalVol(); RebuildBC();         // 6. partition-dependent data
```

The physics class then rebuilds its solver system after these six stages.

## Region Model and Invariants

```cpp
struct Region { std::array<int,3> elem_min, elem_max; };  // inclusive, global element indices
```

- **Tiling invariant:** the `nprocs` regions must be pairwise disjoint and cover the
  whole `xyelemw` mesh exactly. It holds *by construction* (see [Repartitioning](#repartitioning-recursive-coordinate-bisection)) and is re-validated at
  runtime: `UpdateDLBRegions` bounds-checks every region, and the migration builders
  `MPI_Abort` on any particle that is covered by zero or two regions.
- The regions live in `Partition::regions` of the rank's `Grid`. Other code reads
  them there; only `CollectCurrentRegions()` and `UpdateDLBRegions()` write them,
  through collective gather/broadcast paths.
- Each rank holds its own copy; collectives establish consistency.
  `InputParaGriddata()` collects the initial regions from `griddata`, and
  `UpdateDLBRegions()` refreshes them after DLB.

## Sampling

The sampling-rate target is `kSamplesPerRank = 256` times the rank count globally,
distributed approximately in proportion to each rank's particle count. The minimum
rate and integer stride can increase the actual sample count.

```cpp
// ComputeSampleSkip
sample_rate   = max(kMinSampleRate /*1e-4*/, 256 * nprocs / global_count)  // MPI_Allreduce
target        = ceil(local_count * sample_rate)
skip          = max(1, local_count / target)
```

`SelectSamples` then takes every `skip`-th particle. This is deterministic strided
sampling, not random sampling; particle ordering can bias the sample distribution.
`GatherSamplesToRoot` packs `(x,y,z)` into a flat double array and uses
`MPI_Gather`/`MPI_Gatherv`; only the root keeps the global sample set.

## Repartitioning: Recursive Coordinate Bisection

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

## Particle Redistribution

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

## Region and Metadata Update

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

The caller rebuilds mesh/control points, nodal volumes and BCs, then the
physics matrix. `BuildCrsMat()` resets the old PETSc objects, including compact
maps and retained equations; old numbering and KSP state cannot survive a
partition-layout change. See the [shared solver guide](../solver/README.md).

Nodal *state* is intentionally not preserved: the fluid MPM rebuilds all nodal fields
from particles at the next `Particle2Node()` anyway (arrays are re-`VectorAssign`ed
per step with the current `nodec`).

## Boundary-Condition Persistence

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

## Integration

- Ordinary fluid/solid migration reads `Partition::regions`.
- `PairwiseRepulsivePST` intersects particle support balls with those regions
  to select MPI ghosts, including non-uniform partitions.
- `CrsMat` rebuilds ownership from the new `aelemmin/aelemmax` bounds.

## Failure Checks (All `MPI_Abort`)

| site                       | condition                                                                             |
| -------------------------- | ------------------------------------------------------------------------------------- |
| `CutsFromSamples`          | invalid dir/range/cell size, or fewer elements than splits                            |
| `ComputeDLBRegions`        | `nxyr` product ≠ `nprocs`; mesh smaller than rank grid                                |
| `UpdateDLBRegions`         | region count ≠ `nprocs`; region outside global mesh                                   |
| `DetermineParticleRank`    | bad region count/peer rank; overlapping peers; particle crossed a non-neighbor region |
| `DetermineDLBParticleRank` | region overlap; regions don't cover the mesh                                          |
| `CacheGlobalEntries`       | conflicting prescribed values on one global BC node                                   |
| `FilterSamples`            | sample outside the global mesh                                                        |

## Known Limitations

- **Restart:** `*_re.txt` files hold the DLB-adapted per-rank distribution, but a
  restart always boots into the uniform `griddata` partition. The first
  `MoveParticle()` (neighbor-only search) then aborts on particles whose uniform
  owner is not a neighbor. Workaround until this is handled in code: after
  `RestartInput()` + `BuildMesh()`, call
  `RebalanceDLBParticles(grid.part.regions)` once (its global search handles
  any source distribution).
- **Cost:** DLB owner search is O(num × nprocs) per repartition.
- **Balance quality** is limited by the sampling target (256 per rank on average), particle ordering, and the
  element-granular cuts; expect approximate, not exact, load balance.
- **Fixed topology:** `nxyr` cannot change at runtime, and `nprocs` must always equal
  the partition topology stored in `griddata`.
