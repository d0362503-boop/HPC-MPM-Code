# MPM-Code Fluid Solver Quick Reference

> Full documentation: `module/fluid/FEM/README.md` (128 lines).

## Overview

Stabilized FEM fluid solver for incompressible Navier–Stokes with SUPG/PSPG stabilization and generalized-α time integration. Two `CrsMat` objects:

| Matrix | DOF | Purpose |
|--------|-----|---------|
| `NS_` | 4 (u, v, w, p) | Momentum + continuity |
| `PF_` | 1 (φ) | Phase-field / volume fraction (prepared, currently unused in shown source) |

## File Responsibilities

| File | Role |
|------|------|
| `stabilized_fem.h` | Class declaration, inline constructor, global singleton `wfem` |
| `solve_ns_fem.cpp` | Time integration, stabilization coefficients, system assembly, linear solve |
| `phase_field.cpp` | Phase-field init, volume-fraction projection, material blending, liquid-volume calc |
| `fem_data_io.cpp` | Restart I/O, per-step field output, VTK nodal interpolation |

> **Note:** Fluid BC overrides (`BuildPetscBCList`, `BCResidualSet`) are **declared** in `stabilized_fem.h` but **implemented** in `src_fsi/bc_setting.cpp`.

## Class Architecture

```
MaterialPoint
    └── StabilizedFEM
            ├── NS_  : CrsMat (ndof = 4, owner_ = this)
            └── PF_  : CrsMat (ndof = 1, owner_ = this)
```

## Numerical Method Summary

- **Time integration:** Generalized-α with frozen advection velocity per step
- **Stabilization:** Per-element SUPG/PSPG (`tau1`) and shock capturing (`tau2`)
- **System assembly:** `AssembleNSSystem` accumulates mass, advection, diffusion, pressure, continuity, and shock-capturing terms into `NS_.amat`
- **Linear solve:** `SolveNSSystem` dispatches to PETSc KSP or native `GPBiCGAR`
- **Solution unpack:** `NS_.x_lhs` → `nvel` (u,v,w) and `npres` (p)

## Phase-Field / Two-Fluid

- `Particle2NodePhi` projects liquid/gas indicator onto mesh, clips to `[0, 1]`
- `SetPFDomain` blends element-wise `rhoe` / `rmue` from liquid/gas values using nodal `nphi`

## Critical Notes

- `stabilized_fem.h` declares virtual BC hooks, but **implementations live in `src_fsi/bc_setting.cpp`** — `crsmat.cpp` never includes `stabilized_fem.h`.
- `using namespace stabilizedfem;` appears in `.cpp` files only.
- Global singleton `wfem` is declared inline in `stabilized_fem.h`.
