# MPM-Code Solid Mechanics Quick Reference

> Full constitutive model docs: `module/solid/constitutive_model.md` (336 lines).  
> Full implicit solver docs: `module/solid/implicit/README.md` (141 lines).

## Constitutive Model Dispatcher

`ConstitutiveModel::UpdateStress(int model_type, ...)` selects the model:

| `model_type` | Model | `mat_prop` |
|:------------:|:------|:-----------|
| `-1` | Rigid body | — |
| `0` | Linear elasticity | `[E, ν]` |
| `1` | Neo-Hookean | `[E, ν]` |
| `2` | St. Venant-Kirchhoff | `[E, ν]` |
| `3` | Mooney-Rivlin | `[C₁, C₂, K]` |

## Key Signatures

```cpp
class ConstitutiveModel {
public:
    bool implicit_flag = false;                    // true → also assemble stif_mat
    std::array<std::array<double, 6>, 6> stif_mat; // 6x6 tangent stiffness (Voigt)

    void UpdateStress(int model_type,
                      std::vector<double>& stress,  // 6-comp Voigt, in-place update
                      const std::array<std::array<double, 3>, 3>& F,
                      const std::array<std::array<double, 3>, 3>& dF,
                      const std::vector<double>& mat_prop,
                      double jac, double jac_bar);  // jac_bar = F-bar corrected J
};
```

## Critical Rules

- **Voigt stress order:** `{σ_xx, σ_yy, σ_zz, σ_yz, σ_xz, σ_xy}`
- **`stif_mat` is per-particle temporary storage.** Overwritten on every `UpdateStress` call. `ComputeTangentModulus` must read it **immediately** after the call.
- **All models assume 3D isotropic materials.**
- **F-bar correction:** `jac_bar` equals `jac` when F-bar is disabled. Neo-Hookean uses `ln(jac_bar)` for the volumetric term.
- **`idm[3][3]` maps tensor indices to Voigt:**
  - `idm[0][0]=0, idm[1][1]=1, idm[2][2]=2`
  - `idm[1][2]=3, idm[0][2]=4, idm[0][1]=5`

## Implicit vs Explicit

| Mode | Behavior |
|------|----------|
| `implicit_flag == false` | Stress update only; `stif_mat` not computed |
| `implicit_flag == true` | Stress update + tangent stiffness assembly into `stif_mat` |

## Typical Call Flow

```cpp
// Inside SolidMaterialPointBase::UpdateConstitutiveModel
cm_.implicit_flag = this->implicit_flag;
cm_.UpdateStress(model_type, sts, F, dF, mat_prop, jac, jac_bar);
// sts updated in-place
// If implicit, cm_.stif_mat holds the particle's tangent stiffness
```

## Solid Solver Entry Points

| Physics | Key Files |
|---------|-----------|
| Explicit solid | `module/solid/explicit/solve_solid_explicit.cpp` |
| Implicit solid | `module/solid/implicit/solve_solid_implicit.cpp` |
| F-bar projection | `module/solid/explicit/fbar_projection.cpp` |

## Implicit Solver

| File | Role |
|------|------|
| `implicit_mpm_solid.h` | `SolidMaterialPoint` class declaration, wires `solid_.owner_ = this`, inline singleton `sp` |
| `solve_solid_implicit.cpp` | NR loop, system assembly, tangent-modulus computation |
| `var_trans_solid_implicit.cpp` | P2G, G2P, def-grad update, particle shifting |

### Solver Characteristics

- **Time integration:** Newmark-β / Generalized-α (`PredictNewmarkBetaVelAndAccel`, `CommitNodalKinematics`)
- **Nonlinear solve:** Newton–Raphson with `GPBiCGSafe` inner linear solve
- **Matrix:** `solid_` is a `CrsMat` with `ndof = 3`, `owner_ = this`, `FEM_flag = false`
- **Constructor disables PETSc** (`use_petsc = false`) for the native solver path
- **BC hooks** injected through virtuals: `BuildPetscBCList`, `BCResidualSet`, `BCNRSet`
