# ConstitutiveModel Class Documentation

## 1. Overview

`ConstitutiveModel` encapsulates all constitutive model implementations for the MPM solid mechanics solver. It performs per-particle stress updates and, when running in implicit mode, assembles the material stiffness matrix. The class does not own particle state arrays; instead, it operates on a single particle's stress in-place.

## 2. Class Structure

```cpp
class ConstitutiveModel {
public:
    bool implicit_flag = false;
    std::array<std::array<double, 6>, 6> stif_mat;

    void UpdateStress(int model_type, std::vector<double>& stress,
                      const std::array<std::array<double, 3>, 3>& F,
                      const std::array<std::array<double, 3>, 3>& dF,
                      const std::vector<double>& mat_prop,
                      double jac, double jac_bar);

    void MatRigid(std::vector<double>& stress);
    void MatLinElast(std::vector<double>& stress,
                     const std::array<std::array<double, 3>, 3>& delta_def_grad,
                     const std::vector<double>& mat_prop);
    void MatNeoHookean(std::vector<double>& stress,
                       const std::array<std::array<double, 3>, 3>& FF,
                       const std::vector<double>& mat_prop,
                       double jac, double jac_bar);
    void MatStVk(std::vector<double>& stress,
                 const std::array<std::array<double, 3>, 3>& FF,
                 const std::vector<double>& mat_prop, double jac);
    void MatMooneyRivlin(std::vector<double>& stress,
                         const std::array<std::array<double, 3>, 3>& FF,
                         const std::vector<double>& mat_prop, double jac);
};
```

### 2.1 Member Variables

| Variable | Type | Description |
|----------|------|-------------|
| `implicit_flag` | `bool` | Implicit solver flag. When `true`, each constitutive model assembles `stif_mat` after updating stress; when `false`, stiffness computation is skipped for efficiency. |
| `stif_mat` | `std::array<std::array<double,6>,6>` | 6x6 material tangent stiffness matrix $\mathbf{C}$ in Voigt notation, used for tangent modulus assembly in implicit Newton-Raphson iterations. The index mapping follows the global array `idm[3][3]`. |

### 2.2 Public Interface

#### `UpdateStress`

Dispatcher that selects the constitutive model based on `model_type`.

| `model_type` | Constitutive Model |
|:------------:|:-------------------|
| `-1` | Rigid body |
| `0` | Linear elasticity |
| `1` | Neo-Hookean hyperelasticity |
| `2` | St. Venant-Kirchhoff hyperelasticity |
| `3` | Mooney-Rivlin hyperelasticity |

**Parameters**:
- `stress` — Input/output, 6-component Cauchy stress vector $\{\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz}, \sigma_{xy}\}$
- `F` — Current deformation gradient $\mathbf{F}$
- `dF` — Incremental deformation gradient $\Delta\mathbf{F}$
- `mat_prop` — Material parameter array; meaning varies by model
- `jac` — Determinant of current deformation gradient, $J = \det(\mathbf{F})$
- `jac_bar` — F-bar corrected determinant $\bar{J}$ (equals $J$ when F-bar is disabled)

---

## 3. Constitutive Model Details

### 3.1 Rigid Body (`MatRigid`)

**Usage**: Rigid boundary conditions or inactive particles.

**Stress update**:
$$
\boldsymbol{\sigma} = \mathbf{0}
$$

**Implicit stiffness** (if `implicit_flag == true`):
$$
\mathbf{C} = \mathbf{0}
$$

---

### 3.2 Linear Elasticity (`MatLinElast`)

**Material parameters**:
- `mat_prop[0]` — Young's modulus $E$
- `mat_prop[1]` — Poisson's ratio $\nu$

**Elasticity parameters**:
$$
C_0 = \frac{E(1-\nu)}{(1-2\nu)(1+\nu)}, \quad
C_1 = \frac{E\nu}{(1-2\nu)(1+\nu)}, \quad
C_2 = \frac{E}{2(1+\nu)}
$$

**Incremental strain** (engineering strain from $\Delta\mathbf{F}$):
$$
\begin{aligned}
\Delta\varepsilon_{xx} &= \Delta F_{xx} - 1 \\
\Delta\varepsilon_{yy} &= \Delta F_{yy} - 1 \\
\Delta\varepsilon_{zz} &= \Delta F_{zz} - 1 \\
\Delta\varepsilon_{yz} &= \Delta F_{yz} + \Delta F_{zy} \\
\Delta\varepsilon_{xz} &= \Delta F_{xz} + \Delta F_{zx} \\
\Delta\varepsilon_{xy} &= \Delta F_{xy} + \Delta F_{yx}
\end{aligned}
$$

**Stress increment** (Voigt notation):
$$
\begin{aligned}
\Delta\sigma_{xx} &= C_0\Delta\varepsilon_{xx} + C_1\Delta\varepsilon_{yy} + C_1\Delta\varepsilon_{zz} \\
\Delta\sigma_{yy} &= C_1\Delta\varepsilon_{xx} + C_0\Delta\varepsilon_{yy} + C_1\Delta\varepsilon_{zz} \\
\Delta\sigma_{zz} &= C_1\Delta\varepsilon_{xx} + C_1\Delta\varepsilon_{yy} + C_0\Delta\varepsilon_{zz} \\
\Delta\sigma_{yz} &= C_2\Delta\varepsilon_{yz} \\
\Delta\sigma_{xz} &= C_2\Delta\varepsilon_{xz} \\
\Delta\sigma_{xy} &= C_2\Delta\varepsilon_{xy}
\end{aligned}
$$

**Stress update**:
$$
\boldsymbol{\sigma}^{n+1} = \boldsymbol{\sigma}^{n} + \Delta\boldsymbol{\sigma}
$$

**Implicit material stiffness matrix** (6x6):
$$
\mathbf{C} = \begin{bmatrix}
C_0 & C_1 & C_1 & 0 & 0 & 0 \\
C_1 & C_0 & C_1 & 0 & 0 & 0 \\
C_1 & C_1 & C_0 & 0 & 0 & 0 \\
0 & 0 & 0 & C_2 & 0 & 0 \\
0 & 0 & 0 & 0 & C_2 & 0 \\
0 & 0 & 0 & 0 & 0 & C_2
\end{bmatrix}
$$

---

### 3.3 Neo-Hookean Hyperelasticity (`MatNeoHookean`)

**Material parameters**:
- `mat_prop[0]` — Young's modulus $E$
- `mat_prop[1]` — Poisson's ratio $\nu$

**Lame constants**:
$$
\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}, \quad
\mu = \frac{E}{2(1+\nu)}
$$

**Left Cauchy-Green tensor**:
$$
\mathbf{b}^e = \mathbf{F}\mathbf{F}^\mathrm{T}
$$

**Cauchy stress** (using F-bar corrected volumetric deformation):
$$
\boldsymbol{\sigma} = \frac{1}{J}\left[\mu\left(\mathbf{b}^e - \mathbf{I}\right) + \lambda\ln(\bar{J})\,\mathbf{I}\right]
$$

**Implicit material stiffness** (fourth-order tensor $\mathbf{C}$ in Voigt form):
$$
C_{ijkl} = \lambda'\,\delta_{ij}\delta_{kl} + \mu'\left(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}\right)
$$

where:
$$
\lambda' = \frac{\lambda}{J}, \quad
\mu' = \frac{\mu - \lambda\ln(\bar{J})}{J}
$$

---

### 3.4 St. Venant-Kirchhoff Hyperelasticity (`MatStVk`)

**Material parameters**:
- `mat_prop[0]` — Young's modulus $E$
- `mat_prop[1]` — Poisson's ratio $\nu$

**Lame constants** (same as above):
$$
\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}, \quad
\mu = \frac{E}{2(1+\nu)}
$$

**Left Cauchy-Green tensor and its square**:
$$
\mathbf{b}^e = \mathbf{F}\mathbf{F}^\mathrm{T}, \quad
(\mathbf{b}^e)^2 = \mathbf{b}^e\mathbf{b}^e
$$

**Green-Lagrange strain trace**:
$$
\mathrm{tr}(\mathbf{E}_\mathrm{GL}) = \frac{1}{2}\left(\mathrm{tr}(\mathbf{b}^e) - 3\right)
$$

**Cauchy stress**:
$$
\boldsymbol{\sigma} = \frac{1}{J}\left[\mu(\mathbf{b}^e)^2 + \left(\lambda\,\mathrm{tr}(\mathbf{E}_\mathrm{GL}) - \mu\right)\mathbf{b}^e\right]
$$

**Implicit material stiffness** (fourth-order tensor):
$$
C_{ijkl} = \lambda'\,b^e_{ij}\,b^e_{kl} + \mu'\left(b^e_{ik}\,b^e_{jl} + b^e_{il}\,b^e_{jk}\right)
$$

where:
$$
\lambda' = \frac{\lambda}{J}, \quad
\mu' = \frac{\mu}{J}
$$

---

### 3.5 Mooney-Rivlin Hyperelasticity (`MatMooneyRivlin`)

**Material parameters**:
- `mat_prop[0]` — $C_1$
- `mat_prop[1]` — $C_2$
- `mat_prop[2]` — Bulk modulus $K$

**Isochoric deformation gradient**:
$$
\mathbf{F}_\mathrm{iso} = J^{-1/3}\,\mathbf{F}
$$

**Isochoric left Cauchy-Green tensor**:
$$
\mathbf{b}_\mathrm{iso} = \mathbf{F}_\mathrm{iso}\mathbf{F}_\mathrm{iso}^\mathrm{T}
$$

**First and second principal invariants**:
$$
I_1 = \mathrm{tr}(\mathbf{b}_\mathrm{iso}), \quad
I_2 = \frac{1}{2}\left(I_1^2 - \mathrm{tr}(\mathbf{b}_\mathrm{iso}^2)\right)
$$

**Volumetric contribution**:
$$
h_p = \frac{K\ln J}{J}
$$

**Cauchy stress** (deviatoric + volumetric):
$$
\sigma_{ij} = \frac{2}{J}\left[(C_1 + C_2 I_1)\,b_{\mathrm{iso},ij} - C_2\,b_{\mathrm{iso},ij}^2\right] + \delta_{ij}\left[h_p - \frac{2}{3J}\left(C_1 I_1 + 2C_2 I_2\right)\right]
$$

**Implicit material stiffness** (fourth-order tensor decomposed into three parts):

$$
\mathbf{C} = \mathbf{C}^\mathrm{vol} + \mathbf{C}^\mathrm{dev} + \mathbf{C}^\mathrm{w}
$$

**Volumetric part**:
$$
C^\mathrm{vol}_{ijkl} = \frac{K}{J}\delta_{ij}\delta_{kl} - 2h_p\,\mathbb{I}^\mathrm{sym}_{ijkl}
$$

**Deviatoric geometric part**:
$$
C^\mathrm{dev}_{ijkl} = -\frac{2}{3}\left(\sigma^\mathrm{dev}_{ij}\delta_{kl} + \delta_{ij}\sigma^\mathrm{dev}_{kl}\right)
$$

**Mooney-Rivlin material part $\mathbf{C}^\mathrm{w}$**:

Let:
$$
\mathbf{A} = I_1\,\mathbf{b}_\mathrm{iso} - \mathbf{b}_\mathrm{iso}^2
$$

Then:
$$
\begin{aligned}
C^\mathrm{w}_{ijkl} &=
\underbrace{\frac{4}{3J}(C_1 I_1 + 2C_2 I_2)\left(\mathbb{I}^\mathrm{sym}_{ijkl} - \frac{1}{3}\delta_{ij}\delta_{kl}\right)}_{\text{term1}} \\
&+ \underbrace{\frac{4C_2}{J}\left(b_{\mathrm{iso},ij}\,b_{\mathrm{iso},kl} - \frac{1}{2}\left(b_{\mathrm{iso},ik}\,b_{\mathrm{iso},jl} + b_{\mathrm{iso},il}\,b_{\mathrm{iso},jk}\right)\right)}_{\text{term2}} \\
&+ \underbrace{-\frac{4C_2}{3J}\left(A_{ij}\delta_{kl} + \delta_{ij}A_{kl}\right)}_{\text{term3}} \\
&+ \underbrace{\frac{8C_2 I_2}{9J}\delta_{ij}\delta_{kl}}_{\text{term4}}
\end{aligned}
$$

where the fourth-order symmetric identity tensor is:
$$
\mathbb{I}^\mathrm{sym}_{ijkl} = \frac{1}{2}\left(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}\right)
$$

---

## 4. Storage Format

### 4.1 Stress Vector (Voigt Notation)

$$
\boldsymbol{\sigma} = \begin{bmatrix} \sigma_{xx} & \sigma_{yy} & \sigma_{zz} & \sigma_{yz} & \sigma_{xz} & \sigma_{xy} \end{bmatrix}^\mathrm{T}
$$

### 4.2 Stiffness Matrix Index Mapping

The fourth-order tensor $C_{ijkl}$ is mapped to a $6\times6$ Voigt matrix via the global array `idm[3][3]`:

| $(i,j)$ | `idm[i][j]` |
|:-----:|:---------:|
| $(0,0)$ | 0 |
| $(1,1)$ | 1 |
| $(2,2)$ | 2 |
| $(1,2)$ | 3 |
| $(0,2)$ | 4 |
| $(0,1)$ | 5 |

That is:
$$
\mathbf{C}_\mathrm{Voigt}\big[\,\mathrm{idm}[i][j]\,\big]\big[\,\mathrm{idm}[k][l]\,\big] = C_{ijkl}
$$

---

## 5. Usage Example

Typical call flow (implicit MPM):

```cpp
// Inside SolidMaterialPointBase::UpdateConstitutiveModel
cm_.implicit_flag = this->implicit_flag;   // Sync implicit flag
cm_.UpdateStress(model_type, sts, F, dF, mat_prop, jac, jac_bar);
// sts has been updated in-place
// If implicit_flag == true, cm_.stif_mat holds the current particle's tangent stiffness
```

**Important notes**:
- `stif_mat` is temporary per-particle storage; it is overwritten on each `UpdateStress` call.
- `ComputeTangentModulus` must read `cm_.stif_mat` immediately after `UpdateStress`; otherwise the data becomes stale.
- All models assume 3D isotropic materials.
