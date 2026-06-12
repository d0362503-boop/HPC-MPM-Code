---
name: mpm-code-workflow
description: MPM-Code solver run procedure — change MPI topology, regenerate partitions, adjust steps, compile, execute, and verify.
---

# MPM-Code Workflow

## 1. Adjust MPI Topology

Edit `data/para_input_data.txt`:

```
<nxyr>
4  1  2     # nx, ny, nz decomposition (example; replace with your own)
```

- **Total MPI ranks** = `nx * ny * nz`.
- These three numbers are **examples** — replace them according to your actual core count and domain shape.
- After every change to `nxyr`, Step 2 (regenerate partitions) **must** be re-run.

## 2. Regenerate Partitions

```bash
cd data/
./divide/divide.exe
```

Generates `griddata0.txt ... griddataN.txt` and `pointdata0.txt ... pointdataN.txt` under `data/fsi/myrank_data_Turek/`.

> Must re-run after every `nxyr` change.

## 3. Adjust Steps

Edit `data/fsi/input.txt`:

```
<ista, iend, iout>
1  100  1     # start, end, output interval
```

## 4. Compile

```bash
cd MPM-Code/
make -j8
```

## 5. Run

```bash
cd MPM-Code/build
sh run.sh
```

- `run.sh` launches `./MPM` via `mpiexec`.
- Solver results are written to **stdout**; errors and warnings go to **stderr**.
- If you need a different rank count, edit `NP` inside `run.sh` before executing.  
  **Remember:** `NP` must equal the product of `nxyr` from Step 1 (`nx * ny * nz`).
- Always run from `build/` (contains `file.dat`).

## 6. Verify

There is **no automated unit-test framework** (no GTest, Catch2, etc.). Correctness is verified against standard benchmark problems whose inputs live in `data/` and whose reference outputs live in `res/`:

| Benchmark | Type | Location |
|-----------|------|----------|
| Turek FSI | FSI | `data/fsi/`, `res/Turek/` |
| Lid-driven cavity | Fluid | `data/fluid/`, `res/cavity/` |
| Dam break | Fluid | `data/fluid/`, `res/disk/` |
| Cook's membrane | Solid | `data/solid/`, `res/cook/` |
| Cantilever beam | Solid | `data/solid/`, `res/beam/` |

### Verification Workflow

1. Ensure steps 1–5 above completed successfully.
2. Convert raw text output to VTK with `vtk/makvtk.exe`.
3. Compare against reference solutions in `res/`.
