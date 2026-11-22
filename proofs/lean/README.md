# SplitNewton Lean proofs

Local Mathlib project for `SplitNewton.lean` (Lean `v4.37.0-rc1`).

## Layout

Standard background results (independent of the Split-Newton algorithm) live in their own modules
and are imported by the main file:

| Module | Contents |
| --- | --- |
| `SplitNewton/Spectral.lean` | `specRad`, Gelfand: `ρ(T) < 1 ⇒ ‖T^N‖ < 1` |
| `SplitNewton/Switched.lean` | switched iterations, `switched_local_convergence`, Ostrowski's `local_convergence` |
| `SplitNewton/JointSpectralRadius.lean` | `jointSpecRad`, `ρ̂ < 1 ⇒` uniform `N`-step contraction, `ρ̂({T}) ≤ ρ(T)` |
| `SplitNewton/Derivative.lean` | derivative of a Newton-type map at a root |
| `SplitNewton.lean` | `limit_is_root`, block/Schur structure, the `split_newton_*` theorems |

## Check the file

From this directory (`proofs/lean`):

```bash
lake build
```

No output and exit code 0 means every declaration typechecked. To check only the main module after edits (its imports must already be built):

```bash
lake env lean SplitNewton.lean
```

The Lean 4 editor extension only finds Mathlib if this folder is a workspace root (or listed in `SplitNewton.code-workspace`). Open `File → Open Workspace from File…` and choose that file, then open `SplitNewton.lean` and **Restart Server**. The status bar should show Lean `v4.37.0-rc1`, not `v4.36.0`. The Infoview spinner on first open is Lake loading Mathlib and can take a minute.
