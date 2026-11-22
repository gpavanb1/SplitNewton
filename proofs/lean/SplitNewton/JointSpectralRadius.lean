import Mathlib
import SplitNewton.Spectral
import SplitNewton.Switched

/-!
# Joint spectral radius (Rota–Strang)

Reference-level background, independent of Split-Newton.
* `jointSpecRad`                  : `ρ̂(T, S) = limsup_N (sup ‖N-fold products‖)^{1/N}`.
* `exists_prod_norm_lt_one_of_jsr`: `ρ̂ < 1 ⇒` uniform `N`-step contraction (what
  `switched_local_convergence` consumes).
* `jointSpecRad_singleton_le`     : `ρ̂({T}) ≤ ρ(T)` (Gelfand), so `ρ̂` extends `ρ`.
-/

open Filter Topology
open scoped ENNReal NNReal Matrix

namespace SplitNewton

/-! ### Joint spectral radius

For a family `T_i` (`i ∈ S`) of operators the **joint spectral radius** is
`ρ̂ = limsup_N (sup over all N-fold products ‖T_{i_{N-1}} ⋯ T_{i_0}‖)^{1/N}`
(Rota–Strang).  It does not depend on the norm, and for a single operator it is the spectral
radius (`jointSpecRad_singleton_le`).  Below we prove that `ρ̂ < 1` gives a uniform `N`-step
contraction, which is what the convergence engine consumes. -/

section JSR

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] {ι : Type*}

/-- The largest norm of an `N`-fold product `T_{i_{N-1}} ⋯ T_{i_0}` over all switching
patterns `i` with values in `S`. -/
noncomputable def maxProdNorm (T : ι → E →L[ℝ] E) (S : Set ι) (N : ℕ) : ℝ≥0∞ :=
  ⨆ (i : ℕ → ι) (_ : ∀ k, i k ∈ S), (‖prodFrom T i N‖₊ : ℝ≥0∞)

/-- **Joint spectral radius** of the family `{T_i : i ∈ S}`. -/
noncomputable def jointSpecRad (T : ι → E →L[ℝ] E) (S : Set ι) : ℝ≥0∞ :=
  limsup (fun N : ℕ => maxProdNorm T S N ^ (1 / N : ℝ)) atTop

/-- **`ρ̂ < 1` ⇒ uniform `N`-step contraction.**  If the joint spectral radius is `< 1`, then for
some `N ≥ 1` every product of `N` operators from the family (any switching pattern in `S`) has
norm `≤ q < 1`. -/
theorem exists_prod_norm_lt_one_of_jsr {T : ι → E →L[ℝ] E} {S : Set ι}
    (h : jointSpecRad T S < 1) :
    ∃ N : ℕ, 0 < N ∧ ∃ q : ℝ, q < 1 ∧
      ∀ i : ℕ → ι, (∀ k, i k ∈ S) → ‖prodFrom T i N‖ ≤ q := by
  have hev := eventually_lt_of_limsup_lt h
  obtain ⟨N0, hN0⟩ := eventually_atTop.1 hev
  have hNpos : 0 < max N0 1 := lt_max_of_lt_right one_pos
  have h1 := hN0 (max N0 1) (le_max_left _ _)
  have hpos : (0 : ℝ) < 1 / ((max N0 1 : ℕ) : ℝ) := by
    have : (0 : ℝ) < ((max N0 1 : ℕ) : ℝ) := by exact_mod_cast hNpos
    positivity
  -- the N-th root is < 1, hence the N-th maximal product norm is < 1
  have h2 : maxProdNorm T S (max N0 1) < 1 := by
    by_contra hcon
    exact absurd h1 (not_lt.2 (ENNReal.one_le_rpow (not_lt.1 hcon) hpos))
  have hfin : maxProdNorm T S (max N0 1) ≠ ⊤ := ne_top_of_lt h2
  refine ⟨max N0 1, hNpos, (maxProdNorm T S (max N0 1)).toReal, ?_, fun i hi => ?_⟩
  · have := (ENNReal.toReal_lt_toReal hfin ENNReal.one_ne_top).2 h2
    simpa using this
  · have hle : (‖prodFrom T i (max N0 1)‖₊ : ℝ≥0∞) ≤ maxProdNorm T S (max N0 1) :=
      le_iSup₂ (f := fun (i : ℕ → ι) (_ : ∀ k, i k ∈ S) => (‖prodFrom T i (max N0 1)‖₊ : ℝ≥0∞))
        i hi
    have := ENNReal.toReal_mono hfin hle
    simpa using this

end JSR

section JSRSpectral

attribute [local instance] Matrix.linftyOpNormedRing Matrix.linftyOpNormedAlgebra

/-- **The joint spectral radius extends the spectral radius**: for a single operator
`ρ̂({T}) ≤ ρ(T)` (Gelfand's formula).  So `ρ(T) < 1` implies `ρ̂({T}) < 1`, and the variable-pseudo-time
theorem with a one-point set of pseudo-times recovers the fixed-pseudo-time theorem. -/
theorem jointSpecRad_singleton_le {n : ℕ} (T0 : Vec n →L[ℝ] Vec n) :
    jointSpecRad (fun _ : Unit => T0) Set.univ ≤ specRad T0 := by
  set M := (toMat T0).map (algebraMap ℝ ℂ) with hM
  have hG := spectrum.pow_nnNorm_pow_one_div_tendsto_nhds_spectralRadius M
  have hmax : ∀ N : ℕ, maxProdNorm (fun _ : Unit => T0) Set.univ N = (‖T0 ^ N‖₊ : ℝ≥0∞) := by
    intro N
    unfold maxProdNorm
    simp [prodFrom_const]
  unfold jointSpecRad
  calc limsup (fun N : ℕ => maxProdNorm (fun _ : Unit => T0) Set.univ N ^ (1 / N : ℝ)) atTop
      ≤ limsup (fun N : ℕ => (‖M ^ N‖₊ : ℝ≥0∞) ^ (1 / N : ℝ)) atTop := by
        refine limsup_le_limsup (Eventually.of_forall fun N => ?_)
        rw [hmax N]
        refine ENNReal.rpow_le_rpow ?_ (by positivity)
        have := opNorm_le_matrix T0 N
        exact_mod_cast this
    _ = specRad T0 := hG.limsup_eq

end JSRSpectral

end SplitNewton
