import Mathlib

/-!
# Spectral radius of real operators on `ℝⁿ` (Gelfand)

Reference-level background, independent of Split-Newton.
* `Vec n`, `toMat`, `specRad` : spectral radius of a real operator (complex eigenvalues included).
* `exists_pow_norm_lt_one`    : `ρ(T) < 1 ⇒ ‖T^N‖ < 1` for some `N` (Gelfand's formula).
* `opNorm_le_matrix`          : `‖T^k‖ ≤ ‖(M_ℂ)^k‖`, used to compare with Mathlib's `spectralRadius`.
-/

open Filter Topology
open scoped ENNReal NNReal Matrix

namespace SplitNewton

/-! ## 2. Finite-dimensional setting and the spectral radius -/

/-- Coordinate space `ℝⁿ` (sup norm).  Every finite-dimensional real vector space is linearly
isomorphic to one of these, so working here loses no generality for finite systems; it is only
needed to talk about the *spectral radius* of a real operator (via its complexified matrix). -/
abbrev Vec (n : ℕ) := Fin n → ℝ

/-- Matrix of a continuous linear endomorphism of `ℝⁿ` in the standard basis. -/
noncomputable def toMat {n : ℕ} (T : Vec n →L[ℝ] Vec n) : Matrix (Fin n) (Fin n) ℝ :=
  LinearMap.toMatrix' (T : Vec n →ₗ[ℝ] Vec n)

/-- `toMat` is multiplicative (composition of operators = product of matrices). -/
theorem toMat_mul {n : ℕ} (S T : Vec n →L[ℝ] Vec n) : toMat (S * T) = toMat S * toMat T := by
  unfold toMat
  rw [← LinearMap.toMatrix'_mul]
  rfl

/-- `toMat` sends the identity operator to the identity matrix. -/
theorem toMat_one {n : ℕ} : toMat (1 : Vec n →L[ℝ] Vec n) = 1 := by
  unfold toMat
  rw [← LinearMap.toMatrix'_one]
  rfl

/-- `toMat` commutes with powers. -/
theorem toMat_pow {n : ℕ} (T : Vec n →L[ℝ] Vec n) (k : ℕ) : toMat (T ^ k) = toMat T ^ k := by
  induction k with
  | zero => simpa using toMat_one
  | succ k ih => rw [pow_succ, toMat_mul, ih, pow_succ]

/-- **Spectral radius** `ρ(T)` of a real operator `T` on `ℝⁿ`: the largest modulus of an
eigenvalue of `T`, *complex eigenvalues included*.  It is Mathlib's `spectralRadius ℂ` of the
entrywise complexification of the matrix of `T` (a real operator such as a rotation has no real
spectrum, so the complex spectrum is the correct notion). -/
noncomputable def specRad {n : ℕ} (T : Vec n →L[ℝ] Vec n) : ℝ≥0∞ :=
  spectralRadius ℂ ((toMat T).map (algebraMap ℝ ℂ))

/-! ## 3. Spectral radius < 1 ⇒ some power is a strict contraction (Gelfand's formula) -/

section
-- `L∞` operator norm on matrices: a complete submultiplicative norm, which is all that
-- Gelfand's formula needs (the spectral radius itself does not depend on the norm chosen).
attribute [local instance] Matrix.linftyOpNormedRing Matrix.linftyOpNormedAlgebra

/-- entrywise complexification preserves the `L∞` operator norm. -/
theorem nnnorm_map_ofReal {n : ℕ} (R : Matrix (Fin n) (Fin n) ℝ) :
    ‖R.map (algebraMap ℝ ℂ)‖₊ = ‖R‖₊ := by
  rw [Matrix.linfty_opNNNorm_def, Matrix.linfty_opNNNorm_def]
  simp [Matrix.map_apply]

/-- Complexification commutes with matrix powers. -/
theorem map_pow_ofReal {n : ℕ} (R : Matrix (Fin n) (Fin n) ℝ) (k : ℕ) :
    (R.map (algebraMap ℝ ℂ)) ^ k = (R ^ k).map (algebraMap ℝ ℂ) := by
  have := map_pow ((algebraMap ℝ ℂ).mapMatrix : Matrix (Fin n) (Fin n) ℝ →+* _) R k
  simpa using this.symm

/-- Real-valued version of `nnnorm_map_ofReal`. -/
theorem norm_map_ofReal {n : ℕ} (R : Matrix (Fin n) (Fin n) ℝ) :
    ‖R.map (algebraMap ℝ ℂ)‖ = ‖R‖ := by
  have := congrArg ((↑) : ℝ≥0 → ℝ) (nnnorm_map_ofReal R)
  simpa using this

/-- The operator norm of `T^k` (sup norm on `ℝⁿ`) is bounded by the `L∞` norm of the
complexified matrix power. -/
theorem opNorm_le_matrix {n : ℕ} (T : Vec n →L[ℝ] Vec n) (k : ℕ) :
    ‖T ^ k‖ ≤ ‖((toMat T).map (algebraMap ℝ ℂ)) ^ k‖ := by
  rw [map_pow_ofReal, norm_map_ofReal]
  refine ContinuousLinearMap.opNorm_le_bound _ (norm_nonneg _) (fun x => ?_)
  have h1 : (T ^ k) x = ((toMat T) ^ k) *ᵥ x := by
    rw [← toMat_pow]
    unfold toMat
    rw [LinearMap.toMatrix'_mulVec]
    rfl
  rw [h1]
  exact Matrix.linfty_opNorm_mulVec _ _

/-- **Spectral radius ⇒ power contraction.**  If `ρ(T) < 1` then `‖T^N‖ < 1` for some `N ≥ 1`.
This is the only place where the spectral radius enters; the proof is Gelfand's formula
`‖Mⁿ‖^{1/n} → ρ(M)` applied to the complexified matrix `M`. -/
theorem exists_pow_norm_lt_one {n : ℕ} (T : Vec n →L[ℝ] Vec n) (h : specRad T < 1) :
    ∃ N : ℕ, 0 < N ∧ ‖T ^ N‖ < 1 := by
  set M := (toMat T).map (algebraMap ℝ ℂ) with hM
  have hG := spectrum.pow_nnNorm_pow_one_div_tendsto_nhds_spectralRadius M
  have hev := hG.eventually (gt_mem_nhds (show spectralRadius ℂ M < 1 from h))
  obtain ⟨N0, hN0⟩ := eventually_atTop.1 hev
  refine ⟨max N0 1, lt_max_of_lt_right one_pos, ?_⟩
  have h1 := hN0 (max N0 1) (le_max_left _ _)
  have hpos : (0 : ℝ) < 1 / ((max N0 1 : ℕ) : ℝ) := by
    have : (0 : ℝ) < ((max N0 1 : ℕ) : ℝ) := by
      exact_mod_cast lt_max_of_lt_right one_pos
    positivity
  have h2 : (‖M ^ (max N0 1)‖₊ : ℝ≥0∞) < 1 := by
    by_contra hcon
    exact absurd h1 (not_lt.2 (ENNReal.one_le_rpow (not_lt.1 hcon) hpos))
  have h3 : ‖M ^ (max N0 1)‖ < 1 := by
    have := ENNReal.coe_lt_one_iff.1 h2
    exact_mod_cast this
  exact lt_of_le_of_lt (opNorm_le_matrix T _) h3

end

end SplitNewton
