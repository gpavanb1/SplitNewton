import Mathlib

/-!
# Derivative of a Newton-type map at a root

Reference-level calculus, independent of Split-Newton.
* `continuousAt_ringInverse` : operator inversion is continuous at an invertible operator.
* `hasFDerivAt_newtonMap`    : `y ↦ y - (M̂ y)⁻¹ g(y)` has derivative `I - M̂(ys)⁻¹ g'` at a root.
-/

open Filter Topology
open scoped ENNReal NNReal Matrix

namespace SplitNewton

/-! ## 5. Derivative of the outer map at the root -/

section Derivative

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]

/-- Inversion of operators is continuous at an invertible operator: if `f` is continuous at
`x` and `f x` is invertible then `y ↦ (f y)⁻¹` is continuous at `x`. -/
theorem continuousAt_ringInverse {X : Type*} [TopologicalSpace X] {f : X → E →L[ℝ] E} {x : X}
    (hf : ContinuousAt f x) (hu : IsUnit (f x)) :
    ContinuousAt (fun y => Ring.inverse (f y)) x := by
  have h := NormedRing.inverse_continuousAt hu.unit
  rw [hu.unit_spec] at h
  exact ContinuousAt.comp (g := Ring.inverse) h hf

/-- **Derivative of a (generalised) Newton map at a root.**  Let `g` be differentiable at `ys`
with `g ys = 0`, and let `Mhat` be *any* operator-valued map that is merely continuous at
`ys` with `Mhat ys` invertible (exact Jacobian, frozen Jacobian, `δ⁻¹ I + J`, quasi-Newton, ...).
Then `Φ y = y - (Mhat y)⁻¹ g y` has derivative `I - (Mhat ys)⁻¹ g'` at `ys`: the variation of
`Mhat` drops out because it multiplies `g y`, which vanishes at the root. -/
theorem hasFDerivAt_newtonMap {g : E → E} {g' : E →L[ℝ] E} {ys : E}
    (hg : HasFDerivAt g g' ys) (hg0 : g ys = 0)
    {Mhat : E → E →L[ℝ] E} (hM : ContinuousAt Mhat ys) (hu : IsUnit (Mhat ys)) :
    HasFDerivAt (fun y => y - Ring.inverse (Mhat y) (g y))
      (1 - Ring.inverse (Mhat ys) * g') ys := by
  set P : E → E →L[ℝ] E := fun y => Ring.inverse (Mhat y) with hPdef
  set P0 := P ys with hP0
  have hPc : ContinuousAt P ys := by
    have h := NormedRing.inverse_continuousAt hu.unit
    rw [hu.unit_spec] at h
    exact ContinuousAt.comp (g := Ring.inverse) h hM
  have hr : (fun y => g y - g ys - g' (y - ys)) =o[𝓝 ys] (fun y => y - ys) := hg.isLittleO
  have hgO : (fun y => g y) =O[𝓝 ys] (fun y => y - ys) := by
    simpa [hg0] using hg.isBigO_sub
  have hPd : (fun y => ‖P y - P0‖) =o[𝓝 ys] (fun _ => (1 : ℝ)) := by
    rw [Asymptotics.isLittleO_one_iff]
    exact tendsto_iff_norm_sub_tendsto_zero.1 hPc.tendsto
  have e1 : (fun y => ‖P0‖ * ‖g y - g ys - g' (y - ys)‖) =o[𝓝 ys] (fun y => ‖y - ys‖) :=
    hr.norm_norm.const_mul_left _
  have e2 : (fun y => ‖P y - P0‖ * ‖g y‖) =o[𝓝 ys] (fun y => ‖y - ys‖) := by
    have := hPd.mul_isBigO hgO.norm_norm
    simpa using this
  have hbig : (fun y => (y - Ring.inverse (Mhat y) (g y)) - (ys - Ring.inverse (Mhat ys) (g ys))
        - (1 - Ring.inverse (Mhat ys) * g') (y - ys))
      =O[𝓝 ys] (fun y => ‖P0‖ * ‖g y - g ys - g' (y - ys)‖ + ‖P y - P0‖ * ‖g y‖) := by
    refine Asymptotics.IsBigO.of_bound 1 (Filter.Eventually.of_forall fun y => ?_)
    have heq : (y - Ring.inverse (Mhat y) (g y)) - (ys - Ring.inverse (Mhat ys) (g ys))
        - (1 - Ring.inverse (Mhat ys) * g') (y - ys)
        = -(P0 (g y - g ys - g' (y - ys)) + (P y - P0) (g y)) := by
      simp only [hg0, hPdef, hP0, map_zero, sub_apply,
        mul_apply_eq_comp, one_apply_eq_self, map_sub]
      abel
    rw [heq, one_mul, Real.norm_of_nonneg (by positivity), norm_neg]
    calc ‖P0 (g y - g ys - g' (y - ys)) + (P y - P0) (g y)‖
        ≤ ‖P0 (g y - g ys - g' (y - ys))‖ + ‖(P y - P0) (g y)‖ := norm_add_le _ _
      _ ≤ ‖P0‖ * ‖g y - g ys - g' (y - ys)‖ + ‖P y - P0‖ * ‖g y‖ :=
        add_le_add (P0.le_opNorm _) ((P y - P0).le_opNorm _)
  have := hbig.trans_isLittleO (e1.add e2).of_norm_right
  rw [hasFDerivAt_iff_isLittleO]
  exact this

end Derivative

end SplitNewton
