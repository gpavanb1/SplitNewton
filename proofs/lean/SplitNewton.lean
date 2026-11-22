import Mathlib
import SplitNewton.Spectral
import SplitNewton.Switched
import SplitNewton.JointSpectralRadius
import SplitNewton.Derivative

/-!
# Split-Newton: local convergence for general nonlinear systems

Model: x = (y, z), F = (G, H) with G, H NONLINEAR (only differentiability at the root is used).
* inner solve: for the current y, z = ζ(y) is solved exactly (H (y, ζ y) = 0)
* outer step : one (optionally pseudo-transient) Newton step in y, using the
  PARTIAL Jacobian  A = ∂_y G  (or any operator continuous at the root):

      y ↦ y - (Ahat y + c I)⁻¹ · G (y, ζ y),     c = δ⁻¹  (c = 0 : plain Newton).

Blocks at the root:  A = ∂_y G,  B = ∂_z G,  C = ∂_y H,  D = ∂_z H.

Main points
1. `limit_is_root`        : if the iteration converges, the limit is a root of G
                            (partial correctness; any continuous G, any normed space,
                            any pseudo-time sequence bounded below).

Reference-level background lives in separate modules, imported here:
  `SplitNewton.Spectral`            : `specRad`, `exists_pow_norm_lt_one` (Gelfand's formula),
  `SplitNewton.Switched`            : `switched_local_convergence`, `local_convergence` (Ostrowski),
  `SplitNewton.JointSpectralRadius` : `jointSpecRad`, `exists_prod_norm_lt_one_of_jsr`,
                                      `jointSpecRad_singleton_le`,
  `SplitNewton.Derivative`          : `continuousAt_ringInverse`, `hasFDerivAt_newtonMap`.

2. `exists_pow_norm_lt_one`: ρ(T) < 1 ⇒ ‖T^N‖ < 1 for some N  (Gelfand's formula).
3. `switched_local_convergence`: local geometric convergence of ANY time-varying iteration
                            y_{k+1} = Φ_{i k}(y_k) whose linearisations have uniformly
                            contracting N-fold products (joint spectral radius < 1).
   `local_convergence`   : its single-map case (Ostrowski): ρ(Φ'(ys)) < 1 ⇒ local convergence.
4. `hasFDerivAt_newtonMap`: derivative of y ↦ y - (Mhat y)⁻¹ g(y) at a root is I - Mhat(ys)⁻¹ g'.
5. `inner_derivative`     : the exact inner solve has sensitivity ζ' = -D⁻¹ C, hence the
                            reduced Jacobian is the Schur complement  A - B D⁻¹ C.
6. `split_newton_pure`    : **ρ(A⁻¹ B D⁻¹ C) < 1  ⇒  local geometric convergence** of plain
                            split Newton (outer y_k → ys and inner z_k → zs), for nonlinear F.
   `split_newton_pseudo_transient`: fixed pseudo-time δ, with
                            ρ((A + δ⁻¹ I)⁻¹ (δ⁻¹ I + B D⁻¹ C)) < 1.
   `split_newton_variable_pseudo_time`: pseudo-time sequence δ_k VARYING at every step (inverse
                            pseudo-times in a compact set Δ): converges when the family
                            {(A + c I)⁻¹ (c I + B D⁻¹ C) : c ∈ Δ} has JOINT SPECTRAL RADIUS < 1,
                            uniformly in the sequence.  (`jointSpecRad` is defined in `SplitNewton.JointSpectralRadius`;
                            `exists_prod_norm_lt_one_of_jsr` turns ρ̂ < 1 into a uniform N-step
                            contraction; `jointSpecRad_singleton_le` : ρ̂({T}) ≤ ρ(T).)
   `split_newton_variable_pseudo_time_of_norm`: the special case N = 1 (a common contraction).

Scope: the convergence theorems are LOCAL (start near the root) and finite-dimensional.

STATUS: all declarations are checked with `lake build` (Mathlib, Lean v4.37.0-rc1).
-/

open Filter Topology
open scoped ENNReal NNReal Matrix

namespace SplitNewton

/-! ## 1. Partial correctness: a convergent iteration converges to a root -/

theorem limit_is_root {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {Z : Type*} [TopologicalSpace Z]
    (G : E × Z → E) (hG : Continuous G)
    (y : ℕ → E) (z : ℕ → Z) (ys : E) (zs : Z)
    (hx : Tendsto (fun k => (y k, z k)) atTop (𝓝 (ys, zs)))
    (δ : ℕ → ℝ) (δmin M : ℝ) (hδmin : 0 < δmin) (hδ : ∀ k, δmin ≤ δ k)
    (A : ℕ → E →L[ℝ] E) (hA : ∀ k, ‖A k‖ ≤ M)
    -- the outer step:  (δ⁻¹ I + A) (y_{k+1} - y_k) = - G(y_k, z_k)
    (hstep : ∀ k, (δ k)⁻¹ • (y (k+1) - y k) + A k (y (k+1) - y k) = - G (y k, z k)) :
    G (ys, zs) = 0 := by
  have hy : Tendsto y atTop (𝓝 ys) := (continuous_fst.tendsto (ys, zs)).comp hx
  have hΔ : Tendsto (fun k => y (k+1) - y k) atTop (𝓝 0) := by
    have := (hy.comp (tendsto_add_atTop_nat 1)).sub hy
    simpa using this
  have hbound : ∀ k, ‖G (y k, z k)‖ ≤ (δmin⁻¹ + M) * ‖y (k+1) - y k‖ := by
    intro k
    have h1 : ‖G (y k, z k)‖
        = ‖(δ k)⁻¹ • (y (k+1) - y k) + A k (y (k+1) - y k)‖ := by
      rw [hstep k, norm_neg]
    have hpos : 0 ≤ (δ k)⁻¹ := inv_nonneg.2 (le_trans hδmin.le (hδ k))
    have h2 : ‖(δ k)⁻¹ • (y (k+1) - y k)‖ ≤ δmin⁻¹ * ‖y (k+1) - y k‖ := by
      rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg hpos]
      have _h : δmin ≤ δ k := hδ k
      have : (δ k)⁻¹ ≤ δmin⁻¹ := by gcongr
      exact mul_le_mul_of_nonneg_right this (norm_nonneg _)
    have h3 : ‖A k (y (k+1) - y k)‖ ≤ M * ‖y (k+1) - y k‖ :=
      le_trans ((A k).le_opNorm _) (mul_le_mul_of_nonneg_right (hA k) (norm_nonneg _))
    calc ‖G (y k, z k)‖
        = ‖(δ k)⁻¹ • (y (k+1) - y k) + A k (y (k+1) - y k)‖ := h1
      _ ≤ ‖(δ k)⁻¹ • (y (k+1) - y k)‖ + ‖A k (y (k+1) - y k)‖ := norm_add_le _ _
      _ ≤ δmin⁻¹ * ‖y (k+1) - y k‖ + M * ‖y (k+1) - y k‖ := add_le_add h2 h3
      _ = (δmin⁻¹ + M) * ‖y (k+1) - y k‖ := by ring
  have hG0 : Tendsto (fun k => G (y k, z k)) atTop (𝓝 0) := by
    refine squeeze_zero_norm hbound ?_
    simpa using hΔ.norm.const_mul (δmin⁻¹ + M)
  have hG1 : Tendsto (fun k => G (y k, z k)) atTop (𝓝 (G (ys, zs))) :=
    (hG.tendsto (ys, zs)).comp hx
  exact tendsto_nhds_unique hG1 hG0

/-! ## 2. Block structure `x = (y, z)`, `F = (G, H)` and the split-Newton theorems -/

section Blocks

variable {n m : ℕ}

/-- `∂_y` block of a derivative acting on the product `Vec n × Vec m`. -/
noncomputable def blockFst {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    (L : Vec n × Vec m →L[ℝ] F) : Vec n →L[ℝ] F :=
  L.comp (ContinuousLinearMap.inl ℝ (Vec n) (Vec m))

/-- `∂_z` block of a derivative acting on the product `Vec n × Vec m`. -/
noncomputable def blockSnd {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    (L : Vec n × Vec m →L[ℝ] F) : Vec m →L[ℝ] F :=
  L.comp (ContinuousLinearMap.inr ℝ (Vec n) (Vec m))

theorem block_apply {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    (L : Vec n × Vec m →L[ℝ] F) (u : Vec n) (v : Vec m) :
    L (u, v) = blockFst L u + blockSnd L v := by
  unfold blockFst blockSnd
  rw [ContinuousLinearMap.comp_apply, ContinuousLinearMap.comp_apply, ← map_add]
  simp

theorem reduced_hasFDerivAt {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    {Φ : Vec n × Vec m → F} {L : Vec n × Vec m →L[ℝ] F} {ys : Vec n} {zs : Vec m}
    (hΦ : HasFDerivAt Φ L (ys, zs)) {ζ : Vec n → Vec m} {ζ' : Vec n →L[ℝ] Vec m}
    (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys) :
    HasFDerivAt (fun y => Φ (y, ζ y)) (blockFst L + (blockSnd L).comp ζ') ys := by
  have h1 : HasFDerivAt (fun y => (y, ζ y)) ((ContinuousLinearMap.id ℝ (Vec n)).prod ζ') ys :=
    (hasFDerivAt_id ys).prodMk hζ
  have hΦ' : HasFDerivAt Φ L (ys, ζ ys) := by rwa [hζ0]
  have h2 := hΦ'.comp ys h1
  refine h2.congr_fderiv ?_
  ext u : 1
  simp [block_apply]

end Blocks


section Main

variable {n m : ℕ}

/-- **Sensitivity of the exact inner solve.**  If `ζ y` solves `H (y, ζ y) = 0` near `ys`, `ζ` is
differentiable at `ys`, and `H` is differentiable at the root with `D = ∂_z H` invertible, then
differentiating `H (y, ζ y) = 0` gives `C + D ζ' = 0`, i.e. `ζ' = -D⁻¹ C` (with `C = ∂_y H`).
(The implicit function theorem is what supplies such a `ζ` in practice; we take it as given
because it is exactly the object the algorithm computes in its inner solve.) -/
theorem inner_derivative
    {H : Vec n × Vec m → Vec m} {H' : Vec n × Vec m →L[ℝ] Vec m} {ys : Vec n} {zs : Vec m}
    (hH : HasFDerivAt H H' (ys, zs)) {ζ : Vec n → Vec m} {ζ' : Vec n →L[ℝ] Vec m}
    (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys) (hsol : ∀ᶠ y in 𝓝 ys, H (y, ζ y) = 0)
    (hD : IsUnit (blockSnd H')) :
    ζ' = -((Ring.inverse (blockSnd H')).comp (blockFst H')) := by
  have h1 := reduced_hasFDerivAt hH hζ0 hζ
  have h2 : HasFDerivAt (fun y => H (y, ζ y)) (0 : Vec n →L[ℝ] Vec m) ys :=
    (hasFDerivAt_const (0 : Vec m) ys).congr_of_eventuallyEq hsol
  have h3 : blockFst H' + (blockSnd H').comp ζ' = 0 := h1.unique h2
  have hinv : Ring.inverse (blockSnd H') * blockSnd H' = 1 := Ring.inverse_mul_cancel _ hD
  have h4 : (Ring.inverse (blockSnd H')).comp (blockFst H' + (blockSnd H').comp ζ') = 0 := by
    rw [h3]; simp
  rw [ContinuousLinearMap.comp_add, ← ContinuousLinearMap.comp_assoc] at h4
  have h5 : (Ring.inverse (blockSnd H')).comp (blockSnd H') = 1 := hinv
  have h6 : (1 : Vec m →L[ℝ] Vec m).comp ζ' = ζ' := ContinuousLinearMap.id_comp ζ'
  rw [h5, h6] at h4
  exact eq_neg_of_add_eq_zero_right h4

/-- The **reduced (Schur-complement) Jacobian** `A - B D⁻¹ C` of the residual
`y ↦ G (y, ζ y)`, where `A, B = ∂_y G, ∂_z G` and `C, D = ∂_y H, ∂_z H`. -/
noncomputable def reducedJac (G' : Vec n × Vec m →L[ℝ] Vec n) (H' : Vec n × Vec m →L[ℝ] Vec m) :
    Vec n →L[ℝ] Vec n :=
  blockFst G' - (blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H'))

/-- **One outer cycle of split Newton.**  The inner variable is solved exactly (`ζ y`), then one
outer step `y ↦ y - (Mhat y)⁻¹ G (y, ζ y)` is taken, where `Mhat y` is the outer step operator:
the partial Jacobian `∂_y G` for plain Newton, `δ⁻¹ I + ∂_y G` for pseudo-transient steps, etc. -/
noncomputable def splitMap (G : Vec n × Vec m → Vec n) (ζ : Vec n → Vec m)
    (Mhat : Vec n → (Vec n →L[ℝ] Vec n)) (y : Vec n) : Vec n :=
  y - Ring.inverse (Mhat y) (G (y, ζ y))

/-- **Local convergence of split Newton (general nonlinear, finite-dimensional).**

Setting: `F = (G, H)` is nonlinear, differentiable at a root `(ys, zs)`; the inner solve
`ζ` satisfies `H (y, ζ y) = 0` near `ys` and is differentiable at `ys`; `D = ∂_z H` is
invertible; the outer step operator `Mhat` is continuous at `ys` and invertible there.
If the linearised outer iteration matrix `I - Mhat(ys)⁻¹ (A - B D⁻¹ C)` has spectral radius
`< 1`, then every start close enough to `ys` converges geometrically to `ys`. -/
theorem split_newton_converges
    (G : Vec n × Vec m → Vec n) (H : Vec n × Vec m → Vec m) (ys : Vec n) (zs : Vec m)
    (G' : Vec n × Vec m →L[ℝ] Vec n) (H' : Vec n × Vec m →L[ℝ] Vec m)
    (hG : HasFDerivAt G G' (ys, zs)) (hH : HasFDerivAt H H' (ys, zs))
    (hG0 : G (ys, zs) = 0)
    (ζ : Vec n → Vec m) (ζ' : Vec n →L[ℝ] Vec m) (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys)
    (hsol : ∀ᶠ y in 𝓝 ys, H (y, ζ y) = 0) (hD : IsUnit (blockSnd H'))
    (Mhat : Vec n → (Vec n →L[ℝ] Vec n)) (hM : ContinuousAt Mhat ys) (hMu : IsUnit (Mhat ys))
    (hρ : specRad (1 - Ring.inverse (Mhat ys) * reducedJac G' H') < 1) :
    ∃ ε > 0, ∃ C κ : ℝ, 0 ≤ κ ∧ κ < 1 ∧ ∀ y0 : Vec n, ‖y0 - ys‖ < ε →
      ∀ k : ℕ, ‖(splitMap G ζ Mhat)^[k] y0 - ys‖ ≤ C * κ ^ k * ‖y0 - ys‖ := by
  have hred := reduced_hasFDerivAt hG hζ0 hζ
  have hζ' := inner_derivative hH hζ0 hζ hsol hD
  have hg' : blockFst G' + (blockSnd G').comp ζ' = reducedJac G' H' := by
    rw [hζ']; unfold reducedJac; simp [sub_eq_add_neg]
  rw [hg'] at hred
  have hg0 : G (ys, ζ ys) = 0 := by rw [hζ0]; exact hG0
  have hΦ := hasFDerivAt_newtonMap hred hg0 hM hMu
  have hfix : splitMap G ζ Mhat ys = ys := by
    unfold splitMap; rw [hg0]; simp
  obtain ⟨N, hN, hT⟩ := exists_pow_norm_lt_one _ hρ
  exact local_convergence hfix hΦ hN hT


/-- A geometric bound `‖uₖ - ys‖ ≤ C κᵏ r` with `0 ≤ κ < 1` forces `uₖ → ys`. -/
theorem tendsto_of_geom {E : Type*} [NormedAddCommGroup E] {u : ℕ → E} {ys : E} {C κ r : ℝ}
    (hκ0 : 0 ≤ κ) (hκ1 : κ < 1) (hb : ∀ k, ‖u k - ys‖ ≤ C * κ ^ k * r) :
    Tendsto u atTop (𝓝 ys) := by
  rw [tendsto_iff_norm_sub_tendsto_zero]
  have h0 : Tendsto (fun k : ℕ => C * κ ^ k * r) atTop (𝓝 0) := by
    simpa using ((tendsto_pow_atTop_nhds_zero_of_lt_one hκ0 hκ1).const_mul C).mul_const r
  exact squeeze_zero (fun k => norm_nonneg _) hb h0

/-- **Plain Newton on the outer variables (library default `dt0 = 0`).**  Here the outer
operator is the partial Jacobian `Ahat y ≈ ∂_y G`, continuous at `ys` with `Ahat ys = A`, and the
linearised iteration matrix collapses to `A⁻¹ B D⁻¹ C`.  Hence:

  **`ρ(A⁻¹ B D⁻¹ C) < 1`  ⇒  the outer iterates `yₖ → ys` and the inner iterates `zₖ → zs`.**

This is the coupling condition; it is *not* implied by nonsingularity of the diagonal blocks
alone, since `A⁻¹ B D⁻¹ C` can have spectral radius `≥ 1` even when `A` and `D` are invertible. -/
theorem split_newton_pure
    (G : Vec n × Vec m → Vec n) (H : Vec n × Vec m → Vec m) (ys : Vec n) (zs : Vec m)
    (G' : Vec n × Vec m →L[ℝ] Vec n) (H' : Vec n × Vec m →L[ℝ] Vec m)
    (hG : HasFDerivAt G G' (ys, zs)) (hH : HasFDerivAt H H' (ys, zs))
    (hG0 : G (ys, zs) = 0)
    (ζ : Vec n → Vec m) (ζ' : Vec n →L[ℝ] Vec m) (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys)
    (hsol : ∀ᶠ y in 𝓝 ys, H (y, ζ y) = 0) (hD : IsUnit (blockSnd H'))
    (Ahat : Vec n → (Vec n →L[ℝ] Vec n)) (hA : ContinuousAt Ahat ys)
    (hAys : Ahat ys = blockFst G') (hAu : IsUnit (blockFst G'))
    (hρ : specRad (Ring.inverse (blockFst G') *
        ((blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H')))) < 1) :
    ∃ ε > 0, ∀ y0 : Vec n, ‖y0 - ys‖ < ε →
      Tendsto (fun k => (splitMap G ζ Ahat)^[k] y0) atTop (𝓝 ys) ∧
      Tendsto (fun k => ζ ((splitMap G ζ Ahat)^[k] y0)) atTop (𝓝 zs) := by
  have hinv : Ring.inverse (blockFst G') * blockFst G' = 1 := Ring.inverse_mul_cancel _ hAu
  have hT : 1 - Ring.inverse (Ahat ys) * reducedJac G' H' =
      Ring.inverse (blockFst G') *
        ((blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H'))) := by
    rw [hAys]; unfold reducedJac
    rw [mul_sub, hinv]; abel
  obtain ⟨ε, hε, C, κ, hκ0, hκ1, hb⟩ := split_newton_converges G H ys zs G' H' hG hH hG0 ζ ζ' hζ0
    hζ hsol hD Ahat hA (hAys ▸ hAu) (by rw [hT]; exact hρ)
  refine ⟨ε, hε, fun y0 hy0 => ?_⟩
  have h1 : Tendsto (fun k => (splitMap G ζ Ahat)^[k] y0) atTop (𝓝 ys) :=
    tendsto_of_geom hκ0 hκ1 (fun k => hb y0 hy0 k)
  refine ⟨h1, ?_⟩
  have := (hζ.continuousAt.tendsto).comp h1
  rwa [hζ0] at this


/-- **Pseudo-transient outer step with a fixed pseudo-time `δ`.**  The outer operator is
`δ⁻¹ I + Ahat y`; the linearised iteration matrix is `(A + δ⁻¹ I)⁻¹ (δ⁻¹ I + B D⁻¹ C)`.
(`δ → ∞` recovers `split_newton_pure`.) -/
theorem split_newton_pseudo_transient
    (G : Vec n × Vec m → Vec n) (H : Vec n × Vec m → Vec m) (ys : Vec n) (zs : Vec m)
    (G' : Vec n × Vec m →L[ℝ] Vec n) (H' : Vec n × Vec m →L[ℝ] Vec m)
    (hG : HasFDerivAt G G' (ys, zs)) (hH : HasFDerivAt H H' (ys, zs))
    (hG0 : G (ys, zs) = 0)
    (ζ : Vec n → Vec m) (ζ' : Vec n →L[ℝ] Vec m) (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys)
    (hsol : ∀ᶠ y in 𝓝 ys, H (y, ζ y) = 0) (hD : IsUnit (blockSnd H'))
    (δ : ℝ) (Ahat : Vec n → (Vec n →L[ℝ] Vec n)) (hA : ContinuousAt Ahat ys)
    (hAys : Ahat ys = blockFst G') (hMu : IsUnit (blockFst G' + δ⁻¹ • (1 : Vec n →L[ℝ] Vec n)))
    (hρ : specRad (Ring.inverse (blockFst G' + δ⁻¹ • (1 : Vec n →L[ℝ] Vec n)) *
        (δ⁻¹ • (1 : Vec n →L[ℝ] Vec n) +
          (blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H')))) < 1) :
    ∃ ε > 0, ∀ y0 : Vec n, ‖y0 - ys‖ < ε →
      Tendsto (fun k => (splitMap G ζ (fun y => Ahat y + δ⁻¹ • (1 : Vec n →L[ℝ] Vec n)))^[k] y0)
        atTop (𝓝 ys) := by
  set M0 := blockFst G' + δ⁻¹ • (1 : Vec n →L[ℝ] Vec n) with hM0
  set X := (blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H')) with hX
  have hinv : Ring.inverse M0 * M0 = 1 := Ring.inverse_mul_cancel _ hMu
  have hJ : reducedJac G' H' = M0 - (δ⁻¹ • (1 : Vec n →L[ℝ] Vec n) + X) := by
    unfold reducedJac; rw [hM0, hX]; abel
  have hT : 1 - Ring.inverse ((fun y => Ahat y + δ⁻¹ • (1 : Vec n →L[ℝ] Vec n)) ys) *
        reducedJac G' H' = Ring.inverse M0 * (δ⁻¹ • (1 : Vec n →L[ℝ] Vec n) + X) := by
    simp only [hAys]
    rw [hJ, mul_sub, ← hM0, hinv]; abel
  have hMc : ContinuousAt (fun y => Ahat y + δ⁻¹ • (1 : Vec n →L[ℝ] Vec n)) ys :=
    hA.add continuousAt_const
  obtain ⟨ε, hε, C, κ, hκ0, hκ1, hb⟩ := split_newton_converges G H ys zs G' H' hG hH hG0 ζ ζ' hζ0
    hζ hsol hD _ hMc (by simpa [hAys] using hMu) (by rw [hT]; exact hρ)
  exact ⟨ε, hε, fun y0 hy0 => tendsto_of_geom hκ0 hκ1 (fun k => hb y0 hy0 k)⟩


/-- **Pseudo-transient outer step with inverse pseudo-time `c = δ⁻¹`**:
`y ↦ y - (Ahat y + c I)⁻¹ G (y, ζ y)`.  (`c = 0` is plain Newton, `δ = ∞`.) -/
noncomputable def ptcStep (G : Vec n × Vec m → Vec n) (ζ : Vec n → Vec m)
    (Ahat : Vec n → (Vec n →L[ℝ] Vec n)) (c : ℝ) (y : Vec n) : Vec n :=
  y - Ring.inverse (Ahat y + c • (1 : Vec n →L[ℝ] Vec n)) (G (y, ζ y))

/-- Linearised outer iteration matrix `(A + c I)⁻¹ (c I + X)` for inverse pseudo-time `c`,
where `X = B D⁻¹ C` is the coupling term. -/
noncomputable def ptcMatrix (A X : Vec n →L[ℝ] Vec n) (c : ℝ) : Vec n →L[ℝ] Vec n :=
  Ring.inverse (A + c • (1 : Vec n →L[ℝ] Vec n)) * (c • (1 : Vec n →L[ℝ] Vec n) + X)

/-- **Split Newton with a VARIABLE pseudo-time sequence, `N`-step product form.**

Step `k` uses inverse pseudo-time `c k = (δ k)⁻¹ ∈ Δ`, where `Δ ⊂ ℝ` is a compact set of
admissible values (e.g. `Δ = [0, 1/δmin]`; `c k = 0` is a pure Newton step) and `A + c I` is
invertible for every `c ∈ Δ`.  Hypothesis: for some `N`, every product of `N` matrices
`(A + c I)⁻¹ (c I + B D⁻¹ C)` with `c₀, …, c_{N-1} ∈ Δ` has norm `≤ q < 1`.
(`split_newton_variable_pseudo_time` below derives this from `ρ̂ < 1`.) -/
theorem split_newton_variable_pseudo_time_of_products
    (G : Vec n × Vec m → Vec n) (H : Vec n × Vec m → Vec m) (ys : Vec n) (zs : Vec m)
    (G' : Vec n × Vec m →L[ℝ] Vec n) (H' : Vec n × Vec m →L[ℝ] Vec m)
    (hG : HasFDerivAt G G' (ys, zs)) (hH : HasFDerivAt H H' (ys, zs))
    (hG0 : G (ys, zs) = 0)
    (ζ : Vec n → Vec m) (ζ' : Vec n →L[ℝ] Vec m) (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys)
    (hsol : ∀ᶠ y in 𝓝 ys, H (y, ζ y) = 0) (hD : IsUnit (blockSnd H'))
    (Δ : Set ℝ) (hΔ : IsCompact Δ)
    (Ahat : Vec n → (Vec n →L[ℝ] Vec n)) (hA : ContinuousAt Ahat ys)
    (hAys : Ahat ys = blockFst G')
    (hres : ∀ c ∈ Δ, IsUnit (blockFst G' + c • (1 : Vec n →L[ℝ] Vec n)))
    (hjsr : ∃ N : ℕ, 0 < N ∧ ∃ q : ℝ, q < 1 ∧ ∀ c : ℕ → ℝ, (∀ k, c k ∈ Δ) →
      ‖prodFrom (ptcMatrix (blockFst G')
        ((blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H')))) c N‖ ≤ q) :
    ∃ ε > 0, ∃ C κ : ℝ, 0 ≤ κ ∧ κ < 1 ∧ ∀ c : ℕ → ℝ, (∀ k, c k ∈ Δ) →
      ∀ y0 : Vec n, ‖y0 - ys‖ < ε → ∀ k : ℕ,
        ‖switchedOrbit (ptcStep G ζ Ahat) c y0 k - ys‖ ≤ C * κ ^ k * ‖y0 - ys‖ := by
  -- notation
  set A := blockFst G' with hAdef
  set X := (blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H')) with hXdef
  set S0 := reducedJac G' H' with hS0def
  set one : Vec n →L[ℝ] Vec n := 1 with hone
  -- the reduced residual g y = G (y, ζ y) has derivative S0 at ys and vanishes there
  have hred := reduced_hasFDerivAt hG hζ0 hζ
  have hζ' := inner_derivative hH hζ0 hζ hsol hD
  have hg' : blockFst G' + (blockSnd G').comp ζ' = S0 := by
    rw [hζ', hS0def]; unfold reducedJac; simp [sub_eq_add_neg]
  rw [hg'] at hred
  have hg0 : G (ys, ζ ys) = 0 := by rw [hζ0]; exact hG0
  -- algebra: for `c ∈ Δ`,  ptcMatrix = 1 - (A + c)⁻¹ S0
  have hTeq : ∀ c ∈ Δ, ptcMatrix A X c = 1 - Ring.inverse (A + c • one) * S0 := by
    intro c hc
    have hinv : Ring.inverse (A + c • one) * (A + c • one) = 1 :=
      Ring.inverse_mul_cancel _ (hres c hc)
    have hS : S0 = (A + c • one) - (c • one + X) := by
      rw [hS0def]; unfold reducedJac; rw [hAdef, hXdef]; abel
    rw [hS, mul_sub, hinv]; unfold ptcMatrix; abel
  -- uniform bound on the resolvents `(A + c)⁻¹`, `c ∈ Δ` (continuity on a compact set)
  have hcont : ContinuousOn (fun c : ℝ => Ring.inverse (A + c • one)) Δ := by
    intro c hc
    refine ContinuousAt.continuousWithinAt ?_
    exact continuousAt_ringInverse (f := fun c : ℝ => A + c • one)
      (continuous_const.add (continuous_id.smul continuous_const)).continuousAt (hres c hc)
  obtain ⟨K, hK⟩ := hΔ.exists_bound_of_continuousOn hcont
  have hK0 : ∀ c ∈ Δ, ‖Ring.inverse (A + c • one)‖ ≤ max K 0 :=
    fun c hc => le_trans (hK c hc) (le_max_left _ _)
  set K' : ℝ := max K 0 with hK'
  have hK'0 : 0 ≤ K' := le_max_right _ _
  -- uniform (in `c ∈ Δ`) continuity of the resolvent in `y` at `ys` (tube lemma)
  have htube : ∀ η : ℝ, 0 < η → ∀ᶠ y in 𝓝 ys, ∀ c ∈ Δ,
      ‖Ring.inverse (Ahat y + c • one) - Ring.inverse (A + c • one)‖ ≤ η := by
    intro η hη
    refine hΔ.eventually_forall_of_forall_eventually
      (P := fun y c => ‖Ring.inverse (Ahat y + c • one) - Ring.inverse (A + c • one)‖ ≤ η)
      (fun c hc => ?_)
    have hf1 : ContinuousAt (fun z : Vec n × ℝ => Ahat z.1 + z.2 • one) (ys, c) :=
      (hA.comp continuousAt_fst).add (continuousAt_snd.smul continuousAt_const)
    have hf2 : ContinuousAt (fun z : Vec n × ℝ => A + z.2 • one) (ys, c) :=
      continuousAt_const.add (continuousAt_snd.smul continuousAt_const)
    have hu : IsUnit (A + c • one) := hres c hc
    have hc1 := continuousAt_ringInverse hf1 (by simpa [hAys] using hu)
    have hc2 := continuousAt_ringInverse hf2 (by simpa using hu)
    have hn := (hc1.sub hc2).norm
    have hval : ‖Ring.inverse (Ahat ys + c • one) - Ring.inverse (A + c • one)‖ = 0 := by
      rw [hAys]; simp
    exact (hn.eventually (gt_mem_nhds (show _ < η by simpa [hAys] using hη))).mono
      (fun z hz => hz.le)
  -- local bound of the residual and its linearisation remainder
  have hlo := hred.isLittleO
  have hgb : ∀ᶠ y in 𝓝 ys, ‖G (y, ζ y)‖ ≤ (‖S0‖ + 1) * ‖y - ys‖ := by
    have h1 := hlo.def (show (0 : ℝ) < 1 by norm_num)
    filter_upwards [h1] with y hy
    rw [hg0] at hy
    calc ‖G (y, ζ y)‖ = ‖(G (y, ζ y) - 0 - S0 (y - ys)) + S0 (y - ys)‖ := by simp
      _ ≤ ‖G (y, ζ y) - 0 - S0 (y - ys)‖ + ‖S0 (y - ys)‖ := norm_add_le _ _
      _ ≤ 1 * ‖y - ys‖ + ‖S0‖ * ‖y - ys‖ := add_le_add hy (S0.le_opNorm _)
      _ = (‖S0‖ + 1) * ‖y - ys‖ := by ring
  -- uniform linearisation estimate for the family of pseudo-transient steps
  have hlin : ∀ η : ℝ, 0 < η → ∃ ε : ℝ, 0 < ε ∧ ∀ c ∈ Δ, ∀ y : Vec n, ‖y - ys‖ < ε →
      ‖ptcStep G ζ Ahat c y - ys - ptcMatrix A X c (y - ys)‖ ≤ η * ‖y - ys‖ := by
    intro η hη
    set Lg : ℝ := ‖S0‖ + 1 with hLg
    have hLgpos : 0 < Lg := by positivity
    have hη1 : 0 < η / (2 * (K' + 1)) := by positivity
    have hη2 : 0 < η / (2 * Lg) := by positivity
    have e1 := hlo.def hη1
    have e3 := htube _ hη2
    obtain ⟨ε, hε, hball⟩ := Metric.eventually_nhds_iff_ball.1 (e1.and (hgb.and e3))
    refine ⟨ε, hε, fun c hc y hy => ?_⟩
    obtain ⟨hr, hg, ht⟩ := hball y (by simpa [dist_eq_norm] using hy)
    rw [hg0] at hr
    set P0 := Ring.inverse (A + c • one) with hP0
    set P := Ring.inverse (Ahat y + c • one) with hP
    have hident : ptcStep G ζ Ahat c y - ys - ptcMatrix A X c (y - ys)
        = -(P0 (G (y, ζ y) - 0 - S0 (y - ys)) + (P - P0) (G (y, ζ y))) := by
      rw [hTeq c hc]
      simp only [ptcStep, hP, hP0, sub_apply, mul_apply_eq_comp, one_apply_eq_self, map_sub,
        sub_zero]
      abel
    rw [hident, norm_neg]
    have t1 : ‖P0 (G (y, ζ y) - 0 - S0 (y - ys))‖ ≤ η / 2 * ‖y - ys‖ := by
      calc ‖P0 (G (y, ζ y) - 0 - S0 (y - ys))‖
          ≤ ‖P0‖ * ‖G (y, ζ y) - 0 - S0 (y - ys)‖ := P0.le_opNorm _
        _ ≤ K' * (η / (2 * (K' + 1)) * ‖y - ys‖) :=
            mul_le_mul (hK0 c hc) hr (norm_nonneg _) hK'0
        _ ≤ η / 2 * ‖y - ys‖ := by
            have : K' * (η / (2 * (K' + 1))) ≤ η / 2 := by
              rw [← mul_div_assoc, div_le_div_iff₀ (by positivity) (by norm_num)]
              nlinarith
            nlinarith [norm_nonneg (y - ys)]
    have t2 : ‖(P - P0) (G (y, ζ y))‖ ≤ η / 2 * ‖y - ys‖ := by
      calc ‖(P - P0) (G (y, ζ y))‖ ≤ ‖P - P0‖ * ‖G (y, ζ y)‖ := (P - P0).le_opNorm _
        _ ≤ η / (2 * Lg) * (Lg * ‖y - ys‖) :=
            mul_le_mul (ht c hc) hg (norm_nonneg _) hη2.le
        _ = η / 2 * ‖y - ys‖ := by field_simp
    calc ‖P0 (G (y, ζ y) - 0 - S0 (y - ys)) + (P - P0) (G (y, ζ y))‖
        ≤ ‖P0 (G (y, ζ y) - 0 - S0 (y - ys))‖ + ‖(P - P0) (G (y, ζ y))‖ := norm_add_le _ _
      _ ≤ η / 2 * ‖y - ys‖ + η / 2 * ‖y - ys‖ := add_le_add t1 t2
      _ = η * ‖y - ys‖ := by ring
  -- bounded linearisations
  have hKT : ∃ K : ℝ, ∀ c ∈ Δ, ‖ptcMatrix A X c‖ ≤ K := by
    refine ⟨1 + K' * ‖S0‖, fun c hc => ?_⟩
    rw [hTeq c hc]
    calc ‖(1 : Vec n →L[ℝ] Vec n) - Ring.inverse (A + c • one) * S0‖
        ≤ ‖(1 : Vec n →L[ℝ] Vec n)‖ + ‖Ring.inverse (A + c • one) * S0‖ := norm_sub_le _ _
      _ ≤ 1 + K' * ‖S0‖ := by
          refine add_le_add ?_ ?_
          · exact ContinuousLinearMap.norm_id_le
          · exact le_trans (norm_mul_le _ _) (mul_le_mul_of_nonneg_right (hK0 c hc) (norm_nonneg _))
  exact switched_local_convergence (Φ := ptcStep G ζ Ahat) (T := ptcMatrix A X) (S := Δ)
    (ys := ys) hKT hlin hjsr


/-- **Split Newton with a VARIABLE pseudo-time sequence: joint-spectral-radius criterion.**

The pseudo-times may change at every outer step: step `k` uses inverse pseudo-time
`c k = (δ k)⁻¹ ∈ Δ`, with `Δ ⊂ ℝ` compact (e.g. `Δ = [0, 1/δmin]`, keeping the steps bounded
below as in `limit_is_root`; `c k = 0` is a pure Newton step), and `A + c I` invertible for all
`c ∈ Δ`.

  **If the joint spectral radius of `{(A + c I)⁻¹ (c I + B D⁻¹ C) : c ∈ Δ}` is `< 1`,**

then for EVERY admissible sequence `c` and every start near `ys` the iteration converges
geometrically to `ys`, with constants independent of the sequence.  For a one-point `Δ` the joint
spectral radius is at most the spectral radius (`jointSpecRad_singleton_le`), recovering the
condition `ρ((A + δ⁻¹ I)⁻¹ (δ⁻¹ I + B D⁻¹ C)) < 1`. -/
theorem split_newton_variable_pseudo_time
    (G : Vec n × Vec m → Vec n) (H : Vec n × Vec m → Vec m) (ys : Vec n) (zs : Vec m)
    (G' : Vec n × Vec m →L[ℝ] Vec n) (H' : Vec n × Vec m →L[ℝ] Vec m)
    (hG : HasFDerivAt G G' (ys, zs)) (hH : HasFDerivAt H H' (ys, zs))
    (hG0 : G (ys, zs) = 0)
    (ζ : Vec n → Vec m) (ζ' : Vec n →L[ℝ] Vec m) (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys)
    (hsol : ∀ᶠ y in 𝓝 ys, H (y, ζ y) = 0) (hD : IsUnit (blockSnd H'))
    (Δ : Set ℝ) (hΔ : IsCompact Δ)
    (Ahat : Vec n → (Vec n →L[ℝ] Vec n)) (hA : ContinuousAt Ahat ys)
    (hAys : Ahat ys = blockFst G')
    (hres : ∀ c ∈ Δ, IsUnit (blockFst G' + c • (1 : Vec n →L[ℝ] Vec n)))
    (hjsr : jointSpecRad (ptcMatrix (blockFst G')
      ((blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H')))) Δ < 1) :
    ∃ ε > 0, ∃ C κ : ℝ, 0 ≤ κ ∧ κ < 1 ∧ ∀ c : ℕ → ℝ, (∀ k, c k ∈ Δ) →
      ∀ y0 : Vec n, ‖y0 - ys‖ < ε → ∀ k : ℕ,
        ‖switchedOrbit (ptcStep G ζ Ahat) c y0 k - ys‖ ≤ C * κ ^ k * ‖y0 - ys‖ :=
  split_newton_variable_pseudo_time_of_products G H ys zs G' H' hG hH hG0 ζ ζ' hζ0 hζ hsol hD Δ hΔ
    Ahat hA hAys hres (exists_prod_norm_lt_one_of_jsr hjsr)

/-- **Sufficient condition with a single norm.**  If every admissible linearised iteration matrix
`(A + c I)⁻¹ (c I + B D⁻¹ C)`, `c ∈ Δ`, is a strict contraction (`‖·‖ ≤ θ < 1`), the joint
condition holds with `N = 1`, so split Newton converges for every variable pseudo-time sequence
in `Δ`. -/
theorem split_newton_variable_pseudo_time_of_norm
    (G : Vec n × Vec m → Vec n) (H : Vec n × Vec m → Vec m) (ys : Vec n) (zs : Vec m)
    (G' : Vec n × Vec m →L[ℝ] Vec n) (H' : Vec n × Vec m →L[ℝ] Vec m)
    (hG : HasFDerivAt G G' (ys, zs)) (hH : HasFDerivAt H H' (ys, zs))
    (hG0 : G (ys, zs) = 0)
    (ζ : Vec n → Vec m) (ζ' : Vec n →L[ℝ] Vec m) (hζ0 : ζ ys = zs) (hζ : HasFDerivAt ζ ζ' ys)
    (hsol : ∀ᶠ y in 𝓝 ys, H (y, ζ y) = 0) (hD : IsUnit (blockSnd H'))
    (Δ : Set ℝ) (hΔ : IsCompact Δ)
    (Ahat : Vec n → (Vec n →L[ℝ] Vec n)) (hA : ContinuousAt Ahat ys)
    (hAys : Ahat ys = blockFst G')
    (hres : ∀ c ∈ Δ, IsUnit (blockFst G' + c • (1 : Vec n →L[ℝ] Vec n)))
    (θ : ℝ) (hθ : θ < 1)
    (hnorm : ∀ c ∈ Δ, ‖ptcMatrix (blockFst G')
      ((blockSnd G').comp ((Ring.inverse (blockSnd H')).comp (blockFst H'))) c‖ ≤ θ) :
    ∃ ε > 0, ∃ C κ : ℝ, 0 ≤ κ ∧ κ < 1 ∧ ∀ c : ℕ → ℝ, (∀ k, c k ∈ Δ) →
      ∀ y0 : Vec n, ‖y0 - ys‖ < ε → ∀ k : ℕ,
        ‖switchedOrbit (ptcStep G ζ Ahat) c y0 k - ys‖ ≤ C * κ ^ k * ‖y0 - ys‖ :=
  split_newton_variable_pseudo_time_of_products G H ys zs G' H' hG hH hG0 ζ ζ' hζ0 hζ hsol hD Δ hΔ
    Ahat hA hAys hres ⟨1, one_pos, θ, hθ, fun c hc => by
      simpa [prodFrom] using hnorm (c 0) (hc 0)⟩

end Main


end SplitNewton
