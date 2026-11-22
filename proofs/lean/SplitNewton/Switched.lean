import Mathlib

/-!
# Local convergence of switched (time-varying) iterations; Ostrowski's theorem

Reference-level background, independent of Split-Newton.
* `switched_local_convergence` : `y_{k+1} = Φ_{i k}(y_k)`, `Φ_i` uniformly differentiable at the
  common fixed point `ys` with derivative `T_i`; if all `N`-fold products of the `T_i` have norm
  `≤ q < 1` then `ys` attracts nearby starts geometrically, uniformly in the switching pattern.
* `local_convergence`          : the single-map case (Ostrowski): `‖Φ'(ys)^N‖ < 1`.
-/

open Filter Topology
open scoped ENNReal NNReal Matrix

namespace SplitNewton

/-! ## 4. Local convergence for switched (time-varying) iterations

Engine of all convergence results below.  We iterate `y_{k+1} = Φ_{i k} (y_k)` where the index
`i k` (e.g. the pseudo-time of step `k`) may change at every step and ranges over a set `S`.
Each `Φ_i` fixes `ys` and is *uniformly* differentiable there with derivative `T_i`.
If all products of `N` consecutive derivatives (any switching pattern in `S`) have norm `≤ q < 1`
-- a norm form of "joint spectral radius `< 1`" -- then `ys` attracts all nearby starts
geometrically, uniformly in the switching pattern.  With a single index this is exactly
Ostrowski's theorem (`local_convergence` below). -/

section Switched

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] {ι : Type*}

/-- The orbit `y₀, Φ_{i 0} y₀, Φ_{i 1} (Φ_{i 0} y₀), …` of a time-varying iteration. -/
def switchedOrbit (Φ : ι → E → E) (i : ℕ → ι) (y0 : E) : ℕ → E
  | 0 => y0
  | k + 1 => Φ (i k) (switchedOrbit Φ i y0 k)

/-- Ordered product `T_{i (j-1)} ⋯ T_{i 1} T_{i 0}` of the first `j` linearisations. -/
noncomputable def prodFrom (T : ι → E →L[ℝ] E) (i : ℕ → ι) : ℕ → (E →L[ℝ] E)
  | 0 => 1
  | j + 1 => T (i j) * prodFrom T i j

/-- Running `m + n` steps = running `m` steps, then `n` steps with the shifted index sequence. -/
theorem switchedOrbit_add (Φ : ι → E → E) (i : ℕ → ι) (y0 : E) (m n : ℕ) :
    switchedOrbit Φ i y0 (m + n) =
      switchedOrbit Φ (fun k => i (m + k)) (switchedOrbit Φ i y0 m) n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    show Φ (i (m + n)) (switchedOrbit Φ i y0 (m + n)) =
      Φ (i (m + n)) (switchedOrbit Φ (fun k => i (m + k)) (switchedOrbit Φ i y0 m) n)
    rw [ih]

/-- **One block of `N` steps.**  Inside the region where the uniform linearisation estimate
holds, the true orbit stays within a factor `(K+2)^j` of the start, and differs from the
linearised orbit `T_{i(j-1)} ⋯ T_{i 0} (y₀ - ys)` by at most `j η (K+2)^j ‖y₀ - ys‖`
(each step contributes a remainder `≤ η‖y_j - ys‖`, propagated by operators of norm `≤ K`). -/
theorem switched_block {Φ : ι → E → E} {T : ι → E →L[ℝ] E} {S : Set ι} {ys : E}
    {K η ε : ℝ} (hK0 : 0 ≤ K) (hK : ∀ i ∈ S, ‖T i‖ ≤ K) (hη0 : 0 ≤ η) (hη1 : η ≤ 1)
    (hlin : ∀ i ∈ S, ∀ y, ‖y - ys‖ < ε → ‖Φ i y - ys - T i (y - ys)‖ ≤ η * ‖y - ys‖)
    (i : ℕ → ι) (hi : ∀ k, i k ∈ S) (y0 : E) (N : ℕ)
    (h0 : ‖y0 - ys‖ * (K + 2) ^ N < ε) :
    ∀ j ≤ N, ‖switchedOrbit Φ i y0 j - ys‖ ≤ (K + 2) ^ j * ‖y0 - ys‖ ∧
      ‖switchedOrbit Φ i y0 j - ys - prodFrom T i j (y0 - ys)‖
        ≤ j * η * (K + 2) ^ j * ‖y0 - ys‖ := by
  intro j
  induction j with
  | zero => intro _; simp [switchedOrbit, prodFrom]
  | succ j ih =>
    intro hj
    obtain ⟨ih1, ih2⟩ := ih (by omega)
    set w : E := switchedOrbit Φ i y0 j - ys with hw
    set W : ℝ := ‖y0 - ys‖ with hW
    have hW0 : 0 ≤ W := norm_nonneg _
    have hL1 : (1 : ℝ) ≤ K + 2 := by linarith
    -- the current point is still in the region where the linearisation estimate holds
    have hwε : ‖w‖ < ε := by
      calc ‖w‖ ≤ (K + 2) ^ j * W := ih1
        _ ≤ (K + 2) ^ N * W := mul_le_mul_of_nonneg_right
            (pow_le_pow_right₀ hL1 (by omega)) hW0
        _ < ε := by rw [mul_comm]; exact h0
    have he := hlin (i j) (hi j) _ (by simpa [hw, dist_eq_norm] using hwε)
    -- the one-step recursion
    set e : E := Φ (i j) (switchedOrbit Φ i y0 j) - ys - T (i j) w with hedef
    have hstep : switchedOrbit Φ i y0 (j + 1) - ys = T (i j) w + e := by
      simp only [switchedOrbit, hedef]; abel
    have hTw : ‖T (i j) w‖ ≤ K * ‖w‖ :=
      le_trans ((T (i j)).le_opNorm w) (mul_le_mul_of_nonneg_right (hK _ (hi j)) (norm_nonneg _))
    have he' : ‖e‖ ≤ η * ‖w‖ := by simpa [hw] using he
    constructor
    · rw [hstep]
      calc ‖T (i j) w + e‖ ≤ ‖T (i j) w‖ + ‖e‖ := norm_add_le _ _
        _ ≤ K * ‖w‖ + η * ‖w‖ := add_le_add hTw he'
        _ ≤ (K + 2) * ‖w‖ := by nlinarith [norm_nonneg w]
        _ ≤ (K + 2) * ((K + 2) ^ j * W) := mul_le_mul_of_nonneg_left ih1 (by linarith)
        _ = (K + 2) ^ (j + 1) * W := by ring
    · have hsplit : switchedOrbit Φ i y0 (j + 1) - ys - prodFrom T i (j + 1) (y0 - ys)
          = T (i j) (w - prodFrom T i j (y0 - ys)) + e := by
        rw [hstep]; simp only [prodFrom, ContinuousLinearMap.mul_apply, map_sub]; abel
      rw [hsplit]
      have hd : ‖T (i j) (w - prodFrom T i j (y0 - ys))‖
          ≤ K * (j * η * (K + 2) ^ j * W) :=
        le_trans ((T (i j)).le_opNorm _)
          (le_trans (mul_le_mul_of_nonneg_right (hK _ (hi j)) (norm_nonneg _))
            (mul_le_mul_of_nonneg_left ih2 hK0))
      have key : ((j : ℝ) + 1) * η * (K + 2) ^ (j + 1) * W
          - (K * (j * η * (K + 2) ^ j * W) + η * ((K + 2) ^ j * W))
          = η * (K + 2) ^ j * W * (2 * j + K + 1) := by ring
      have hpos : 0 ≤ η * (K + 2) ^ j * W * (2 * j + K + 1) := by positivity
      calc ‖T (i j) (w - prodFrom T i j (y0 - ys)) + e‖
          ≤ ‖T (i j) (w - prodFrom T i j (y0 - ys))‖ + ‖e‖ := norm_add_le _ _
        _ ≤ K * (j * η * (K + 2) ^ j * W) + η * ((K + 2) ^ j * W) :=
            add_le_add hd (le_trans he' (mul_le_mul_of_nonneg_left ih1 hη0))
        _ ≤ ((j + 1 : ℕ) : ℝ) * η * (K + 2) ^ (j + 1) * W := by push_cast; linarith

/-- **Local convergence of a switched iteration (uniform in the switching pattern).**

Hypotheses: the `Φ_i` (`i ∈ S`) are uniformly differentiable at `ys` with derivatives `T_i` of
bounded norm, and some `N`-fold product of the `T_i` (along *any* sequence in `S`) has norm
`≤ q < 1`.  Conclusion: from any start close enough to `ys`, and for *every* sequence `i` in `S`,
the orbit converges to `ys` geometrically, with constants independent of `i`. -/
theorem switched_local_convergence {Φ : ι → E → E} {T : ι → E →L[ℝ] E} {S : Set ι} {ys : E}
    (hK : ∃ K : ℝ, ∀ i ∈ S, ‖T i‖ ≤ K)
    (hlin : ∀ η : ℝ, 0 < η → ∃ ε : ℝ, 0 < ε ∧ ∀ i ∈ S, ∀ y : E, ‖y - ys‖ < ε →
      ‖Φ i y - ys - T i (y - ys)‖ ≤ η * ‖y - ys‖)
    (hN : ∃ N : ℕ, 0 < N ∧ ∃ q : ℝ, q < 1 ∧
      ∀ i : ℕ → ι, (∀ k, i k ∈ S) → ‖prodFrom T i N‖ ≤ q) :
    ∃ ε : ℝ, 0 < ε ∧ ∃ C κ : ℝ, 0 ≤ κ ∧ κ < 1 ∧ ∀ i : ℕ → ι, (∀ k, i k ∈ S) →
      ∀ y0 : E, ‖y0 - ys‖ < ε → ∀ k : ℕ,
        ‖switchedOrbit Φ i y0 k - ys‖ ≤ C * κ ^ k * ‖y0 - ys‖ := by
  obtain ⟨K, hK⟩ := hK
  obtain ⟨N, hN0, q, hq1, hq⟩ := hN
  set K' : ℝ := max K 0 with hK'
  have hK'0 : 0 ≤ K' := le_max_right _ _
  have hKK' : ∀ i ∈ S, ‖T i‖ ≤ K' := fun i hi => le_trans (hK i hi) (le_max_left _ _)
  set L : ℝ := K' + 2 with hL
  have hL1 : (1 : ℝ) ≤ L := by linarith
  have hLpos : 0 < L := by linarith
  have hLN : 0 < L ^ N := pow_pos hLpos N
  -- q₀ = max q 0 ∈ [0,1), contraction factor per block q' = (1 + q₀)/2 ∈ (0,1)
  set q0 : ℝ := max q 0 with hq0
  have hq00 : 0 ≤ q0 := le_max_right _ _
  have hq01 : q0 < 1 := max_lt hq1 one_pos
  set q' : ℝ := (1 + q0) / 2 with hq'
  have hq'pos : 0 < q' := by linarith
  have hq'1 : q' < 1 := by linarith
  -- tolerance of the linearisation: N η L^N ≤ (1 - q₀)/2
  set η : ℝ := min 1 ((1 - q0) / (2 * ((N : ℝ) + 1) * L ^ N)) with hη
  have hNpos : (0 : ℝ) < (N : ℝ) + 1 := by positivity
  have hηpos : 0 < η := lt_min one_pos (by
    have : 0 < 1 - q0 := by linarith
    positivity)
  have hη1 : η ≤ 1 := min_le_left _ _
  have hηN : (N : ℝ) * η * L ^ N ≤ (1 - q0) / 2 := by
    have h1 : η ≤ (1 - q0) / (2 * ((N : ℝ) + 1) * L ^ N) := min_le_right _ _
    calc (N : ℝ) * η * L ^ N
        ≤ (N : ℝ) * ((1 - q0) / (2 * ((N : ℝ) + 1) * L ^ N)) * L ^ N := by
          gcongr
      _ = (N : ℝ) * (1 - q0) / (2 * ((N : ℝ) + 1)) := by field_simp
      _ ≤ (1 - q0) / 2 := by
          rw [div_le_div_iff₀ (by positivity) (by norm_num)]
          have : 0 ≤ 1 - q0 := by linarith
          nlinarith [(Nat.cast_nonneg N : (0 : ℝ) ≤ N)]
  obtain ⟨ε, hε, hlin'⟩ := hlin η hηpos
  -- starting radius for which a whole block stays inside the linearisation region
  set ε0 : ℝ := ε / L ^ N with hε0
  have hε0pos : 0 < ε0 := by positivity
  -- one block of N steps contracts by q' and never grows by more than L^N
  have block : ∀ i : ℕ → ι, (∀ k, i k ∈ S) → ∀ y0 : E, ‖y0 - ys‖ < ε0 →
      (∀ j ≤ N, ‖switchedOrbit Φ i y0 j - ys‖ ≤ L ^ j * ‖y0 - ys‖) ∧
      ‖switchedOrbit Φ i y0 N - ys‖ ≤ q' * ‖y0 - ys‖ := by
    intro i hi y0 hy0
    have h0 : ‖y0 - ys‖ * (K' + 2) ^ N < ε := by
      have := mul_lt_mul_of_pos_right hy0 hLN
      rw [hε0, div_mul_cancel₀ _ hLN.ne'] at this
      exact this
    have hb := switched_block hK'0 hKK' hηpos.le hη1 hlin' i hi y0 N h0
    refine ⟨fun j hj => (hb j hj).1, ?_⟩
    obtain ⟨_, hb2⟩ := hb N le_rfl
    have hprod : ‖prodFrom T i N (y0 - ys)‖ ≤ q0 * ‖y0 - ys‖ :=
      le_trans ((prodFrom T i N).le_opNorm _)
        (mul_le_mul_of_nonneg_right (le_trans (hq i hi) (le_max_left _ _)) (norm_nonneg _))
    have hsplit : switchedOrbit Φ i y0 N - ys
        = prodFrom T i N (y0 - ys) + (switchedOrbit Φ i y0 N - ys - prodFrom T i N (y0 - ys)) := by
      abel
    calc ‖switchedOrbit Φ i y0 N - ys‖
        ≤ ‖prodFrom T i N (y0 - ys)‖ + ‖switchedOrbit Φ i y0 N - ys - prodFrom T i N (y0 - ys)‖ := by
          exact (congrArg norm hsplit).le.trans (norm_add_le _ _)
      _ ≤ q0 * ‖y0 - ys‖ + N * η * L ^ N * ‖y0 - ys‖ := add_le_add hprod hb2
      _ ≤ q0 * ‖y0 - ys‖ + (1 - q0) / 2 * ‖y0 - ys‖ := by
          have := mul_le_mul_of_nonneg_right hηN (norm_nonneg (y0 - ys))
          linarith
      _ = q' * ‖y0 - ys‖ := by rw [hq']; ring
  -- geometric rate per step
  set κ : ℝ := q' ^ ((N : ℝ)⁻¹) with hκ
  have hκ0 : 0 < κ := Real.rpow_pos_of_pos hq'pos _
  have hκ1 : κ < 1 := Real.rpow_lt_one hq'pos.le hq'1 (inv_pos.2 (by exact_mod_cast hN0))
  have hκN : κ ^ N = q' := Real.rpow_inv_natCast_pow hq'pos.le (by omega)
  refine ⟨ε0, hε0pos, L ^ N / q', κ, hκ0.le, hκ1, ?_⟩
  intro i hi y0 hy0 k
  -- after `a` blocks the distance has shrunk by `q'^a`
  have hblocks : ∀ a : ℕ, ‖switchedOrbit Φ i y0 (a * N) - ys‖ ≤ q' ^ a * ‖y0 - ys‖ := by
    intro a
    induction a with
    | zero => simp [switchedOrbit]
    | succ a ih =>
      have hlt : ‖switchedOrbit Φ i y0 (a * N) - ys‖ < ε0 := by
        refine lt_of_le_of_lt ih (lt_of_le_of_lt ?_ hy0)
        calc q' ^ a * ‖y0 - ys‖ ≤ 1 * ‖y0 - ys‖ :=
              mul_le_mul_of_nonneg_right (pow_le_one₀ hq'pos.le hq'1.le) (norm_nonneg _)
          _ = ‖y0 - ys‖ := one_mul _
      have hb := (block (fun k => i (a * N + k)) (fun k => hi _) _ hlt).2
      rw [Nat.succ_mul, switchedOrbit_add]
      calc _ ≤ q' * ‖switchedOrbit Φ i y0 (a * N) - ys‖ := hb
        _ ≤ q' * (q' ^ a * ‖y0 - ys‖) := mul_le_mul_of_nonneg_left ih hq'pos.le
        _ = q' ^ (a + 1) * ‖y0 - ys‖ := by ring
  -- write k = a * N + j with j < N
  have hk : k = (k / N) * N + k % N := (Nat.div_add_mod' k N).symm
  have hjN : k % N < N := Nat.mod_lt _ hN0
  set a := k / N with ha
  set j := k % N with hj
  have hlt : ‖switchedOrbit Φ i y0 (a * N) - ys‖ < ε0 := by
    refine lt_of_le_of_lt (hblocks a) (lt_of_le_of_lt ?_ hy0)
    calc q' ^ a * ‖y0 - ys‖ ≤ 1 * ‖y0 - ys‖ :=
          mul_le_mul_of_nonneg_right (pow_le_one₀ hq'pos.le hq'1.le) (norm_nonneg _)
      _ = ‖y0 - ys‖ := one_mul _
  have hbj := (block (fun k => i (a * N + k)) (fun k => hi _) _ hlt).1 j hjN.le
  have hkeq : switchedOrbit Φ i y0 k = switchedOrbit Φ (fun k => i (a * N + k))
      (switchedOrbit Φ i y0 (a * N)) j := by
    conv_lhs => rw [hk]
    exact switchedOrbit_add _ _ _ _ _
  rw [hkeq]
  have hLj : L ^ j ≤ L ^ N := pow_le_pow_right₀ hL1 hjN.le
  -- q'^a * q' ≤ κ^k
  have hqa : q' ^ a * q' ≤ κ ^ k := by
    have h1 : q' ^ a * q' = κ ^ (N * a + N) := by
      rw [pow_add, pow_mul, hκN]
    rw [h1]
    exact pow_le_pow_of_le_one hκ0.le hκ1.le (by
      have : a * N ≤ N * a := by rw [mul_comm]
      omega)
  have hqa' : q' ^ a ≤ κ ^ k / q' := by
    rw [le_div_iff₀ hq'pos]; exact hqa
  calc ‖switchedOrbit Φ (fun k => i (a * N + k)) (switchedOrbit Φ i y0 (a * N)) j - ys‖
      ≤ L ^ j * ‖switchedOrbit Φ i y0 (a * N) - ys‖ := hbj
    _ ≤ L ^ N * (q' ^ a * ‖y0 - ys‖) :=
        mul_le_mul hLj (hblocks a) (norm_nonneg _) hLN.le
    _ ≤ L ^ N * ((κ ^ k / q') * ‖y0 - ys‖) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hqa' (norm_nonneg _)) hLN.le
    _ = L ^ N / q' * κ ^ k * ‖y0 - ys‖ := by ring

end Switched

/-! ### Special case: a single map (Ostrowski's theorem) -/

section Ostrowski

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- With a constant index sequence the product of `j` linearisations is the power `T ^ j`. -/
theorem prodFrom_const (T : E →L[ℝ] E) (i : ℕ → Unit) (j : ℕ) :
    prodFrom (fun _ : Unit => T) i j = T ^ j := by
  induction j with
  | zero => simp [prodFrom]
  | succ j ih => rw [prodFrom, ih, pow_succ']

/-- With a constant index sequence the switched orbit is the usual iterate `Φ^[k] y₀`. -/
theorem switchedOrbit_const (Φ : E → E) (i : ℕ → Unit) (y0 : E) (k : ℕ) :
    switchedOrbit (fun _ : Unit => Φ) i y0 k = Φ^[k] y0 := by
  induction k with
  | zero => rfl
  | succ k ih => rw [Function.iterate_succ_apply', ← ih]; rfl

/-- **Ostrowski-type local convergence (general nonlinear map, any real normed space).**
If `Φ ys = ys`, `Φ` is Fréchet differentiable at `ys` with derivative `T`, and some power of
`T` has norm `< 1` (e.g. `ρ(T) < 1`, see `exists_pow_norm_lt_one`), then every start close
enough to `ys` converges to `ys` *geometrically*:
`‖Φᵏ y₀ - ys‖ ≤ C κᵏ ‖y₀ - ys‖` with `κ < 1`. -/
theorem local_convergence {Φ : E → E} {T : E →L[ℝ] E} {ys : E}
    (hfix : Φ ys = ys) (hΦ : HasFDerivAt Φ T ys)
    {N : ℕ} (hN : 0 < N) (hT : ‖T ^ N‖ < 1) :
    ∃ ε > 0, ∃ C : ℝ, ∃ κ : ℝ, 0 ≤ κ ∧ κ < 1 ∧ ∀ y0 : E, ‖y0 - ys‖ < ε →
      ∀ k : ℕ, ‖Φ^[k] y0 - ys‖ ≤ C * κ ^ k * ‖y0 - ys‖ := by
  -- the linearisation estimate, in epsilon-delta form
  have hlin : ∀ η : ℝ, 0 < η → ∃ ε : ℝ, 0 < ε ∧ ∀ y : E, ‖y - ys‖ < ε →
      ‖Φ y - ys - T (y - ys)‖ ≤ η * ‖y - ys‖ := by
    intro η hη
    have hlo := hΦ.isLittleO.def hη
    rw [hfix] at hlo
    obtain ⟨ε, hε, hball⟩ := Metric.eventually_nhds_iff_ball.1 hlo
    exact ⟨ε, hε, fun y hy => by
      simpa using hball y (by simpa [dist_eq_norm] using hy)⟩
  obtain ⟨ε, hε, C, κ, hκ0, hκ1, hb⟩ := switched_local_convergence
    (Φ := fun _ : Unit => Φ) (T := fun _ : Unit => T) (S := Set.univ) (ys := ys)
    ⟨‖T‖, fun _ _ => le_rfl⟩
    (fun η hη => by
      obtain ⟨ε, hε, h⟩ := hlin η hη
      exact ⟨ε, hε, fun _ _ y hy => h y hy⟩)
    ⟨N, hN, ‖T ^ N‖, hT, fun i _ => by rw [prodFrom_const]⟩
  refine ⟨ε, hε, C, κ, hκ0, hκ1, fun y0 hy0 k => ?_⟩
  have := hb (fun _ => ()) (fun _ => trivial) y0 hy0 k
  rwa [switchedOrbit_const] at this

end Ostrowski

end SplitNewton
