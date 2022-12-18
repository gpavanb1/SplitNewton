**Theorem:** 

Let $\mathbf{x} := (\mathbf{y}, \mathbf{z})$ and $F(\mathbf{x}) := (G(\mathbf{y}, \mathbf{z}), H(\mathbf{y}, \mathbf{z}))$. Define a sequence

$\mathbf{x}_{k} := (\mathbf{y}_{k}, \mathbf{z}_{k})$ and $\mathbf{s}_{k} := \mathbf{x}_{k+1} - \mathbf{x}_{k} = (\mathbf{y}_{k+1} - \mathbf{y}_{k}, \mathbf{z}_{k+1} - \mathbf{z}_{k}) := (\mathbf{t}_{k}, \mathbf{u}_{k})$

where 

$\mathbf{t}_{k} = \{\delta_{k} \mathbf{I} - G'(\mathbf{y}_{k}, \mathbf{z}_{k})\}G(\mathbf{y}_{k}, \mathbf{z}_{k})$

$\mathbf{u}_{k} = \{\delta_{k} \mathbf{I} - H'(\mathbf{y}_{k}, \mathbf{z}_{k})\}H(\mathbf{y}_{k}, \mathbf{z}_{k})$

where $G'$ and $H'$ exists and invertible $\forall (\mathbf{y}_{k}, \mathbf{z}_{k})$

and $\delta_{k} = min(\delta_{0} ||F_{0}||/||F_{k}||, \delta_{max})$

Then, $\lim\limits_{k \to \infty} (\mathbf{y}_{k}, \lim\limits_{l \to \infty} z^{l}_{k}) = \mathbf{x}^{*}$, where $F(\mathbf{x}^{*}) = 0$

<br> 

**Proof:**

From convergence result in [1] for pseudo-transient continuation, on applying to $H(\mathbf{y}_{k}, \mathbf{z}^{l}_{k}) := H_{\mathbf{y}_{k}}(\mathbf{z}^{l}_{k})$

$\lim\limits_{k \to \infty} (\mathbf{y}_{k}, \lim\limits_{l \to \infty} \mathbf{z}^{l}_{k}) := \lim\limits_{k \to \infty} (\mathbf{y}_{k}, \mathbf{z}^{*}_{k})$, it follows that $\lim\limits_{k \to \infty} H_{\mathbf{y}_{k}}(\mathbf{z}_{k}^{*}) = 0$ and thus $\lim\limits_{k \to \infty} H(\mathbf{y}_{k}, \mathbf{z}^{*}_{k}) = 0$

Now, $\lim\limits_{k \to \infty}(G(\mathbf{y}_{k}, \mathbf{z}^{*}_{k}), H(\mathbf{y}_{k}, \mathbf{z}^{*}_{k})) = \lim\limits_{k \to \infty} (G(\mathbf{y}_{k}, \mathbf{z}^{*}_{k}), 0)$

Let $G(\mathbf{y}_{k}, \mathbf{z}^{*}_{k}) := G_{\mathbf{z}^{*}_{k}}(\mathbf{y}_{k})$

From convergence result in [1] applied to $G_{\mathbf{z}^{*}_{k}}(\mathbf{y}_{k})$, we have
$\lim\limits_{k \to \infty} G_{\mathbf{z}^{*}_{k}}(\mathbf{y}_{k}) := G_{\mathbf{z}^{*}}(\mathbf{y}^{*}) = 0$

Thus, $G(\mathbf{y}^{*}, \mathbf{z}^{*}) = 0$ and $H(\mathbf{y}^{*}, \mathbf{z}^{*}) = 0$ as $\lim\limits_{k \to \infty}H_{\mathbf{y}_{k}}(\mathbf{z}^{*}) = 0$




### References

[1] Kelley, Carl Timothy, and David E. Keyes. "Convergence analysis of pseudo-transient continuation." SIAM Journal on Numerical Analysis 35.2 (1998): 508-523.