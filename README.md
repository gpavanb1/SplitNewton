# SplitNewton

[![Downloads](https://pepy.tech/badge/splitnewton)](https://pepy.tech/project/splitnewton)

Unbounded SPLIT [Newton](https://en.wikipedia.org/wiki/Newton%27s_method) with [pseudo-transient continuation
](https://ctk.math.ncsu.edu/TALKS/Purdue.pdf) and [backtracking](https://en.wikipedia.org/wiki/Backtracking_line_search)

Good for ill-conditioned problems where there are two different sets of systems

Particular applications include
* [Fast-Slow Reaction-Diffusion systems](https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system)
* [CFD](https://en.wikipedia.org/wiki/Computational_fluid_dynamics) - Pressure-Velocity coupling


## How to install and execute?

Just run 
```
pip install splitnewton
```

There is an [examples](https://github.com/gpavanb1/SplitNewton/examples) folder that contains a test function and driver program

## How good is this?

Consider the test problem

$\lambda_{a} = 10^{6}$, 
$\lambda_{b} = 10^{2}$

and the second system
$\lambda_{c} = 10^{-1}$
$\lambda_{d} = 10^{-4}$

$$
u' = 0.25 * \begin{bmatrix}
\lambda_a & \dots & 0 & 0 & \dots & 0 \newline
\vdots & \ddots & \vdots & \vdots & \ddots & \vdots \newline
0 & \dots & \lambda_b & 0 & \dots & 0 \newline
0 & \dots & 0 & \lambda_c & \dots & 0 \newline
\vdots & \ddots & \vdots & \vdots & \ddots & \vdots \newline
0 & \dots & 0 & 0 & \dots & \lambda_d
\end{bmatrix} u^{4}
$$

For N=5000 (with no backtracking and pseudo-transient continuation), 

| Method    | Time       | Iterations    |
|-----------|------------|---------------|
| Split Newton    |    9 seconds |  32   |
| Newton |  not converged > 1 min  | NA  |


## Whom to contact?

Please direct your queries to [gpavanb1](http://github.com/gpavanb1)
for any questions.