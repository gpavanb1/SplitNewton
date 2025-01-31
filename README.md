# SplitNewton

[![Downloads](https://pepy.tech/badge/splitnewton)](https://pepy.tech/project/splitnewton)
![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14782293.svg)](https://doi.org/10.5281/zenodo.14782293)

Bounded, SPLIT [Newton](https://en.wikipedia.org/wiki/Newton%27s_method) with [pseudo-transient continuation
](https://ctk.math.ncsu.edu/TALKS/Purdue.pdf) and [backtracking](https://en.wikipedia.org/wiki/Backtracking_line_search)

Good for ill-conditioned problems where there are two different sets of systems

Particular applications include
* [Fast-Slow Reaction-Diffusion systems](https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system)
* [CFD](https://en.wikipedia.org/wiki/Computational_fluid_dynamics) - Pressure-Velocity coupling

## What does 'split' mean?

The system is divided into multiple segments, and for ease of communication, let’s refer to the first segment of variables as "outer" and the remaining as "inner".

* Holding the outer variables fixed, Newton iteration is performed recursively for the inner variables, using the sub-Jacobian associated with them, until convergence is reached.

* One Newton step is then performed for the outer variables, while the inner variables are kept fixed, using the sub-Jacobian for the outer subsystem.

* This process is repeated, alternating between solving the inner and outer subsystems, until the convergence criterion for the entire system (similar to standard Newton) is met.

### Example:

Consider a system of 5 variables, with the split locations at indices [1, 4]. This results in the following segments:

  * `a1` (variables from 0 to 1)
  * `a2 a3 a4` (variables from 1 to 4)
  * `a5` (variable at index 4)

1. First, the innermost segment `a5` is solved recursively using Newton's method while holding the variables `a1` and `a2 a3 a4`) fixed. This step is repeated until the convergence criterion for `a5` is met.

2. Next, one Newton step is taken for the segment `a2 a3 a4`, with `a5` held fixed. This step is followed by solving `a5` again till convergence.

3. This alternating process repeats: solving for `a5` until convergence, then one step for `a2 a3 a4`, and so on, until all subsystems converge.

Finally, one Newton step is performed for `a1`, with the other segments fixed. This completes one cycle of the split Newton process.

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
with the second system,
$\lambda_{c} = 10^{-1}$
$\lambda_{d} = 10^{-4}$
and third system,
$\lambda_{c} = 10^{-6}$
$\lambda_{d} = 10^{-8}$

and using `logspace` for variation in $\lambda_{i}$


$$ F(u) = \lambda_{a} u^{4}_{1} + ... + \lambda_{b} u^{4}_{\lfloor N/3 \rfloor} + \lambda_{c} u^{4}_{\lceil N/3 \rceil} + ... + \lambda_{d} u^{4}_{\lfloor 2N/3 \rfloor} + \lambda_{e} u^{4}_{\lceil 2N/3 \rceil} + ... + \lambda_{f} u^{4}_{N}$$

$$
J(u) = 3 \times \begin{bmatrix}
\lambda_a & \dots & 0 & 0 & \dots & 0 & 0 & \dots & 0 \\
\vdots & \ddots & \vdots & \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\
0 & \dots & \lambda_b & 0 & \dots & 0 & 0 & \dots & 0 \\
0 & \dots & 0 & \lambda_c & \dots & 0 & 0 & \dots & 0 \\
\vdots & \ddots & \vdots & \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & 0 & \dots & \lambda_d & 0 & \dots & 0 \\
0 & \dots & 0 & 0 & \dots & 0 & \lambda_e & \dots & 0 \\
\vdots & \ddots & \vdots & \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & 0 & \dots & 0 & 0 & \dots & \lambda_f
\end{bmatrix} \cdot u^2
$$

For N=5000 (with no backtracking and pseudo-transient continuation), 

| Method    | Time       | Iterations    |
|-----------|------------|---------------|
| Split Newton    |    34 seconds |  33   |
| Newton |  not converged > 1 min  | NA  |

## How to test?
You can run tests with the `pytest` framework using `python -m pytest`

The coverage reports can be generated with `pytest-cov` plugin using `python -m pytest --cov=splitnewton`

## Whom to contact?

Please direct your queries to [gpavanb1](http://github.com/gpavanb1)
for any questions.

## Citing

If you are using `SplitNewton` in any scientific work, please make sure to cite as follows
```
@software{pavan_b_govindaraju_2025_14782293,
  author       = {Pavan B Govindaraju},
  title        = {gpavanb1/SplitNewton: v0.3.1},
  month        = jan,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v0.3.1},
  doi          = {10.5281/zenodo.14782293},
  url          = {https://doi.org/10.5281/zenodo.14782293},
}
```
