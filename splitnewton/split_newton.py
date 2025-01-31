import logging
from copy import deepcopy
from numpy import inf, concatenate
from scipy.sparse import csr_matrix
from splitnewton.newton import newton, criterion


def split_newton(df, J, x0, locs, maxiter=100, sparse=True, dt0=0., dtmax=1e-1, armijo=False, bounds=None):
    """
    Solve a nonlinear system using a hierarchical split Newton method.

    Arguments:
        df: Function to compute the residual.
        J: Function to compute the Jacobian.
        x0: Initial guess (full system state).
        locs: List of split locations.
        maxiter: Maximum number of global iterations.
        sparse: Whether to use sparse matrices.
        dt0, dtmax: Step size control for Newton.
        armijo: Armijo line search flag.
        bounds: Tuple (lower, upper bounds) for constrained optimization.

    Returns:
        x: Solution vector.
        s: Step size vector.
        iter: Number of iterations.
    """

    # Validate split locations
    for i in locs:
        if i < 0 or i > len(x0):
            raise ValueError('Incorrect split location')

    # Base case: if no splits remain, solve with Newton directly
    if len(locs) == 0:
        return newton(df, J, x0, maxiter, sparse, dt0, dtmax, armijo, bounds)

    # Split x0 into xa and xb
    loc = locs[0]
    xa, xb = deepcopy(x0[:loc]), deepcopy(x0[loc:])

    # Define residual and Jacobian for xa (keeping xb fixed)
    def dfa(x): return df(concatenate((x, xb)))[:loc]

    def Ja(x):
        J_matrix = J(concatenate((x, xb)))[:loc, :loc]
        return csr_matrix(J_matrix) if sparse else J_matrix

    # Define residual and Jacobian for xb (keeping xa fixed)
    def dfb(x): return df(concatenate((xa, x)))[loc:]

    def Jb(x):
        J_matrix = J(concatenate((xa, x)))[loc:, loc:]
        return csr_matrix(J_matrix) if sparse else J_matrix

    # Adjust locs for recursion (relative to xb)
    new_locs = [l - loc for l in locs[1:]]

    # Adjust bounds for recursion
    bounds_a = ([b[:loc] for b in bounds] if bounds else None)
    bounds_b = ([b[loc:] for b in bounds] if bounds else None)

    x, s, crit, iter = deepcopy(x0), inf, inf, 0

    while crit >= 1 and iter < maxiter:
        # Solve rightmost subsystem recursively
        xb, _, _ = split_newton(dfb, Jb, xb, new_locs,
                                maxiter, sparse, dt0, dtmax, armijo, bounds_b)

        # One Newton step for left subsystem
        xa, _, _ = newton(dfa, Ja, xa, 1, sparse, dt0, dtmax, armijo, bounds_a)

        # Construct full x and check convergence
        xnew = concatenate((xa, xb))
        s, crit = xnew - x, criterion(x, s)

        logging.info(f"{x}, {s}, {crit}")
        x, iter = xnew, iter + 1

    return x, s, iter
