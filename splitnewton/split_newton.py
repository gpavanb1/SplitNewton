import logging
from copy import deepcopy
from numpy import inf, concatenate
from splitnewton.newton import newton, criterion


def attach(x, y):
    return concatenate([x, y])

# Specify non-zero dt0 for pseudo-transient continuation


def split_newton(df, J, x0, loc, maxiter=inf, sparse=False, dt0=0, dtmax=1., armijo=False, bounds=None, bound_fac=0.8):
    """
    Unbounded SPLIT Newton with pseudo-transient continuation
    (SER criterion for timestep modification) and Armijo rule
    https://ctk.math.ncsu.edu/TALKS/Purdue.pdf

    Operations are preferred to return numpy arrays

    Arguments
    ---------
        df: Function to compute gradient of objective
        J: Function to compute Jacobian of gradient
        x0: ndarray, Seed location
        loc: int, Location at which to split the system
        maxiter: int, Maximum number of iterations
        sparse: bool, Use sparse or dense linear solver
        dt0: float, Initial pseudo-timestep size
        dtmax: float, Maximum pseudo-timestep size
        armijo: bool, Apply Armijo rule to choose step fraction
        bounds: list of 2-lists, each of size len(x) that contain lower and upper bound
        bound_fac: damping factor by which solver must step to avoid crossing limits

    Returns
    -------
        x: ndarray, Final solution
        s: ndarray, Final step
        iter: int, Number of iterations
    """
    if dt0 < 0 or dtmax < 0:
        raise Exception("Must specify positive dt0 and dtmax")
    dt = dt0

    if loc > len(x0):
        raise Exception('Incorrect split location')

    xa = deepcopy(x0[:loc])
    xb = deepcopy(x0[loc:])
    x = deepcopy(x0)

    s = inf

    crit = inf

    iter = 0
    while crit >= 1 and iter < maxiter:
        # B Cycle
        def dfb(x): return df(attach(xa, x))[loc:]
        def Jb(x): return J(attach(xa, x))[loc:, loc:]
        local_bounds = [bounds[0][loc:], bounds[1]
                        [loc:]] if bounds is not None else None
        xb, sb, local_iter = newton(
            dfb, Jb, xb, maxiter, sparse, dt, dtmax, armijo, local_bounds)
        logging.debug(f"B cycle: {xb}, {sb}")
        logging.debug(f"B iterations: {local_iter}")

        # A Cycle
        def dfa(x): return df(attach(x, xb))[:loc]
        def Ja(x): return J(attach(x, xb))[:loc, :loc]
        local_bounds = [bounds[0][:loc], bounds[1]
                        [:loc]] if bounds is not None else None
        xa, sa, local_iter = newton(
            dfa, Ja, xa, 1, sparse, dt, dtmax, armijo, local_bounds)
        logging.debug(f"A cycle: {xa}, {sa}")
        logging.debug(f"A iterations: {local_iter}")

        # Construct new x and step
        xnew = attach(xa, xb)
        s = xnew - x

        # Check convergence
        crit = criterion(x, s)
        logging.info(f"{x}, {s}, {crit}")

        # Update x
        x = xnew
        iter += 1

    return x, s, iter
