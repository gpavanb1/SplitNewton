import logging
from copy import deepcopy
from numpy import Inf, concatenate
from splitnewton.newton import newton, criterion


def attach(x, y):
    return concatenate([x, y])

# Specify non-zero dt0 for pseudo-transient continuation
def split_newton(df, J, x0, loc, maxiter=Inf, sparse=False, dt0=0, dtmax=1., armijo=False):
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

    s = Inf

    crit = Inf

    iter = 0
    while crit >= 1  and iter < maxiter:
        # B Cycle
        dfb = lambda x: df(attach(xa, x))[loc:]
        Jb = lambda x: J(attach(xa, x))[loc:, loc:]
        xb, sb, local_iter = newton(dfb, Jb, xb, maxiter, sparse, dt, dtmax, armijo)
        logging.debug(f"B cycle: {xb}, {sb}")
        logging.debug(f"B iterations: {local_iter}")

        # A Cycle
        dfa = lambda x: df(attach(x, xb))[:loc]
        Ja = lambda x: J(attach(x, xb))[:loc, :loc]
        xa, sa, local_iter = newton(dfa, Ja, xa, 1, sparse, dt, dtmax, armijo)
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
    