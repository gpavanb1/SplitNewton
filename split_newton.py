import logging
from numpy import Inf, concatenate
from newton import newton, criterion
from copy import deepcopy


def attach(x, y):
    return concatenate([x, y])

# Specify non-zero dt0 for pseudo-transient continuation
def split_newton(df, J, x0, loc, maxiter=Inf, sparse=False, dt0=0, dtmax=1., armijo=False):
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
    