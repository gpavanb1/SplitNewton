import logging
import numpy.linalg as la
from copy import deepcopy
from numpy import Inf, identity
from scipy import sparse as sp


def criterion(x, s, abs=1e-5, rel=1e-6):
    return la.norm(s/(x*rel + abs))



def newton(df, J, x0, maxiter=Inf, sparse=False, dt0=0., dtmax=1., armijo=False):
    """
    Unbounded Newton with pseudo-transient continuation
    (SER criterion for timestep modification) and Armijo rule
    https://ctk.math.ncsu.edu/TALKS/Purdue.pdf

    Operations are preferred to return numpy arrays

    Arguments
    ---------
        df: Function to compute gradient of objective
        J: Function to compute Jacobian of gradient
        x0: ndarray, Seed location
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

    x = deepcopy(x0)
    f0 = la.norm(df(x0))
    s = Inf
    crit = Inf

    iter = 0
    while crit >= 1 and iter < maxiter:

        # Update jacobian and function
        # Apply pseudo-transient continuation (diagonal modification) if required
        jac = J(x)
        if dt != 0:
            jac += (1/dt) * identity(len(x))
        dfx = df(x)
        fn = la.norm(dfx)

        # Choose solver based on sparsity
        if sparse:
            s, info = sp.linalg.gmres(sp.csr_matrix(jac), dfx)
            if info != 0:
                logging.warning("GMRES not converged")
        else:
            s = la.solve(jac, dfx)

        # Armijo rule
        if armijo:
            alpha = 1e-4
            for i in range(1, 10):
                fac = 2**-i
                newstep_norm = la.norm(df(x + fac*s))
                logging.debug(f"Armijo iteration: {fac}, {newstep_norm}, {fn}")
                if newstep_norm <= (1-alpha*fac)*fn:
                    break
            s *= fac
            logging.info(f"Armijo factor: {fac}")

        # Check convergence
        crit = criterion(x, s)
        logging.info(f"{x}, {s}, {crit}")
        if dt != 0:
            logging.info(f"Timestep: {dt}")

        # Update x
        x -= s
        iter += 1

        # Update timestep
        # No if condition as much faster directly
        dt = min(dt0*f0/fn, dtmax)
        
    return x, s, iter