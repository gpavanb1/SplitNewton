import logging
import numpy.linalg as la
from copy import deepcopy
from numpy import inf, identity
from scipy import sparse as sp
from math import ulp

EPS = ulp(1.0)


def criterion(x, s, abs=1e-5, rel=1e-6):
    return la.norm(s/(x*rel + abs))


def check_within_bounds(x0, bounds):
    if bounds is None:
        return True

     # Check if valid format of bounds
    if len(bounds) != 2:
        raise Exception(
            "Bounds must be a 2 lists, with lower and upper bounds of each dimension")

    lower, upper = bounds
    if len(lower) != len(x0) or len(upper) != len(x0):
        raise Exception(
            "Each bounds list must be as long as the solution vector")

    for i in range(len(x0)):
        if not (lower[i] <= x0[i] <= upper[i]):
            return False
    return True


def newton(df, J, x0, maxiter=inf, sparse=False, dt0=0., dtmax=1., armijo=False, bounds=None, bound_fac=0.8):
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

    # Check if seed within bounds
    if not check_within_bounds(x0, bounds):
        raise Exception("Seed must be within the provided bounds")

    x = deepcopy(x0)
    f0 = la.norm(df(x0))
    s = inf
    crit = inf

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
        # Note that the final step is -s
        if sparse:
            s, info = sp.linalg.gmres(sp.csr_matrix(jac), dfx, atol='legacy')
            if info != 0:
                logging.warning("GMRES not converged")
        else:
            s = la.solve(jac, dfx)

        # Apply sign convention
        s = -s

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

        ####
        # Apply bounds
        ####
        if bounds is not None:
            lower, upper = bounds
            # Set the scaling factor as the highest damping that must be done
            best_scaling = 1.0
            for i in range(len(x)):
                if s[i] == 0.0:
                    continue
                # Move only atmost `legal_delta` in each dimension
                else:
                    # Calculate the new position and legal delta
                    legal_delta = bound_fac * \
                        (lower[i] - x[i]) if s[i] < 0 else bound_fac * \
                        (upper[i] - x[i])
                    # Both legal_delta and s[i] will be of the same sign
                    # If it exceeds 1, the step is retained, else
                    # the scaling factor is updated
                    scaling = legal_delta / s[i]
                    best_scaling = min(best_scaling, scaling)

            # Show damping factor if invoked
            if best_scaling != 1.0:
                logging.info(f"Bounds hit. Damping factor: {best_scaling}")
            s *= best_scaling

        # Check convergence
        crit = criterion(x, s)
        logging.info(f"{x}, {s}, {crit}")
        if dt != 0:
            logging.info(f"Timestep: {dt}")

        # Update x
        x += s
        iter += 1

        # Update timestep
        # No if condition as much faster directly
        dt = min(dt0*f0/(fn + EPS), dtmax)

    return x, s, iter
