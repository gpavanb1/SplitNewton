from numpy import Inf
import numpy.linalg as la
from scipy import sparse as sp
from copy import deepcopy


def criterion(x, s, abs=1e-5, rel=1e-6):
    return la.norm(s/(x*rel + abs))


# Find df = 0 or extrema of f
def newton(df, J, x0, maxiter=Inf, sparse=False):
    x = deepcopy(x0)
    s = Inf
    crit = Inf

    iter = 0
    while crit >= 1 and iter < maxiter:
        if sparse:
            s, _ = sp.linalg.gmres(sp.csr_matrix(J(x)), df(x))
        else:
            s = la.solve(J(x), df(x))
        crit = criterion(x, s)
        print(x, s, crit)
        x -= s
        iter += 1
        
    return x, s