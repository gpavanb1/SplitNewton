from numpy import Inf
import numpy.linalg as la
from copy import deepcopy


def criterion(x, s, abs=1e-5, rel=1e-6):
    return la.norm(s/(x*rel + abs))


# Find df = 0 or extrema of f
def newton(df, J, x0, maxiter=Inf):
    x = deepcopy(x0)
    s = Inf
    crit = Inf

    iter = 0
    while crit >= 1 and iter < maxiter:
        s = la.solve(J(x), df(x))
        crit = criterion(x, s)
        print(x, s, crit)
        x -= s
        iter += 1
        
    return x, s