from numpy import Inf
import numpy.linalg as la


def criterion(x, s, abs=1e-5, rel=1e-6):
    return la.norm(s/(x*rel + abs))


def newton(f, df, J, x0):
    x = x0
    s = Inf
    crit = Inf

    while crit >= 1:
        s = la.solve(J(x), df(x))
        crit = criterion(x, s)
        print(x, s, crit)
        x -= s
        
    return x