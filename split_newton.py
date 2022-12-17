from numpy import Inf, concatenate
from newton import newton, criterion
from copy import deepcopy


def attach(x, y):
    return concatenate([x, y])

def split_newton(df, J, x0, loc, maxiter=100, sparse=False):
    if loc > len(x0):
        raise Exception('Incorrect split location')

    xa = deepcopy(x0[:loc])
    xb = deepcopy(x0[loc:])
    x = deepcopy(x0)

    s = Inf

    crit = Inf

    while crit >= 1:
        dfb = lambda x: df(attach(xa, x))[loc:]
        Jb = lambda x: J(attach(xa, x))[loc:, loc:]
        xb, sb = newton(dfb, Jb, xb, maxiter, sparse=sparse)
        print("B", xb, sb)

        dfa = lambda x: df(attach(x, xb))[:loc]
        Ja = lambda x: J(attach(x, xb))[:loc, :loc]
        xa, sa = newton(dfa, Ja, xa, 1, sparse=sparse)
        print("A", xa, sa)

        xnew = attach(xa, xb)
        s = xnew - x
        crit = criterion(x, s)
        print(x, s, crit)
        x = xnew

    return x, s
    