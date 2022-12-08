from numpy import Inf, concatenate
from newton import newton, criterion


def attach(x, y):
    return concatenate([x, y])

def split_newton(df, J, x0, loc, maxiter=100):
    if loc > len(x0):
        raise Exception('Incorrect split location')

    xa = x0[:loc]
    xb = x0[loc:]
    x = x0

    sa = Inf
    sb = Inf
    s = Inf

    crita = Inf
    critb = Inf

    while criterion(x, s) >= 1:
        dfa = lambda x: df(attach(x, xb))[:loc]
        Ja = lambda x: J(attach(x, xb))[:loc, :loc]
        xa, sa = newton(dfa, Ja, xa, maxiter)
        crita = criterion(xa, sa)
        print("A", xa, sa, crita)

        dfb = lambda x: df(attach(xa, x))[loc:]
        Jb = lambda x: J(attach(xa, x))[loc:, loc:]
        xb, sb = newton(dfb, Jb, xb, maxiter)
        critb = criterion(xb, sb)
        print("B", xb, sb, critb)

        x = attach(xa, xb)
        s = attach(sa, sb)

    return x, s
    



