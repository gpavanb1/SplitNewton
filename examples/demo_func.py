import numpy as np

lambda_a = 6
lambda_b = 2
lambda_c = -1
lambda_d = -4
lambda_e = -6
lambda_f = -8


def common_init(x0):
    la = len(x0)//3
    lb = len(x0)//3
    lc = len(x0) - la - lb

    a = np.logspace(lambda_a, lambda_b, la)
    b = np.logspace(lambda_c, lambda_d, lb)
    c = np.logspace(lambda_c, lambda_d, lc)

    return a, b, c


def test_func(x0):
    a, b, c = common_init(x0)
    coeff = np.concatenate((a, b, c))
    return 0.25 * sum(coeff * (np.array(x0) ** 4))


def test_der(x0):
    a, b, c = common_init(x0)
    return np.concatenate((a, b, c)) * (np.array(x0) ** 3)


def test_hess(x0):
    a, b, c = common_init(x0)
    diagonal = 3 * np.concatenate((a, b, c)) * (np.array(x0) ** 2)
    return np.diag(diagonal)
