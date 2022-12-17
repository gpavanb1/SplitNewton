import numpy as np

lambda_a = 4
lambda_b = 2
lambda_c = -1
lambda_d = -2

def common_init(x0):
    la = int(len(x0)/2)
    lb = len(x0) - la

    a = np.logspace(lambda_a, lambda_b, la)
    b = np.logspace(lambda_c, lambda_d, lb)

    return a, b

def test_func(x0):
    a, b = common_init(x0)
    coeff = np.concatenate((a, b))
    return 0.25 * sum(coeff * (np.array(x0) ** 4))

def test_der(x0):
    a, b = common_init(x0)
    return np.concatenate((a, b)) * (np.array(x0) ** 3)

def test_hess(x0):
    a, b = common_init(x0)
    diagonal = 3 * np.concatenate((a, b)) * (np.array(x0) ** 2)
    return np.diag(diagonal)

