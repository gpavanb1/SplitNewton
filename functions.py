from scipy.optimize import rosen, rosen_der, rosen_hess
from test_func import test_func, test_der, test_hess

def set_functions(mode):
    if mode == "ROSENBROCK":
        return rosen, rosen_der, rosen_hess
    else:
        return test_func, test_der, test_hess