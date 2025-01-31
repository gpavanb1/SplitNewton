import argparse
import logging
import numpy as np
from scipy import sparse as sp
from timeit import default_timer as time
from splitnewton.newton import newton
from splitnewton.split_newton import split_newton

# Create test functions
from scipy.optimize import rosen, rosen_der, rosen_hess
from examples.demo_func import test_func, test_der, test_hess


def set_functions(mode):
    if mode == "ROSENBROCK":
        return rosen, rosen_der, lambda x: sp.csr_matrix(rosen_hess(x))
    else:
        return test_func, test_der, test_hess


# Set logging level
parser = argparse.ArgumentParser()
parser.add_argument("--log", dest="loglevel",
                    help="Set the loglevel for your solver  (DEBUG, INFO, WARNING, CRITICAL, ERROR)", type=str, default="WARNING")
args = parser.parse_args()
loglevel = getattr(logging, args.loglevel.upper())
logging.basicConfig(level=loglevel)

# Seed
x0 = np.linspace(21.2, 31.2, 50)

mode = "TEST"  # or ROSENBROCK
func, der, hess = set_functions(mode)

dt0 = 0.0
dtmax = 0.1

# Split Newton
start = time()
print('Starting Split-Newton...')
xfa, _, iter = split_newton(
    der, hess, x0, [len(x0)//2], maxiter=100, sparse=True, dt0=dt0, dtmax=dtmax)
print("Final root: ", xfa)
print("Final Residual: ", func(xfa))
print(f"Elapsed time: {time() - start}")
print(f"Total iterations: {iter}")
input('')

print('-' * 20)

# Hierarchical Split Newton
start = time()
print('Starting Hierarchical Split-Newton...')
xfb, _, iter = split_newton(
    der, hess, x0, [len(x0)//4, 3*len(x0)//4, 7*len(x0)//8], maxiter=100, sparse=True, dt0=dt0, dtmax=dtmax)
print("Final root: ", xfb)
print("Final Residual: ", func(xfb))
print(f"Elapsed time: {time() - start}")
print(f"Total iterations: {iter}")
print("Split and Hierarchical agree/difference: ", "Yes" if np.allclose(
    xfa, xfb, atol=1e-6) else np.max(np.abs(xfa - xfb)))
input('')

print('-' * 20)


# Newton
start = time()
print('Starting Newton...')
xfc, _, iter = newton(der, hess, x0, sparse=True, dt0=dt0, dtmax=dtmax)
print("Final root: ", xfc)
print("Final Residual: ", func(xfc))
print(f"Elapsed time: {time() - start}")
print(f"Total iterations: {iter}")
print("Split and Newton agree/difference: ", "Yes" if np.allclose(
    xfa, xfc, atol=1e-6) else np.max(np.abs(xfa - xfc)))
input('')
