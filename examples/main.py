import argparse
import logging
import numpy as np
from timeit import default_timer as time
from splitnewton.newton import newton
from splitnewton.split_newton import split_newton

# Create test functions
from scipy.optimize import rosen, rosen_der, rosen_hess
from examples.demo_func import test_func, test_der, test_hess


def set_functions(mode):
    if mode == "ROSENBROCK":
        return rosen, rosen_der, rosen_hess
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
x0 = np.linspace(21.2, 31.2, 5000)

mode = "TEST"  # or ROSENBROCK
func, der, hess = set_functions(mode)

dt0 = 0.0
dtmax = 0.1

# Split Newton
start = time()
print('Starting Split-Newton...')
xf, _, iter = split_newton(der, hess, x0, int(
    len(x0)/2), sparse=True, dt0=dt0, dtmax=dtmax)
print("Final root: ", xf)
print("Final Residual: ", func(xf))
print(f"Elapsed time: {time() - start}")
print(f"Total iterations: {iter}")
input('')

print('-' * 20)

# Newton
start = time()
print('Starting Newton...')
xf, _, iter = newton(der, hess, x0, sparse=True, dt0=dt0, dtmax=dtmax)
print("Final root: ", xf)
print("Final Residual: ", func(xf))
print(f"Elapsed time: {time() - start}")
print(f"Total iterations: {iter}")
input('')
