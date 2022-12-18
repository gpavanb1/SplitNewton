import argparse
import logging
import numpy as np
from timeit import default_timer as time
from functions import set_functions
from newton import newton
from split_newton import split_newton

# Set logging level
parser = argparse.ArgumentParser()
parser.add_argument("--log", dest="loglevel", help="Set the loglevel for your solver  (DEBUG, INFO, WARNING, CRITICAL, ERROR)", type=str)
args = parser.parse_args()
loglevel = getattr(logging, args.loglevel.upper())
logging.basicConfig(level=loglevel)

# x0 = np.array([21.2, 16.3, 31.3, 9.2])
x0 = np.linspace(21.2, 31.2, 5000)

mode = "TEST"
func, der, hess = set_functions(mode)

dt0 = 0.0
dtmax = 0.1

# Newton
start = time()
print('Starting Newton...')
xf, _, iter = newton(der, hess, x0, sparse=True, dt0=dt0, dtmax=dtmax)
print("Final root: ", xf)
print("Final Residual: ", func(xf))
print(f"Elapsed time: {time() - start}")
print(f"Total iterations: {iter}")
input('')

print('-'* 20)

# Split Newton
start = time()
print('Starting Split-Newton...')
xf, _, iter = split_newton(der, hess, x0, int(len(x0)/2), sparse=True, dt0=dt0, dtmax=dtmax)
print("Final root: ", xf)
print("Final Residual: ", func(xf))
print(f"Elapsed time: {time() - start}")
print(f"Total iterations: {iter}")
input('')