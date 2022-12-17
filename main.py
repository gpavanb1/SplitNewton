import numpy as np
from functions import set_functions
from newton import newton
from split_newton import split_newton

x0 = np.array([21.2, 16.3, 31.3, 9.2])

mode = "ROSENBROCK"
func, der, hess = set_functions(mode)


# Newton
print('Starting Newton...')
xf, _ = newton(der, hess, x0)
print(xf)
print(func(xf))
input('')

print('-'* 20)

# Split Newton
print('Starting Split-Newton...')
xf, _ = split_newton(der, hess, x0, 2)
print(xf)
print(func(xf))
input('')