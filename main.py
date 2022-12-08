import numpy as np
from scipy.optimize import rosen, rosen_der, rosen_hess
from newton import newton
from split_newton import split_newton

# Split Newton
x0 = np.array([1.1, 1.4, 2.3, 1.2])
print('Starting Split-Newton...')
xf = split_newton(rosen_der, rosen_hess, x0, 2)
print(xf)
print(rosen(xf))

print('-'* 20)

# Newton
x0 = np.array([1.1, 1.4, 2.3, 1.2])
print('Starting Newton...')
xf = newton(rosen_der, rosen_hess, x0)
print(xf)
print(rosen(xf))
