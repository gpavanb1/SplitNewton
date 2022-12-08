import numpy as np
from scipy.optimize import rosen, rosen_der, rosen_hess
from newton import newton

x0 = np.array([1.1, 1.4, 2.3, 1.2])
print('Starting Newton...')
xf = newton(rosen, rosen_der, rosen_hess, x0)
print(xf)
print(rosen(xf))
