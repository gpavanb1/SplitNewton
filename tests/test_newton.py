import pytest
import logging
import numpy as np
from scipy.sparse import csr_matrix
import numpy.linalg as la
from scipy.optimize import rosen, rosen_der, rosen_hess
from splitnewton.newton import criterion, check_within_bounds, newton

# Tests for criterion function


def test_criterion():
    x = np.array([1.0, 2.0])
    s = np.array([0.1, 0.2])
    result = criterion(x, s)
    expected = la.norm(s / (x * 1e-6 + 1e-5))
    assert np.isclose(result, expected)

# Tests for check_within_bounds function


def test_check_within_bounds_valid():
    x0 = [1.0, 2.0, 3.0]
    bounds = [[0.0, 1.5, 2.5], [2.0, 3.0, 4.0]]
    assert check_within_bounds(x0, bounds)


def test_check_within_bounds_invalid():
    x0 = [1.0, 2.0, 5.0]
    bounds = [[0.0, 1.5, 2.5], [2.0, 3.0, 4.0]]
    assert not check_within_bounds(x0, bounds)


def test_check_within_bounds_no_bounds():
    x0 = [1.0, 2.0]
    bounds = None
    assert check_within_bounds(x0, bounds)


def test_check_within_bounds_bounds_format():
    with pytest.raises(Exception, match="Bounds must be a 2 lists, with lower and upper bounds of each dimension"):
        check_within_bounds([1.0], [[0.0, 1.0]])

    with pytest.raises(Exception, match="Each bounds list must be as long as the solution vector"):
        check_within_bounds([1.0, 2.0], [[0.0], [1.0]])

# Tests for newton function


def test_negative_dt0():
    x0 = np.array([0.1, 0.2])
    func, der, hess = rosen, rosen_der, rosen_hess
    with pytest.raises(Exception, match="Must specify positive dt0 and dtmax"):
        newton(der, hess, x0, dt0=-1.)


def test_pseudotransient():
    x0 = np.array([4., 4.])
    expected_step = np.array([-2., -2.])
    # Function values and gradients

    def df(x):
        return x - 1

    def J(x):
        return np.eye(len(x))
    # Run the Newton solver with Armijo rule
    x, step, iterations = newton(
        df, J, x0, maxiter=1, sparse=True, dt0=2.)
    assert np.allclose(step, expected_step, atol=1e-5)


def test_invalid_seed():
    x0 = np.array([2., 2.])
    bounds = [[0., 0.], [1., 1.]]
    func, der, hess = rosen, rosen_der, rosen_hess
    with pytest.raises(Exception, match="Seed must be within the provided bounds"):
        newton(der, hess, x0, bounds=bounds)


def test_newton_no_bounds():
    x0 = np.array([0.1, 0.2])
    x_expected = np.array([1.0, 1.0])
    func, der, hess = rosen, rosen_der, rosen_hess
    x, s, iter = newton(der, hess, x0)
    assert np.allclose(x, x_expected, atol=1e-5)
    assert np.allclose(0.0, rosen(x), atol=1e-5)


def test_newton_with_bounds():
    x0 = np.array([0.25, 0.2])
    bounds = [[0., 0.], [2., 2.]]
    x_expected = np.array([1.0, 1.0])
    func, der, hess = rosen, rosen_der, rosen_hess
    x, s, iter = newton(der, hess, x0, bounds=bounds)
    assert np.allclose(x, x_expected, atol=1e-5)
    assert np.allclose(0.0, rosen(x), atol=1e-5)


def test_exact_bounds():
    x0 = np.array([0., 0.])
    bounds = [[0., 0.], [0., 0.]]
    func, der, hess = rosen, rosen_der, rosen_hess
    expected_step = [0., 0.]
    x, s, iter = newton(der, hess, x0, bounds=bounds, maxiter=1)
    assert np.allclose(s, expected_step, atol=1e-5)


def test_sparse_singular(caplog):
    x0 = np.array([0., 0.])
    bounds = [[0., 0.], [0., 0.]]
    # Function values and gradients
    EPS = 1e-16

    def df(x):
        return np.array([1.0, 1.0])  # Simple gradient

    def J(x):
        A = np.array([[1, 2], [2, 4 + EPS]])
        return csr_matrix(A)

    with caplog.at_level(logging.WARNING):
        x, s, iter = newton(df, J, x0, bounds=bounds, sparse=True)
    assert "GMRES not converged" in caplog.text


def test_armijo_rule():
    # Initial guess
    x0 = np.array([0.0])

    # Function values and gradients
    def df(x):
        return x - 1

    def J(x):
        return np.eye(len(x))

    # Expected result is close to the optimal point which is x = 1
    expected_x = np.array([1.0])
    tolerance = 1e-6

    # Run the Newton solver with Armijo rule
    x_opt, step, iterations = newton(
        df, J, x0, maxiter=1, sparse=True, armijo=True)

    # Check if Armijo scaling is working
    # The step should be scaled down, so check if the step is less than initial step size
    initial_step_size = 1.0
    assert np.linalg.norm(
        step) < initial_step_size, "Armijo scaling did not reduce the step size"


def test_newton_sparse_solver():
    x0 = np.array([0.1, 0.2])
    x_expected = np.array([1.0, 1.0])
    func, der, hess = rosen, rosen_der, rosen_hess
    x, s, iter = newton(der, hess, x0, sparse=True)
    assert np.allclose(x, x_expected, atol=1e-5)
    assert np.allclose(0.0, rosen(x), atol=1e-5)


def test_newton_invalid_bounds():
    x0 = np.array([0.0, 0.0])
    bounds = [[0.0, 0.0], [1.0]]
    func, der, hess = rosen, rosen_der, rosen_hess
    with pytest.raises(Exception, match="Each bounds list must be as long as the solution vector"):
        newton(der, hess, x0, maxiter=10, bounds=bounds)


# Run the tests
if __name__ == "__main__":
    pytest.main()
