import numpy as np
import pytest
from splitnewton.newton import newton
from splitnewton.split_newton import split_newton

# Example gradient and Jacobian functions


def df_example(x):
    return 2 * (x - 1)


def J_example(x):
    return 2 * np.eye(len(x))

# Test comparison between newton and split_newton


def test_newton_vs_split_newton_sparse():
    x0 = np.array([0.5, 1.5, 0.5, 1.5], dtype=np.float64)
    bounds = [[0.]*4, [3.]*4]
    loc = [2]  # Split location

    # Run newton solver
    x_opt_newton, step_newton, iterations_newton = newton(
        df_example, J_example, x0, sparse=True, dt0=0.1, dtmax=1.0, armijo=True, bounds=bounds)

    # Run split_newton solver
    x_opt_split_newton, step_split_newton, iterations_split_newton = split_newton(
        df_example, J_example, x0, loc, sparse=True, dt0=0.1, dtmax=1.0, armijo=True, bounds=bounds)

    # Compare the results
    np.testing.assert_allclose(
        x_opt_newton, x_opt_split_newton, rtol=1e-6, atol=1e-4)
    np.testing.assert_allclose(
        step_newton, step_split_newton, rtol=1e-6, atol=1e-4)
    assert iterations_newton <= iterations_split_newton


def test_newton_vs_split_newton_dense():
    x0 = np.array([0.5, 1.5, 0.5, 1.5], dtype=np.float64)
    bounds = [[0.]*4, [3.]*4]
    loc = [2]  # Split location

    # Run newton solver
    x_opt_newton, step_newton, iterations_newton = newton(
        df_example, J_example, x0, sparse=False, dt0=0.1, dtmax=1.0, armijo=True, bounds=bounds)

    # Run split_newton solver
    x_opt_split_newton, step_split_newton, iterations_split_newton = split_newton(
        df_example, J_example, x0, loc, sparse=False, dt0=0.1, dtmax=1.0, armijo=True, bounds=bounds)

    # Compare the results
    np.testing.assert_allclose(
        x_opt_newton, x_opt_split_newton, rtol=1e-6, atol=1e-4)
    np.testing.assert_allclose(
        step_newton, step_split_newton, rtol=1e-6, atol=1e-4)
    assert iterations_newton <= iterations_split_newton


def test_negative_dt_exception():
    x0 = np.array([0.5, 1.5, 0.5, 1.5], dtype=np.float64)
    loc = [2]

    # Test negative dt0
    with pytest.raises(Exception, match="Must specify positive dt0 and dtmax"):
        split_newton(df_example, J_example, x0, loc, dt0=-0.1)

    # Test negative dtmax
    with pytest.raises(Exception, match="Must specify positive dt0 and dtmax"):
        split_newton(df_example, J_example, x0, loc, dtmax=-0.1)


def test_incorrect_split_location_exception():
    x0 = np.array([0.5, 1.5, 0.5, 1.5], dtype=np.float64)
    loc = [5]  # Incorrect location, greater than length of x0

    with pytest.raises(Exception, match="Incorrect split location"):
        split_newton(df_example, J_example, x0, loc)


if __name__ == "__main__":
    pytest.main([__file__])
