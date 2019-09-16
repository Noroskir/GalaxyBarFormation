import numpy as np
import matplotlib.pyplot as plt
import galaxy as g


def simulate_errors(x, y, yerr, N=1000):
    """Monte Carlo simulation to estimate the errors.
    Args:
        x (np.array): array of the x values
        y (np.array): array of the y values
        yerr (np.array): array of the errors of the y values
        N (int): number of curves
    Returns: 
        np.array: array of arrays with simulated curves."""
    mc_y = []
    for i in range(N):
        yNew = np.array([np.random.normal(y[j], yerr[j], 1)[0]
                         for j in range(len(x))])
        mc_y.append(yNew)

    return np.array(mc_y)


def simulate_errors_fast(x, y, yerr, N=1000):
    """Monte Carlo simulation to estimate the errors.
    Args:
        x (np.array): array of the x values
        y (np.array): array of the y values
        yerr (np.array): array of the errors of the y values
        N (int): number of curves
    Returns: 
        np.array: array of arrays with simulated curves."""
    mc_y = []
    for i in range(len(y)):
        yNew = np.random.normal(y[i], yerr[i], N)
        mc_y.append(yNew)
    return np.transpose(np.array(mc_y))


def mean_2d_array(x):
    """Calculate the mean value and the standard deviation of the columns of a 2d np.array.
    Args:
        x (np.array): array of arrays of the values
    Returns:
        tuple (np.array, np.array): mean value and standard deviation."""
    N = len(x)
    mean = []
    std = []
    for i in range(len(x[0])):
        m = np.sum(x[:, i]/N)
        s = np.sqrt(1 / (N-1) * np.sum((x[:, i] - m)**2))
        mean.append(m)
        std.append(s)
    return np.array(mean), np.array(std)
