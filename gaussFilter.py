import numpy as np
import matplotlib.pyplot as plt


def gaussian(x, mu, sig):
    """Non normalized Gaussian Function.
    Args:
        x (float or np.array): the x-value(s)
        mu (float): center of gaussian
        sig (float): with of gaussian
    Returns:
        float or np.array: gaussian evaluated at x-value(s)"""
    return np.exp(-1/(2*sig**2)*(x-mu)**2)


def smooth_data(x, y, sig):
    """Smooth data with a gaussian with given with 'sig'.
    Args:
        x (np.array): the x values of the data
        y (np.array): the y values of the data
        sig (float):  the width of the gaussian
    Returns:
        np.array: the smoothed data at the x-values."""
    # clear data from np.nan's
    ind = np.where(np.isnan(y))
    x = np.delete(x, ind)
    y = np.delete(y, ind)
    sm_y = np.zeros(x.shape)
    for i in range(len(x)):
        gVal = gaussian(x, x[i], sig)
        gVal = gVal / np.sum(gVal)
        y = np.nan_to_num(y)
        ySm = gVal.dot(y)
        sm_y[i] = ySm
    return x, sm_y


if __name__ == "__main__":
    np.random.seed(5)
    n_points = 40
    x_vals = np.arange(n_points)
    y_vals = np.random.normal(size=n_points)
    plt.bar(x_vals, y_vals)
    plt.show()
    sig = 4 / np.sqrt(8 * np.log(2))
    y_pos = gaussian(x_vals, 13, sig)
    y_smooth = smooth_data(x_vals, y_vals, sig)
    plt.bar(x_vals, y_smooth)
    plt.show()
    plt.close()
