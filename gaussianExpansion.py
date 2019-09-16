import numpy as np
import matplotlib.pyplot as plt


def gaussian(x,y, L, sig, q):
    return L *  np.exp(-(x**2+ y**2/q**2)/(2*sig**2))

def get_gaussian(L, sig, q):
    """Create a 2D gaussian with fixed parameters.
    Args:
        L (float): Luminosity
        sig (float): with of the gaussian
        q (float): intrinsic flattening
    Returns:
        function g(x, y)."""
    # return lambda x,y: L /(2*np.pi*sig**2 *q) * np.exp(-(x**2+ y**2/q**2)/(2*sig**2))
    return lambda x,y: L * np.exp(-(x**2 + y**2/q**2)/(2*sig**2))


def read_mass(filename):
    """Read the galaxy_mass_mge.txt file.
    Args:
        filename (str): name of the file
    Returns:
        None."""
    f = open(filename, 'r')
    lines = f.readlines()
    line = lines[0].split()
    massLight = float(line[0])
    line = lines[1].split()
    inclAngle = float(line[0])
    line = lines[2].split()
    nGaussLM = int(line[0])
    line = lines[3].split()
    nGaussDM = int(line[0])
    paramGaussLM = []
    paramGaussDM = []
    for i in range(6, 6 + nGaussLM):
        l = lines[i].split()
        paramGaussLM.append((float(l[0]), float(l[1]), float(l[2])) )
    for i in range(6 + nGaussLM + 2, 6 + nGaussLM + 2 + nGaussDM):
        l = lines[i].split()
        paramGaussDM.append((float(l[0]), float(l[1]), float(l[2])))

    return massLight, inclAngle, paramGaussLM, paramGaussDM

def multi_gauss_exp(paramGauss):
    """Creates a functions which consists of the sum of 2D Gaussians with the
    parameters given.
    Args:
        paramGaussLM (list): list with tupels of the parameters, (L, sigma, q)
    Returns:
        function mge(x, y): Multi-Gaussian-Expansion with x and y coordinates."""
    gaussians = [ get_gaussian(*p) for p in paramGauss]
    return lambda x, y: sum(g(x,y) for g in gaussians)
    

if __name__=="__main__":
    galaxy = "NGC4210"
    massLight, inclAngle, paramGaussLM, paramGaussDM = read_mass('data/'+galaxy+'_mass_mge.txt')
    r = np.linspace(0, 10, 500)
    for p in paramGaussLM:
        plt.plot(r, gaussian(r, np.zeros(len(r)), *p), label='part')
    mge = multi_gauss_exp(paramGaussLM)
    plt.semilogy(r, mge(r, np.zeros(len(r))), label='Sum')
    plt.ylim(1, 1e4)
    plt.xlim(0, 50)
    plt.legend()
    plt.show()
