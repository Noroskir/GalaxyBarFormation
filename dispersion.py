import numpy as np
import matplotlib.pyplot as plt
import galaxy as gal

def read_spherical(filename):
    '''Read the columns of the data-file in the format:
    Radius[arcsec], betaR, err, vrr, error, vpp, error,  vtt, error, vp, error
    Args:
        filename (str): the name of the file
    Returns:
        tuple of lists:  '''
    radius = []
    betaR, e_betaR = [], []
    vrr, e_vrr = [], []
    vpp, e_vpp = [], []
    vtt, e_vtt = [], []
    vp, e_vp = [], []
    file = open(filename, 'r')
    lines = file.readlines()
    for line in lines:
        line = line.split()
        if len(line) == 11:
            try:
                radius.append(float(line[0]))
                betaR.append(float(line[1]))
                e_betaR.append(float(line[2]))
                vrr.append(float(line[3]))
                e_vrr.append(float(line[4]))
                vpp.append(float(line[5]))
                e_vpp.append(float(line[6]))
                vtt.append(float(line[7]))
                e_vtt.append(float(line[8]))
                vp.append(float(line[9]))
                e_vp.append(float(line[10]))
            except:
                print("Error: Converting string into int:")
                print(line)

    return radius, betaR, e_betaR, vrr, e_vrr, vpp, e_vpp, vtt, e_vtt, vp, e_vp


def read_cylindrical(filename):
    '''Read the columns of the data-file in the format:
    Radius, betaZ, e_betaZ, vzz, error, vRR, error, vpp, error, vp, error
    Args:
        filename (str): name of the file
    Returns:
        tuple of lists
    '''
    radius = []
    betaZ, e_betaZ = [], []
    vzz, e_vzz = [], []
    vrr, e_vrr = [], []
    vpp, e_vpp = [], []
    vp, e_vp = [], []
    file = open(filename, 'r')
    lines = file.readlines()
    for line in lines:
        line = line.split()
        if len(line) == 11:
            try:
                radius.append(float(line[0]))
                betaZ.append(float(line[1]))
                e_betaZ.append(float(line[2]))
                vzz.append(float(line[3]))
                e_vzz.append(float(line[4]))
                vrr.append(float(line[5]))
                e_vrr.append(float(line[6]))
                vpp.append(float(line[7]))
                e_vpp.append(float(line[8]))
                vp.append(float(line[9]))
                e_vp.append(float(line[10]))
            except:
                print("Error: Converting string into float:")
                print(line)

    return radius, betaZ, e_betaZ, vzz, e_vzz, vrr, e_vrr, vpp, e_vpp, vp, e_vp


def velocity_disp_sph(radius, betaR, e_betaR, vrr, e_vrr, vpp, e_vpp, vtt, e_vtt, vp, e_vp):
    """Calculate the fraction between the velocity dispersion vrr and vtt.
    Args:
        lists
    Returns:
        tuple (radius, vrr/vtt): tuple of two np.arrays, with x and y values"""
    vpp = np.array(vpp)
    vtt = np.array(vtt)
    vrr = np.array(vrr)
    vp = np.array(vp)
    radius = np.array(radius)
    sig_t = np.sqrt((vpp - vp**2 + vtt)/2)
    sig_r = np.sqrt(vrr)
    print("Spherical:\n vpp:")
    print(vpp, "\n vtt:")
    print(vtt, "\n vrr:")
    print(vrr, "\n Radius:")
    print(radius)

    return radius, sig_r/sig_t


def velocity_disp_cyl(radius, betaZ, e_betaZ, vzz, e_vzz, vRR, e_vRR, vpp, e_vpp, vp, e_vp):
    """Calculate the fraction between the velocity dispersion vzz and vRR.
    Args:
        lists
    Returns:
        tuple (radius, vzz/vRR): tuple of np.arrays, with x and y values"""
    sig_z = np.sqrt(np.array(vzz))
    sig_R = np.sqrt(np.array(vRR))
    sig_z = sig_z[sig_z != 0]
    sig_R = sig_R[sig_R != 0]
    radius = radius[:len(sig_R)]

    s_sig_z = 0.5 * 1 / \
        np.sqrt(np.array(vzz)[:len(sig_z)]) * np.array(e_vzz)[:len(sig_z)]
    s_sig_R = 0.5 * 1 / \
        np.sqrt(np.array(vRR)[:len(sig_R)]) * np.array(e_vRR)[:len(sig_R)]
    radius = np.array(radius)

    error = sig_z/sig_R * np.sqrt((s_sig_z/sig_z)**2 + (s_sig_R/sig_R)**2)

    # print("Cylindrical:\n vzz:")
    # print(vzz, "\n vRR:")
    # print(vRR, "\n Radius")
    # print(radius)
    return radius, sig_z/sig_R, error


def plot(x1, y1, e1, x2, y2, e2, title):
    plt.plot(x1, y1, '-x', color='black', label=r'$\sigma_r/\sigma_t (r)$')
    plt.plot(x2, y2, '-x', color='red', label=r'$\sigma_z/\sigma_R (R)$')
    plt.errorbar(x2, y2, yerr=e2, fmt='none', color='red', capsize=2)
    plt.errorbar(x1, y1, yerr=e1, fmt='none', color='red', capsize=2)
    plt.ylim(0, 2)
    plt.xlim(0, 40)
    plt.grid('--')
    plt.title(title)
    plt.legend()
    plt.show()


def derivative(x, y):
    """Differentiate a function with discrete values.
    Args:
        x (np.array): the x-values
        y (np.array): the y-values
    Returns:
        tuple (x, y): np.arrays"""
    xm = x[:-1]
    xp = x[1:]
    ym = y[:-1]
    yp = y[1:]
    Dy = (yp - ym)/(xp-xm)
    Dx = xm + (xp-xm)/2
    dx = xp-xm
    return Dx, Dy, dx


def test_derivative():
    x = np.linspace(0, 10, 20)
    y = np.sin(x)
    yprime = np.cos(x)
    Dx, Dy, dx = derivative(x, y)

    x2 = Dx + dx/2
    y2 = (Dy[1:] + Dy[:-1])/2

    y3 = np.gradient(y, x)
    
    #plt.plot(x, y, '-x', label='f(x)')
    plt.plot(Dx, Dy, '*', label='num derivative')
    plt.plot(x2[:-1], y2, '-s', label='shifted Num deri')
    plt.plot(x, yprime, 'o', label='Analytical deri')
    plt.plot(x, y3, '^', label='np gradient')
    plt.legend()
    plt.show()

def kappa_sq(radius, vp):
    R = np.array(radius)
    R = R[R != 0]
    vc = np.array(vp)[:len(R)]
    omega = vc/R
    Domega = np.gradient(omega, R)

    # Plot the derivative of Omega(R)
    # plt.figure(figsize=(10, 6))
    # plt.plot(R, omega, '-x', label=r'$\Omega (R)$')
    # plt.plot(R, Domega, '-x', label=r'$d\Omega(R)/dR$')
    # # plt.plot(DR, Domega_my, '-x', label='My derive')
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    kappa_squared = 4*omega**2 + 2*R* omega*Domega
    print(kappa_squared)
    plt.plot(R, kappa_squared)
    plt.grid()
    plt.show()

    
    

if __name__ == "__main__":
    galaxy = "NGC6278"
    #data1 = read_spherical("dapta/"+galaxy+"_anisotropy.dat")
    #  velocity_disp(*data)
    data2 = read_cylindrical("data/"+galaxy+"_betaz.dat")

    g4210 = gal.Galaxy('data/', 'NGC6278')

    # x1, y1 = velocity_disp_sph(*data2)
    x2, y2, e2 = velocity_disp_cyl(*data2)
    x1 ,y1, e1 = g4210.velocity_disp_cyl()
    #plot(x1, y1, e1, x2, y2, e2, galaxy)


    
    kappa_sq(data2[0], data2[9])

    #test_derivative()
    
    # DERIVATIVE TESTING
    # x = np.linspace(0, 10, 100)
    # y = np.sin(x)
    # Dx, Dy = derivative(x, y)
    # plt.plot(x, y, 'o', label='sin(x)')
    # plt.plot(Dx, Dy, label='Derivative')
    # plt.plot(x, np.cos(x), label='cos(x)')
    # plt.legend()
    # plt.show()
