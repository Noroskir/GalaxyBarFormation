import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import tqdm
import subprocess
import galaxy as g
import photometrics as ph


P = ph.Photometrics("data/photometric_decomposition.fits")


def read_master_file(filename):
    """Read the CALIFA master file.
    Args:
        filename (str): location and name of file
    Returns:
        list: of names"""
    f = open(filename, 'r')
    lines = f.readlines()
    names = list()
    hType = list()
    for l in lines:
        l = l.split()
        if l[0][0] != '#':
            names.append(l[1])
            hType.append(l[10])
    return names


def read_filenames(directory, ending):
    """Read the galaxy names of a directory with a specified file ending.
    Args:
        directory (str): name of directory
        ending (str): common part of the filenames
    Returns:
        list: list with filenames."""
    ls = subprocess.Popen(["ls", "-p", directory],
                          stdout=subprocess.PIPE).stdout
    names = list()
    for line in ls:
        name = str(line)[2:-len("_mass_mge.txt\n'") - 1]
        names.append(name)

    return names


def analyze_galaxies(names):
    """Calculate the Toomre Parameter for the selected galaxies.
    Generate plots and write into files.
    Args:
        location (str): location of the files 'betaZ' and 'massmge'
    Returns:
        None."""
    length = len(names)
    with tqdm.tqdm(total=length, desc="Galaxies") as pbar:
        for i in range(length):
            try:
                pbar.write("%s" % names[i])
                galaxy = g.Galaxy("data/sRsz_data/", "data/massmge/", names[i])
                galaxy.plot_total_cyl()
            except Exception as E:
                tqdm.tqdm.write('Exception', E)
            pbar.update(1)


def read_toomre(directory, names):
    """Read the files with the Toomre parameter and return it and the radius.
    Args:
        directory (str): location of files
        names (list): list of galaxy names
    Returns:
        np.array: np.arrays of R, Q, e_Q"""
    toomre = []
    e_toomre = []
    radius = []  # arcsec
    for n in names:
        f = open(directory + n + '_data.txt', 'r')
        lines = f.readlines()
        R = []
        Q = []
        e_Q = []
        for l in lines:
            l = l.split()
            if len(l) > 4:
                try:
                    R.append(float(l[0]))
                    Q.append(float(l[2]))
                    if float(l[2]) > 200:
                        print("Galaxie:")
                        print(n)
                    e_Q.append(float(l[3]))
                except:
                    pass
        toomre.append(np.array(Q))
        e_toomre.append(np.array(e_Q))
        radius.append(np.array(R))
    R = np.array(radius)  # arcsec
    Q = np.array(toomre)
    e_Q = np.array(e_toomre)
    return R, Q, e_Q


def read_vcirc_toomre(directory, names):
    """Read the files with the Toomre parameter calculated with the vcirc
    velocity profile and return it and the radius.
    Args:
        directory (str): location of files
        names (list): list of galaxy names
    Returns:
        np.array: np.arrays of R, Q, e_Q"""
    toomre = []
    radius = []  # arcsec
    for n in names:
        bNan = False
        f = open(directory + n + '_data.txt', 'r')
        lines = f.readlines()
        R = []
        sigma_R = []
        kappa_vcirc = []
        Sigma = []
        dist = float(lines[1].split()[0])
        for l in lines:
            l = l.split()
            if len(l) > 4:
                try:
                    R.append(float(l[0]))
                    sigma_R.append(float(l[7]))
                    kappa_vcirc.append(float(l[-1]))
                    if np.isnan(kappa_vcirc[-1]):
                        bNan = True
                        print(n)
                        break
                    Sigma.append(float(l[6]))
                except:
                    pass
        if not bNan:
            R = np.array(R)
            sigma_R = np.array(sigma_R)
            Sigma = np.array(Sigma) * const.M_sun.value / const.pc.value**2
            kappa = np.array(kappa_vcirc)
            kappa = kappa / (dist * const.kpc.value *
                             np.pi / (3600*180))  # 1/km
            Q = kappa * sigma_R * 1e3 / (3.36 * Sigma * const.G.value)
            radius.append(R)
            toomre.append(Q)
    R = np.array(radius)  # arcsec
    Q = np.array(toomre)
    return R, Q


def read_swing_ampl(names):
    """Calculate the swing amplification parameter X for the galaxies.
    Args:
        names (list): list of galaxies
    Returns:
        np.array, np.array: R and X.
    """
    X = []
    R = []
    for i in range(len(names)):
        f = open("data/swing/"+names[i]+"_X_data.txt", 'r')
        lines = f.readlines()
        lR, lX = [], []
        for l in lines[2:]:
            l = l.split()
            lR.append(float(l[0]))
            lX.append(float(l[1]))
        R.append(np.array(lR))
        X.append(np.array(lX))
    return R, X


def get_effectiveR(names):
    """Get effective radii of galaxies in [arcsec].
    Args:
        names (list): list of galaxie names.
    Returns:
        list: list with effective radii."""
    f = open("data/CALIFA_master_selection.txt", 'r')
    lines = f.readlines()
    galaxy = list()
    reff = list()
    for l in lines:
        l = l.split()
        if l[0][0] != '#':
            galaxy.append(l[1])
            reff.append(float(l[-4]))
    d = {galaxy[i]: reff[i] for i in range(len(galaxy))}

    return [d[n] for n in names]


def get_median_mean_std(x):
    """Calculate the median, mean and the standard deviation for every row
    of a two dimensional array x.
    Args:
        x (np.array): 2D np.array
    Returns:
        (np.array, np.array, np.array).
    """
    median = np.zeros(len(x))
    mean = np.zeros(len(x))
    std = np.zeros(len(x))
    for i in range(len(x)):
        median[i] = np.median(x[i])
        mean[i] = np.mean(x[i])
        std[i] = np.std(x[i], ddof=1)
    return median, mean, std


def get_sigma_Z(names, R=None):
    """Get the velocity dispersion in z direction.
    Args:
        names (list): galaxies
        R (np.array or None): if not none restrict to radial values
    Returns:
        np.array, np.array: sigma_Z, e_sigma_Z"""
    sigZ = []
    e_sigZ = []
    for i in range(len(names)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", names[i])
        if R == None:
            s = np.sqrt(G.data_cyl['vzz'])
            e_s = 0.5/s * G.data_cyl['e_vzz']
            sigZ.append(s)
            e_sigZ.append(e_s)
        else:
            r = G.data_cyl['R']
            s = np.sqrt(G.data_cyl['vzz'])
            e_s = 0.5/s * G.data_cyl['e_vzz']
            SZ = []
            e_SZ = []
            for j in range(len(r)):
                if r[j] in R[i]:
                    SZ.append(s[j])
                    e_SZ.append(e_s[j])
            sigZ.append(np.array(SZ))
            e_sigZ.append(e_SZ)
    return sigZ, e_sigZ


def get_type(names):
    """Get Hubble types of galaxies.
    Args:
        names (list): list of galaxie names.
    Returns:
        list: list with Hubble types."""
    f = open("data/CALIFA_master_selection.txt", 'r')
    lines = f.readlines()
    galaxy = list()
    hType = list()
    for l in lines:
        l = l.split()
        if l[0][0] != '#':
            galaxy.append(l[1])
            hType.append(l[10])
    d = {galaxy[i]: hType[i] for i in range(len(galaxy))}

    return [d[n] for n in names]


def barred_unbarred(names, title='', vcirc=False):
    """Compare barred and unbarred galaxies in their Toomre Parameter
    Args:
        None.
    Returns:
        None."""
    barred = P.is_barred(names)
    disked = P.is_disked(names)
    dGalax = []
    bGalax = []
    for i in range(len(disked)):
        if disked[i] and not barred[i]:
            dGalax.append(names[i])
        elif barred[i]:
            bGalax.append(names[i])
    # from data
    if not vcirc:
        Rd, Qd, e_Qd = read_toomre("data/toomre/", dGalax)
        Rb, Qb, e_Qb = read_toomre("data/toomre/", bGalax)
    else:
        # from vcirc modelling
        Rd, Qd = read_vcirc_toomre("data/toomre/", dGalax)
        Rb, Qb = read_vcirc_toomre("data/toomre/", bGalax)
    Rdeff = np.array(get_effectiveR(dGalax))
    Rbeff = np.array(get_effectiveR(bGalax))
    Rd = Rd / Rdeff
    Rb = Rb / Rbeff
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].grid()
    axes[1].grid()
    for i in range(len(Rd)):
        if i == 0:
            axes[0].plot(Rd[i], Qd[i], 'x', color='black', label='Disk')
        else:
            axes[0].plot(Rd[i], Qd[i], 'x', color='black')
    for i in range(len(Rb)):
        if i == 0:
            axes[1].plot(Rb[i], Qb[i], 'x', color='red', label='Bar')
        else:
            axes[1].plot(Rb[i], Qb[i], 'x', color='red')
    axes[0].set_title('No Bar')
    axes[0].set_ylim(0, 3)
    axes[0].set_xlim(0, 8)
    axes[0].set_xlabel(r'R in $[R_e]$')
    axes[0].set_ylabel('Q')
    axes[1].set_title('With Bar')
    axes[1].set_ylim(0, 3)
    axes[1].set_xlim(0, 8)
    axes[1].set_xlabel(r'R in $[R_e]$')
    axes[1].set_ylabel('Q')
    if title != '':
        fig.suptitle(title)
    plt.show()
    plt.close()


def compare_types(names, hType, Rcut=False, title='', vcirc=False):
    """Compare two types of galaxies in their Toomre parameter
    Args:
        names (list): names of galaxies
        hType (str): Hubble type
    Returns:
        None."""
    lTypes = get_type(names)
    l1 = []
    for i in range(len(lTypes)):
        if lTypes[i] == hType:
            l1.append(names[i])
    if Rcut:
        compare_dominated_regions(l1, title, vcirc)
    else:
        barred_unbarred(l1, title, vcirc)


def get_disk_galaxies(names):
    """Filter out the disk galaxies without a bar.
    Args:
        names (list): galaxy names
    Returns:
        list: list of disked galaxies"""
    disked = P.is_disked(names)
    barred = P.is_barred(names)
    dGalax = []
    for i in range(len(disked)):
        if disked[i] and not barred[i]:
            dGalax.append(names[i])
    return dGalax


def get_barred_galaxies(names):
    """Filter out the galaxies with a bar.
    Args:
        names (list): galaxy names
    Returns:
        list: list of barred galaxies.
    """
    barred = P.is_barred(names)
    bGalax = []
    for i in range(len(barred)):
        if barred[i]:
            bGalax.append(names[i])
    return bGalax


def get_elliptical_galaxies(names):
    """Filter out the galaxies with a single sersic component.
    Args:
        names (list): galaxy names
    Returns:
        list: list of elliptical galaxies.
    """
    elli = []
    # hTypes = get_type(names)
    # for i in range(len(hTypes)):
    #     if hTypes[i][0] == "E":
    #         elli.append(names[i])
    for i in range(len(names)):
        if P.is_bulged(names[i]) and not P.is_barred(names[i]) and not P.is_disked(names[i]):
            elli.append(names[i])
    return elli


def stack_curves(x, y, step, threshold=1):
    """Stack curves in an interval with a certain stepsize.
    Args:
        x (np.array): 2D array with the x values
        y (np.array): 2D array with the y values
        step (float): size of interval in which the y values will be stacked
        threshold (int): minimal number of entries in bin to give an output value
    Returns:
        (x, y): np.arrays x-centers of stacked values, y-stacked values.
    """
    if len(x) == 0:
        return np.array([]), np.array([])
    s = step
    bins = [[]]
    valmax = max([np.max(x[i]) for i in range(len(x))])
    while s < valmax:
        s += step
        bins.append([])
    for i in range(len(x)):
        for j in range(len(x[i])):
            if y[i][j] != 0.:
                ind = int(x[i][j] / step)
                bins[ind].append(y[i][j])
    stacked = []
    xnew = []
    for i in range(len(bins)):
        if len(bins[i]) != 0 and len(bins[i]) >= threshold:
            stacked.append(sum(bins[i])/len(bins[i]))
            xnew.append(i*step + step/2)
        else:
            pass
    return np.array(xnew), np.array(stacked)


def mass_bin(R, Q, M, mbins):
    """Bin Q with respect to the mass.
    Args:
        R (np.array): 2D array - radius values
        Q (np.array): 2D array - Q values
        M (np.array): 1D array - Mass values
        mbins (list): list with the lowest to highest mass value
    Returns:
        R (list), Q (list), M (list).
    """
    Rnew = [[] for i in range(len(mbins)-1)]
    Qnew = [[] for i in range(len(mbins)-1)]
    step = mbins[1] - mbins[0]
    for i in range(len(R)):
        for j in range(len(mbins)):
            if M[i] == mbins[0]:
                Rnew[0].append(R[i])
                Qnew[0].append(Q[i])
                break
            elif M[i] <= mbins[0] + j*step:
                Rnew[j-1].append(R[i])
                Qnew[j-1].append(Q[i])
                break
    return Rnew, Qnew


def compare_dominated_regions(names, title='', vcirc=False):
    """Compare the disk and bar dominated regions of a galaxy type
    in their Toomre parameter.
    Args:
        names (list): names of galaxies
        hType (str): Hubble type
    Returns:
        None.
    """
    barred = P.is_barred(names)
    disked = P.is_disked(names)
    dGalax = []
    bGalax = []
    for i in range(len(disked)):
        if disked[i] and not barred[i]:
            dGalax.append(names[i])
        elif barred[i]:
            bGalax.append(names[i])
    if not vcirc:
        Rd, Qd, e_Qd = read_toomre("data/toomre/", dGalax)
        Rb, Qb, e_Qb = read_toomre("data/toomre/", bGalax)
    else:
        Rd, Qd = read_vcirc_toomre("data/toomre/", dGalax)
        Rb, Qb = read_vcirc_toomre("data/toomre/", bGalax)
    Rdeff = np.array(get_effectiveR(dGalax))
    Rbeff = np.array(get_effectiveR(bGalax))
    Rd, Qd = filter_dom_regions(dGalax, Rd, Qd, filt='disk')
    Rb, Qb = filter_dom_regions(bGalax, Rb, Qb, filt='bar')
    Rd = Rd / Rdeff
    Rb = Rb / Rbeff
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].grid()
    axes[1].grid()
    for i in range(len(Rd)):
        if i == 0:
            axes[0].plot(Rd[i], Qd[i], 'x', color='black', label='Disk')
        else:
            axes[0].plot(Rd[i], Qd[i], 'x', color='black')
    for i in range(len(Rb)):
        if i == 0:
            axes[1].plot(Rb[i], Qb[i], 'x', color='red', label='Bar')
        else:
            axes[1].plot(Rb[i], Qb[i], 'x', color='red')
    axes[0].set_title('No Bar')
    axes[0].set_ylim(0, 3)
    axes[0].set_xlim(0, 8)
    axes[0].set_xlabel(r'R in $[R_e]$')
    axes[0].set_ylabel('Q')
    axes[1].set_title('With Bar')
    axes[1].set_ylim(0, 3)
    axes[1].set_xlim(0, 8)
    axes[1].set_xlabel(r'R in $[R_e]$')
    axes[1].set_ylabel('Q')
    if title != '':
        fig.suptitle(title)
    plt.show()
    plt.close()


def toomre_ellipticals(names):
    """Plot the toomre parameter for elliptical galaxies.
    Args:
        names (list): names of the galaxies
    Returns:
        None.
    """
    E0, E1, E2, E3, E4, E5, E6, E7 = [], [], [], [], [], [], [], []
    hTypes = get_type(names)
    for i in range(len(hTypes)):
        if hTypes[i] == "E0":
            E0.append(names[i])
        elif hTypes[i] == "E1":
            E1.append(names[i])
        elif hTypes[i] == "E2":
            E2.append(names[i])
        elif hTypes[i] == "E3":
            E3.append(names[i])
        elif hTypes[i] == "E4":
            E4.append(names[i])
        elif hTypes[i] == "E5":
            E5.append(names[i])
        elif hTypes[i] == "E6":
            E6.append(names[i])
        elif hTypes[i] == "E7":
            E7.append(names[i])

    E = [E0, E1, E2, E3, E4, E5, E6, E7]
    color = ['yellow', 'green', 'blue', 'red',
             'black', 'lightblue', 'brown', 'orange']
    labels = ['E0', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7']
    for i in range(len(E)):
        print(len(E[i]))
        R, Q = read_vcirc_toomre("data/toomre/", E[i])
        # R, Q, e_Q = read_toomre("data/toomre/", E[i])
        Reff = np.array(get_effectiveR(E[i]))
        R = R / Reff
        for j in range(len(E[i])):
            if j == 0:
                plt.plot(R[j], Q[j], 'x', color=color[i], label=labels[i])
            else:
                plt.plot(R[j], Q[j], 'x', color=color[i])
    plt.grid()
    plt.legend()
    plt.show()


def read_mass(names):
    """Get the total galaxy mass.
    Args:
        names (list): names of galaxies
    Returns:
        list: list of galaxy masses."""
    mass = []
    for n in names:
        f = open("data/toomre/"+n+"_data.txt")
        lines = f.readlines()
        m = float(lines[2].split()[0])
        mass.append(m)
    return mass


def read_stellar_mass(names):
    """Get the stellar mass from the 'stellar_masses.txt' file.
    Args:
        names (list): names of the galaxies
    Returns:
        np.array: array with the masses.
    """
    f = open("data/stellar_masses.txt", 'r')
    lines = f.readlines()
    massFile = []
    nameFile = []
    masses = np.zeros(len(names))
    for l in lines:
        l = l.split()
        if l[0] != '#':
            nameFile.append(l[1])
            massFile.append(float(l[3]))
    for i in range(len(names)):
        ind = nameFile.index(names[i])
        masses[i] = massFile[ind]
    return masses


def read_lambda_re(names):
    """Get the lambda_re parameter from the 'stellar_masses.txt' files.
    Args:
        names (list): names of the galaxies
    Returns:
        np.array: array with the lambda_re parameter.
    """
    f = open("data/stellar_masses.txt", 'r')
    lines = f.readlines()
    lambFile = []
    nameFile = []
    lambRe = np.zeros(len(names))
    for l in lines:
        l = l.split()
        if l[0] != '#':
            nameFile.append(l[1])
            lambFile.append(float(l[2]))
    for i in range(len(names)):
        ind = nameFile.index(names[i])
        lambRe[i] = lambFile[ind]
    return lambRe


def mass_overview(names):
    mass = read_mass(names)
    for i in range(len(mass)):
        plt.plot(i, mass[i], 'o')
    plt.grid()
    plt.show()
    plt.close()


def compare_by_mass(names, m_low, m_upp, vcirc=False, title=''):
    """Compare barred and unbarred galaxies by their mass.
    Args:
        names (list): names of galaxies
        m_low (float): lower bound of galaxy mass
        m_upp (float): upper bound of galaxy mass
    Returns:
        None.
    """
    barred = P.is_barred(names)
    disked = P.is_disked(names)
    dGalax = []
    bGalax = []
    for i in range(len(disked)):
        if disked[i] and not barred[i]:
            dGalax.append(names[i])
        elif barred[i]:
            bGalax.append(names[i])
    Md = read_mass(dGalax)
    Mb = read_mass(bGalax)
    lbar = []
    ldisk = []
    for i in range(len(Md)):
        if Md[i] >= m_low and Md[i] <= m_upp:
            ldisk.append(dGalax[i])
    for i in range(len(Mb)):
        if Mb[i] >= m_low and Mb[i] <= m_upp:
            lbar.append(bGalax[i])
    dGalax = ldisk
    bGalax = lbar
    if not vcirc:
        Rd, Qd, e_Qd = read_toomre("data/toomre/", dGalax)
        Rb, Qb, e_Qb = read_toomre("data/toomre/", bGalax)
    else:
        Rd, Qd = read_vcirc_toomre("data/toomre/", dGalax)
        Rb, Qb = read_vcirc_toomre("data/toomre/", bGalax)
    Rd, Qd = filter_dom_regions(dGalax, Rd, Qd, filt='disk')
    Rb, Qb = filter_dom_regions(bGalax, Rb, Qb, filt='bar')
    Rdeff = np.array(get_effectiveR(dGalax))
    Rbeff = np.array(get_effectiveR(bGalax))
    Rd = Rd / Rdeff
    Rb = Rb / Rbeff
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].grid()
    axes[1].grid()
    for i in range(len(Rd)):
        if i == 0:
            axes[0].plot(Rd[i], Qd[i], 'x', color='black', label='Disk')
        else:
            axes[0].plot(Rd[i], Qd[i], 'x', color='black')
    for i in range(len(Rb)):
        if i == 0:
            axes[1].plot(Rb[i], Qb[i], 'x', color='red', label='Bar')
        else:
            axes[1].plot(Rb[i], Qb[i], 'x', color='red')
    axes[0].set_title('No Bar')
    axes[0].set_ylim(0, 3)
    axes[0].set_xlim(0, 8)
    axes[0].set_xlabel(r'R in $[R_e]$')
    axes[0].set_ylabel('Q')
    axes[1].set_title('With Bar')
    axes[1].set_ylim(0, 3)
    axes[1].set_xlim(0, 8)
    axes[1].set_xlabel(r'R in $[R_e]$')
    axes[1].set_ylabel('Q')
    if title != '':
        fig.suptitle(title)
    plt.show()
    plt.close()


def compare_by_redshift(names, z_int1, z_int2):
    """Compare the Toomre parameter by the redshift of the galaxies.
    Args:
        names (list): names of the galaxies
        z_int1 (tuple): first redshift interval (lower_z, upper_z)
        z_int2 (tuple): second redshift interval (lower_z, upper_z)
    Returns:
        None."""
    redz = P.get_redshift(names)
    g1, g2 = [], []
    for i in range(len(names)):
        if redz[i] >= z_int1[0] and redz[i] <= z_int1[1]:
            g1.append(names[i])
        elif redz[i] >= z_int2[0] and redz[i] <= z_int2[1]:
            g2.append(names[i])
    x1, y1, e1 = read_toomre("data/toomre/", g1)
    x2, y2, e2 = read_toomre("data/toomre/", g2)
    x1eff = get_effectiveR(g1)
    x2eff = get_effectiveR(g2)
    x1 = x1 / x1eff
    x2 = x2 / x2eff
    title1 = "Galaxies in redshift interval {:2.3f} to {:2.3f}".format(
        z_int1[0], z_int1[1])
    title2 = "Galaxies in redshift interval {:2.3f} to {:2.3f}".format(
        z_int2[0], z_int2[1])
    plot_against(x1, y1, title1, x2, y2, title2, r'R in $[R_e]$', "Q")


def plot_against(x1, y1, title1, x2, y2, title2, xlabel, ylabel, title=''):
    """Plot two data sets against each other."""
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for i in range(len(x1)):
        if i == 0:
            ax[0].plot(x1[i], y1[i], 'x', color='black')
        else:
            ax[0].plot(x1[i], y1[i], 'x', color='black')
    for i in range(len(x2)):
        if i == 0:
            ax[1].plot(x2[i], y2[i], 'x', color='red')
        else:
            ax[1].plot(x2[i], y2[i], 'x', color='red')
    ax[0].grid()
    ax[1].grid()
    ax[0].set_title(title1)
    ax[0].set_ylim(0, 3)
    ax[0].set_xlim(0, 8)
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax[1].set_title(title2)
    ax[1].set_ylim(0, 3)
    ax[1].set_xlim(0, 8)
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(ylabel)
    if title != '':
        fig.suptitle(title)
    plt.show()
    plt.close()


def filter_dom_regions(names, R, Q, filt='bar'):
    """Get dominated regions of bar or bulge. Return shortened R and Q values.
    Args:
        names (list): names of galaxies
        R (list): list of lists with radius values in [arcs]
        Q (list): list of Toomre parameter profile
        filt (str): either 'bar' or 'disk'
    Returns:
        (list, list): tuple of R and Q."""
    R_ret = []
    Q_ret = []
    for i in range(len(R)):
        if filt == 'bar':
            rdom = P.bar_dominating(names[i])
        elif filt == 'disk':
            rdom = P.disk_dominating(names[i])
        else:
            raise Exception(
                "Wrong argument 'filt={:}', allowed values: 'bar' and 'disk'".format(filt))
        ind = 0
        for r in R[i]:
            if rdom > r:
                ind += 1
            else:
                break
        R_ret.append(R[i][ind:])
        Q_ret.append(Q[i][ind:])
    return R_ret, Q_ret


def filter_radius(R, Q, rmax):
    """Return the values up to a max radius.
    Args:
        R (list): 2D array of the radius values in arcsec
        Q (list): 2D array the y- values
        rmax (float): cut of in arcsec
    Returns:
        np.array, np.array: R, Q - the shortened arrays.
    """
    R = list(R)
    Q = list(Q)
    for i in range(len(R)):
        for j in range(len(R[i])):
            if R[i][j] > rmax:
                R[i] = R[i][:j]
                Q[i] = Q[i][:j]
                break
    return R, Q


if __name__ == "__main__":
    names = read_filenames("data/massmge/", "_mass_mge.txt")
#    mass_overview(names)
    # not in fits file:
    not_there = ['IC1755', 'IC1652', 'IC0540', 'IC0480', 'IC2101', 'IC2247', 'IC2487', 'IC5376', 'MCG-01-54-016', 'MCG-02-02-040', 'MCG-02-03-015', 'NGC0169', 'NGC0177', 'NGC0216', 'NGC0217', 'NGC0429', 'NGC0444', 'NGC0504', 'NGC0529', 'NGC0681', 'NGC0741', 'NGC0781', 'NGC0810', 'NGC0825', 'NGC1056', 'NGC1542', 'NGC2480', 'NGC2481', 'NGC2554', 'NGC3160', 'NGC3303', 'NGC4149', 'NGC4676A', 'NGC4676B', 'NGC4841A', 'NGC4874', 'NGC5216', 'NGC5218', 'NGC5614', 'NGC5797', 'NGC5908', 'NGC5930', 'NGC5934', 'NGC5953', 'NGC5966', 'NGC5987', 'NGC6081', 'NGC6168', 'NGC6310', 'NGC6338', 'NGC7025', 'NGC7436B', 'NGC7608', 'NGC7625', 'NGC7684',
                 'NGC7711', 'NGC7783NED01', 'NGC7800', 'SDSS-C4-DR33247', 'UGC00148', 'UGC00335NED02', 'UGC00809', 'UGC00841', 'UGC01057', 'UGC03151', 'UGC03539', 'UGC03899', 'UGC03969', 'UGC04029', 'UGC04197', 'UGC04280', 'UGC04722', 'UGC05113', 'UGC05498NED01', 'UGC05598', 'UGC05990', 'UGC06036', 'UGC08107', 'UGC08778', 'UGC09537', 'UGC09665', 'UGC09873', 'UGC09892', 'UGC10123', 'UGC10205', 'UGC10257', 'UGC10297', 'UGC10331', 'UGC10380', 'UGC10384', 'UGC10650', 'UGC10710', 'UGC10972', 'UGC11680NED01', 'UGC11717', 'UGC12054', 'UGC12308', 'UGC12518', 'UGC12519', 'UGC12723', 'UGC12857', 'VV488NED02']

    # not in CALIFA_master_selection.txt:
    not_there += ['IC2402', 'LSBCF560-04', 'MCG-02-02-086', 'NGC2484']

    # also not there:
    not_there += ['IC2378', 'NGC0647']

    # for vcirc Q
    not_there += ['NGC0171', 'NGC3381', "NGC5631"]

    for n in not_there:
        names.remove(n)

    # compare_dominated_regions(names)

    # toomre_ellipticals(names)
