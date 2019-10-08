import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import tqdm
import subprocess
import galaxy as g
import photometrics as ph
from scipy import integrate
from astropy.io import fits

P = ph.Photometrics("data/cut_photometrics.fits")


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
        np.array: np.arrays of R, Q"""
    toomre = []
    radius = []  # arcsec
    for n in names:
        f = open(directory + n + '_Q_data.txt', 'r')
        lines = f.readlines()
        R = []
        Q = []
        e_Q = []
        for l in lines:
            l = l.split()
            if len(l) > 4 and l[0] != '#':
                try:
                    R.append(float(l[0]))
                    Q.append(float(l[1]))
                except Exception as e:
                    print(e)
                    print("Error reading Toomre parameter")
        toomre.append(np.array(Q))
        radius.append(np.array(R))
    R = np.array(radius)  # arcsec
    Q = np.array(toomre)
    return R, Q


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


def get_weighted_mean(x, y):
    """Calculate the mean weighted by the density of data points.
    Data points at high density weigh less then those with low density.
    Density is calculated by the cumulative density to all other datapoints.
    Args:
        x (np.array): 2D array of the x values
        y (np.array): 2D array of the y values
    Returns:
        float: the weighted average.
    """
    weights = []
    mean = np.zeros(len(x))
    for i in range(len(x)):
        w = []
        for j in range(len(x[i])):
            w.append(np.sum(np.abs(x[i]-x[i][j])))
        weights.append(np.array(w))
        mean[i] = np.sum(y[i] * weights[i]) / np.sum(weights[i])
    return mean


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
        G = g.Galaxy(names[i])
        if R == None:
            s = np.sqrt(G.data['vzz'])
            e_s = 0.5/s * G.data['e_vzz']
            sigZ.append(s)
            e_sigZ.append(e_s)
        else:
            r = G.data['R']
            s = np.sqrt(G.data['vzz'])
            e_s = 0.5/s * G.data['e_vzz']
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

    dataMS = fits.open("data/sample/CALIFA_2_MS_class.fits")[1].data
    dataES = fits.open("data/sample/CALIFA_2_ES_class.fits")[1].data
    hType = []
    for i in range(len(names)):
        t = ''
        ind = np.where(dataMS['REALNAME'] == names[i])[0]
        if len(ind) != 0:
            ind = ind[0]
            t = dataMS['hubtyp'][ind] + dataMS['hubsubtyp'][ind]
        else:
            ind = np.where(dataES['realname'] == names[i])[0][0]
            t = dataES['hubtyp'][ind] + dataES['hubsubtyp'][ind]
        hType.append(t)
    return hType
    ##############################
    # f = open("data/CALIFA_master_selection.txt", 'r')
    # lines = f.readlines()
    # galaxy = list()
    # hType = list()
    # for l in lines:
    #     l = l.split()
    #     if l[0][0] != '#':
    #         galaxy.append(l[1])
    #         hType.append(l[10])
    # d = {galaxy[i]: hType[i] for i in range(len(galaxy))}

    # return [d[n] for n in names]


def get_redshift(names):
    """Get redshift of galaxies.
    Args:
        names (list): list of galaxie names.
    Returns:
        list: list with redshifts."""
    redshift = np.zeros(len(names))
    for i in range(len(names)):
        z = P.data['REDSHIFT'][P.index[names[i]]]
        redshift[i] = z
    return redshift


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
    for i in range(len(R)):
        for j in range(1, len(mbins)):
            if M[i] <= mbins[j]:
                Rnew[j-1].append(R[i])
                Qnew[j-1].append(Q[i])
                break
    return Rnew, Qnew


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


def filter_type(R, Q, names, hType):
    """Return the R and Q values of galaxies of a specific Hubble Type.
    Args:
        R (np.array): radial values
        Q (np.array): Toomre parameter
        names (list): galaxies to filter
        hType (str or list): hubble type
    Returns:
        np.arrays: R and Q values of galaxies with Hubble type.
    """
    if type(hType) == str:
        hType = [hType]
    rnames = []
    rR = []
    rQ = []
    t = get_type(names)
    for i in range(len(names)):
        if t[i] in hType:
            rnames.append(names[i])
            rR.append(R[i])
            rQ.append(Q[i])
    return np.array(rR), np.array(rQ)


def clearly_classified(names, Q, classif, barred):
    """Return the Q values only for galaxies classified
    as classif.
    Args:
        Q (np.array): the Q values
        classif (str): the classification to return
        barred (dict): dict {'name': 'barredness'}
    Returns:
        np.array.
    """
    rQ = []
    for i in range(len(names)):
        if barred[names[i]] == classif:
            rQ.append(Q[i])
    return np.array(rQ)
