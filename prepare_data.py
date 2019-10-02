import numpy as np
import tool_analyze as ta
import galaxy as g
import photometrics as pm


P = pm.Photometrics("data/photometric_decomposition.fits")


def get_names():
    """Get the filtered samples names, return total, 
    barred, non-barred and ellipticals.
    Args:
        names (list): total names
    Returns:
        total, gElli, gDisk, gBar: list of names.
    """
    names = ta.read_filenames("data/massmge/", "_mass_mge.txt")
    # overlap photometrics and orbital decomposition
    not_there = ['IC1755', 'IC1652', 'IC0540', 'IC0480', 'IC2101', 'IC2247', 'IC2487', 'IC5376', 'MCG-01-54-016', 'MCG-02-02-040', 'MCG-02-03-015', 'NGC0169', 'NGC0177', 'NGC0216', 'NGC0217', 'NGC0429', 'NGC0444', 'NGC0504', 'NGC0529', 'NGC0681', 'NGC0741', 'NGC0781', 'NGC0810', 'NGC0825', 'NGC1056', 'NGC1542', 'NGC2480', 'NGC2481', 'NGC2554', 'NGC3160', 'NGC3303', 'NGC4149', 'NGC4676A', 'NGC4676B', 'NGC4841A', 'NGC4874', 'NGC5216', 'NGC5218', 'NGC5614', 'NGC5797', 'NGC5908', 'NGC5930', 'NGC5934', 'NGC5953', 'NGC5966', 'NGC5987', 'NGC6081', 'NGC6168', 'NGC6310', 'NGC6338', 'NGC7025', 'NGC7436B', 'NGC7608',
                 'NGC7625', 'NGC7684', 'NGC7711', 'NGC7783NED01', 'NGC7800', 'SDSS-C4-DR33247', 'UGC00148', 'UGC00335NED02', 'UGC00809', 'UGC00841', 'UGC01057', 'UGC03151', 'UGC03539', 'UGC03899', 'UGC03969', 'UGC04029', 'UGC04197', 'UGC04280', 'UGC04722', 'UGC05113', 'UGC05498NED01', 'UGC05598', 'UGC05990', 'UGC06036', 'UGC08107', 'UGC08778', 'UGC09537', 'UGC09665', 'UGC09873', 'UGC09892', 'UGC10123', 'UGC10205', 'UGC10257', 'UGC10297', 'UGC10331', 'UGC10380', 'UGC10384', 'UGC10650', 'UGC10710', 'UGC10972', 'UGC11680NED01', 'UGC11717', 'UGC12054', 'UGC12308', 'UGC12518', 'UGC12519', 'UGC12723', 'UGC12857', 'VV488NED02']

    print("Analysis with the modelled vcirc\n" + '*'*32)
    print("Galaxy infos:")

    for n in not_there:
        names.remove(n)
    print('Overlap orbital and photometric decomposition:', len(names))

    double_galax = P.get_doubles()

    for n in double_galax:
        try:
            names.remove(n)
        except:
            # print(n)
            pass
    print("Excluding ambiguous photometric components:", len(names))
    # not in master file
    not_there = ['IC2402', 'LSBCF560-04',
                 'MCG-02-02-086', 'NGC2484', 'IC2378', 'NGC0647']
    for n in not_there:
        try:
            names.remove(n)
        except:
            pass
    print("No Hubble Type classification", len(names))
    # because inclination angle to low for deprojection for vcirc modelling
    not_there = ['NGC0171', 'NGC3381', 'NGC5631']
    for n in not_there:
        names.remove(n)
    print("Vcirc modelling cut:", len(names))
    gBar = ta.get_barred_galaxies(names)
    print("Barred Galaxies:", len(gBar))
    gDisk = ta.get_disk_galaxies(names)
    print("Non Barred Galaxies:", len(gDisk))
    gElli = ta.get_elliptical_galaxies(names)
    print("Ellipticals:", len(gElli))
    return names, gElli, gDisk, gBar


def read_interpolated_data(names, extension, directory, c=1):
    """Read the interpolated data from file. Return the first two columns.
    Args:
        names (list): names of the galaxy
        extension (str): name of the file extension
        directory (str): location of the file
        c (int): second column, 1 == Q or X, -1 == sigZ
    Returns:
        2D np.arrays.
    """
    R = np.zeros((len(names), 300))
    Y = np.zeros((len(names), 300))
    for i in range(len(names)):
        f = open(directory + names[i] + extension, 'r')
        lines = f.readlines()
        l = lines[2:]
        for j in range(len(l)):
            l[j] = l[j].split()
            R[i][j] = l[j][0]
            Y[i][j] = l[j][c]
    return R, Y


names, gElli, gDisk, gBar = get_names()

#Rbar, Qbar = ta.read_vcirc_toomre("data/toomre/", gBar)
Rbar, Qbar = ta.read_toomre("data/toomre/", gBar)
Rbar, Xbar = ta.read_swing_ampl(gBar)
Mbar = ta.read_stellar_mass(gBar)
#Rdisk, Qdisk = ta.read_vcirc_toomre("data/toomre/", gDisk)
Rdisk, Qdisk = ta.read_toomre("data/toomre/", gDisk)
Rdisk, Xdisk = ta.read_swing_ampl(gDisk)
Mdisk = ta.read_stellar_mass(gDisk)
#Relli, Qelli = ta.read_vcirc_toomre("data/toomre/", gElli)
Relli, Qelli = ta.read_toomre("data/toomre/", gElli)
Relli, Xelli = ta.read_swing_ampl(gElli)
Melli = ta.read_stellar_mass(gElli)


# cut at arcsec radius of 35
rmax = 35  # arcsec
tempR, Qbar = ta.filter_radius(Rbar, Qbar, rmax)
Rbar, Xbar = ta.filter_radius(Rbar, Xbar, rmax)
tempR, Qdisk = ta.filter_radius(Rdisk, Qdisk, rmax)
Rdisk, Xdisk = ta.filter_radius(Rdisk, Xdisk, rmax)
tempR, Qelli = ta.filter_radius(Relli, Qelli, rmax)
Relli, Xelli = ta.filter_radius(Relli, Xelli, rmax)

# no filter for ellipticals needed
tempR, Qbar = ta.filter_dom_regions(gBar, Rbar, Qbar, filt='bar')
Rbar, Xbar = ta.filter_dom_regions(gBar, Rbar, Xbar, filt='bar')
tempR, Qdisk = ta.filter_dom_regions(gDisk, Rdisk, Qdisk, filt='disk')
Rdisk, Xdisk = ta.filter_dom_regions(gDisk, Rdisk, Xdisk, filt='disk')


# before converting to effective radius
sigZbar, e_sigZbar = ta.get_sigma_Z(gBar, Rbar)
sigZdisk, e_sigZdisk = ta.get_sigma_Z(gDisk, Rdisk)
sigZelli, e_sigZelli = ta.get_sigma_Z(gElli, Relli)

# convert to effective radius
REbar = ta.get_effectiveR(gBar)
REdisk = ta.get_effectiveR(gDisk)
REelli = ta.get_effectiveR(gElli)

Rbar = np.array(Rbar) / np.array(REbar)
Rdisk = np.array(Rdisk) / np.array(REdisk)

for i in range(len(Relli)):
    Relli[i] = Relli[i] / REelli[i]
Relli = np.array(Relli)


# integrated quantities

# interpolated data for mean
iRbar, iQbar = read_interpolated_data(
    gBar, '_Qintpol_data.txt', "data/toomre/")
iRbar, iXbar = read_interpolated_data(gBar, '_Xintpol_data.txt', "data/swing/")
iRbar, iSigZbar = read_interpolated_data(
    gBar, '_Qintpol_data.txt', "data/toomre/", -1)

iRdisk, iQdisk = read_interpolated_data(
    gDisk, '_Qintpol_data.txt', "data/toomre/")
iRdisk, iXdisk = read_interpolated_data(
    gDisk, '_Xintpol_data.txt', "data/swing/")
iRdisk, iSigZdisk = read_interpolated_data(
    gDisk, '_Qintpol_data.txt', "data/toomre/", -1)

iRelli, iQelli = read_interpolated_data(
    gElli, '_Qintpol_data.txt', "data/toomre/")
iRelli, iXelli = read_interpolated_data(
    gElli, '_Xintpol_data.txt', "data/swing/")
iRelli, iSigZelli = read_interpolated_data(
    gElli, '_Qintpol_data.txt', "data/toomre/", -1)


# integrated
tempR, iQbar = ta.filter_radius(iRbar, iQbar, rmax)
tempR, iSigZbar = ta.filter_radius(iRbar, iSigZbar, rmax)
iRbar, iXbar = ta.filter_radius(iRbar, iXbar, rmax)

tempR, iQdisk = ta.filter_radius(iRdisk, iQdisk, rmax)
tempR, iSigZdisk = ta.filter_radius(iRdisk, iSigZdisk, rmax)
iRdisk, iXdisk = ta.filter_radius(iRdisk, iXdisk, rmax)

tempR, iQelli = ta.filter_radius(iRelli, iQelli, rmax)
tempR, iSigZelli = ta.filter_radius(iRelli, iSigZelli, rmax)
iRelli, iXelli = ta.filter_radius(iRelli, iXelli, rmax)

tempR, iQbar = ta.filter_dom_regions(gBar, iRbar, iQbar, filt='bar')
tempR, iSigZbar = ta.filter_dom_regions(gBar, iRbar, iSigZbar, filt='bar')
iRbar, iXbar = ta.filter_dom_regions(gBar, iRbar, iXbar, filt='bar')

tempR, iQdisk = ta.filter_dom_regions(gDisk, iRdisk, iQdisk, filt='disk')
tempR, iSigZdisk = ta.filter_dom_regions(gDisk, iRdisk, iSigZdisk, filt='disk')
iRdisk, iXdisk = ta.filter_dom_regions(gDisk, iRdisk, iXdisk, filt='disk')


iRbar = iRbar / np.array(REbar)
iRdisk = iRdisk / np.array(REdisk)
iRelli = iRelli / np.array(REelli)

Qbar_median, Qbar_mean, Qbar_std = ta.get_median_mean_std(iQbar)
Xbar_median, Xbar_mean, Xbar_std = ta.get_median_mean_std(iXbar)
sigZbar_median, sigZbar_mean, sigZbar_std = ta.get_median_mean_std(iSigZbar)

Qdisk_median, Qdisk_mean, Qdisk_std = ta.get_median_mean_std(iQdisk)
Xdisk_median, Xdisk_mean, Xdisk_std = ta.get_median_mean_std(iXdisk)
sigZdisk_median, sigZdisk_mean, sigZdisk_std = ta.get_median_mean_std(
    iSigZdisk)

Qelli_median, Qelli_mean, Qelli_std = ta.get_median_mean_std(iQelli)
Xelli_median, Xelli_mean, Xelli_std = ta.get_median_mean_std(iXelli)
sigZelli_median, sigZelli_mean, sigZelli_std = ta.get_median_mean_std(
    iSigZelli)
