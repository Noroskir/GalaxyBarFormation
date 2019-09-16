import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tool_analyze as ta
import galaxy as g
import photometrics as pm


names = ta.read_filenames("data/massmge/", "_mass_mge.txt")
# overlap photometrics and orbital decomposition
not_there = ['IC1755', 'IC1652', 'IC0540', 'IC0480', 'IC2101', 'IC2247', 'IC2487', 'IC5376', 'MCG-01-54-016', 'MCG-02-02-040', 'MCG-02-03-015', 'NGC0169', 'NGC0177', 'NGC0216', 'NGC0217', 'NGC0429', 'NGC0444', 'NGC0504', 'NGC0529', 'NGC0681', 'NGC0741', 'NGC0781', 'NGC0810', 'NGC0825', 'NGC1056', 'NGC1542', 'NGC2480', 'NGC2481', 'NGC2554', 'NGC3160', 'NGC3303', 'NGC4149', 'NGC4676A', 'NGC4676B', 'NGC4841A', 'NGC4874', 'NGC5216', 'NGC5218', 'NGC5614', 'NGC5797', 'NGC5908', 'NGC5930', 'NGC5934', 'NGC5953', 'NGC5966', 'NGC5987', 'NGC6081', 'NGC6168', 'NGC6310', 'NGC6338', 'NGC7025', 'NGC7436B', 'NGC7608',
             'NGC7625', 'NGC7684', 'NGC7711', 'NGC7783NED01', 'NGC7800', 'SDSS-C4-DR33247', 'UGC00148', 'UGC00335NED02', 'UGC00809', 'UGC00841', 'UGC01057', 'UGC03151', 'UGC03539', 'UGC03899', 'UGC03969', 'UGC04029', 'UGC04197', 'UGC04280', 'UGC04722', 'UGC05113', 'UGC05498NED01', 'UGC05598', 'UGC05990', 'UGC06036', 'UGC08107', 'UGC08778', 'UGC09537', 'UGC09665', 'UGC09873', 'UGC09892', 'UGC10123', 'UGC10205', 'UGC10257', 'UGC10297', 'UGC10331', 'UGC10380', 'UGC10384', 'UGC10650', 'UGC10710', 'UGC10972', 'UGC11680NED01', 'UGC11717', 'UGC12054', 'UGC12308', 'UGC12518', 'UGC12519', 'UGC12723', 'UGC12857', 'VV488NED02']


print("Analysis with the modelled vcirc\n" + '*'*32)
print("Galaxy infos:")

for n in not_there:
    names.remove(n)
print('Overlap orbital and photometric decomposition:', len(names))

P = pm.Photometrics("data/photometric_decomposition.fits")
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


Rbar, Qbar = ta.read_vcirc_toomre("data/toomre/", gBar)
Rbar, Xbar = ta.read_swing_ampl(gBar)
Mbar = ta.read_stellar_mass(gBar)
Rdisk, Qdisk = ta.read_vcirc_toomre("data/toomre/", gDisk)
Rdisk, Xdisk = ta.read_swing_ampl(gDisk)
Mdisk = ta.read_stellar_mass(gDisk)
Relli, Qelli = ta.read_vcirc_toomre("data/toomre/", gElli)
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


Qbar_median, Qbar_mean, Qbar_std = ta.get_median_mean_std(Qbar)
Xbar_median, Xbar_mean, Xbar_std = ta.get_median_mean_std(Xbar)
sigZbar_median, sigZbar_mean, sigZbar_std = ta.get_median_mean_std(sigZbar)

Qdisk_median, Qdisk_mean, Qdisk_std = ta.get_median_mean_std(Qdisk)
Xdisk_median, Xdisk_mean, Xdisk_std = ta.get_median_mean_std(Xdisk)
sigZdisk_median, sigZdisk_mean, sigZdisk_std = ta.get_median_mean_std(sigZdisk)

Qelli_median, Qelli_mean, Qelli_std = ta.get_median_mean_std(Qelli)
Xelli_median, Xelli_mean, Xelli_std = ta.get_median_mean_std(Xelli)
sigZelli_median, sigZelli_mean, sigZelli_std = ta.get_median_mean_std(sigZelli)


def plot_error_hist():
    bins = np.linspace(0, 1, 10)
    n, bins, patches = plt.hist([Xbar_std, Xdisk_std, Xelli_std], bins, label=[
        'Bar', 'No Bar', 'Ellipticals'])
    plt.title("Distribution of the standard deviation of X")
    plt.legend()
    plt.show()
    plt.close()


def plot_Q_profile(name):
    G = g.Galaxy("data/sRsz_data/", "data/massmge/", name)
    Q, eQ = G.get_toomre_vcirc()
    R = G.data_cyl['R']
    plt.plot(R, Q, color='red')
    plt.errorbar(R, Q, yerr=eQ, fmt='.-', color='red', capsize=2)
    plt.grid()
    plt.title("Toomre Parameter Q " + name)
    plt.xlabel(r'Radius [arcsec]')
    plt.ylabel(r'Q')
    plt.show()


def plot_X_profile(name):
    G = g.Galaxy("data/sRsz_data/", "data/massmge/", name)
    R = G.data_cyl['R']
    X = G.swing_amplif_X()
    plt.plot(R, X)
    plt.title(r'Swing Amplification Parameter X '+name)
    plt.ylabel(r'X')
    plt.xlabel(r'R [arcsec]')
    plt.grid()
    plt.show()
    plt.close()


def plot_Q_X_profile(name):
    G = g.Galaxy("data/sRsz_data/", "data/massmge/", name)
    R = G.data_cyl['R']
    X = G.swing_amplif_X()
    Q, eQ = G.get_toomre_vcirc()
    plt.plot(R, Q, color='red', label='Q')
    plt.errorbar(R, Q, yerr=eQ, fmt='.-', color='red', capsize=2)
    plt.plot(R, X, color='black', label='X')
    plt.title(r'Q and X '+name)
    plt.xlabel(r'R [arcsec]')
    plt.ylabel(r'Q and X')
    plt.grid()
    plt.legend()
    plt.show()
    plt.close()


def plot_sigR_profile(name):
    G = g.Galaxy("data/sRsz_data/", "data/massmge/", name)
    R = G.data_cyl['R']
    sigR = np.sqrt(G.data_cyl['vRR'])
    e_sigR = 0.5 * 1/sigR * G.data_cyl['e_vRR']
    plt.plot(R, sigR, color='red')
    plt.errorbar(R, sigR, yerr=e_sigR, fmt='.-',
                 color='red', capsize=2)
    plt.grid()
    plt.title(r"Radial Velocity Dispersion $\sigma_R$ " + name)
    plt.xlabel(r'Radius [arcsec]')
    plt.ylabel(r'$\sigma_R$ [km/s]')
    plt.show()


def plot_vcirc_profile(name):
    G = g.Galaxy("data/sRsz_data/", "data/massmge/", name)
    R = G.data_cyl['R']
    vcirc = G.velocity_vcirc()
    plt.plot(R, vcirc, color='red')
    plt.grid()
    plt.title(r'Circular Velocity ' + name)
    plt.xlabel(r'Radius [arcsec]')
    plt.ylabel(r'$v_c$ [km/s]')
    plt.show()
    plt.close()


def plot_SMD_profile(name):
    G = g.Galaxy("data/sRsz_data/", "data/massmge/", name)
    R = G.data_cyl['R']
    smd = G.get_surfaceMD('all')
    smd_lm = G.get_surfaceMD('lm')
    smd_dm = G.get_surfaceMD('dm')
    plt.semilogy(R, smd_lm, '--', color='red', label='Stellar Mass')
    plt.semilogy(R, smd_dm, '--', color='black', label='Dark Matter')
    plt.semilogy(R, smd, color='blue', label='Total')
    plt.title(r'Surface Mass Density')
    plt.xlabel('R [arcsec]')
    plt.ylabel(r'$\Sigma$ [M$_{sun}$/pc$^2$]')
    plt.grid()
    plt.legend()
    plt.show()
    plt.close()


def plot_Q_total_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].plot(Relli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].plot(Rdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].plot(Rbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, 3.5)
        ax[i].set_ylim(0, 3.5)
        ax[i].set_xlabel(r'R/R$_e$')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/plot_Q_total_mass.pdf")
    plt.show()
    plt.close()


def plot_X_total_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].plot(Relli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].plot(Rdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].plot(Rbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, 2.5)
        ax[i].set_ylim(0, 2.5)
        ax[i].set_xlabel(r'R/R$_e$')
        ax[i].set_ylabel(r'X')
    plt.suptitle('Swing Amplification Parameter X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/plot_X_total_mass.pdf")
    plt.show()
    plt.close()


def plot_Q_total_type():
    tBar = ta.get_type(gBar)
    tDisk = ta.get_type(gDisk)
    types = ['E', '0', 'a', 'b', 'c', 'd']
    cmap = plt.cm.get_cmap('jet_r')
    valmax = len(types)
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    for i in range(len(Relli)):
        c = cmap(0)
        ax[0].plot(Relli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].plot(Rdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].plot(Rbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, 3.5)
        ax[i].set_ylim(0, 3.5)
        ax[i].set_xlabel(r'R/R$_e$')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=0, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ticks=np.linspace(
        0.5, valmax+0.5, valmax+1), ax=ax.ravel().tolist(), boundaries=np.arange(0, valmax+1, 1))
    clb.ax.set_yticklabels(['E']+['S' + t for t in types[1:]])
    clb.set_label('Hubble Type')
    plt.savefig("figures/report/plot_Q_total_type.pdf")
    plt.show()
    plt.close()


def plot_X_total_type():
    tBar = ta.get_type(gBar)
    tDisk = ta.get_type(gDisk)
    types = ['E', '0', 'a', 'b', 'c', 'd']
    cmap = plt.cm.get_cmap('jet_r')
    valmax = len(types)
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    for i in range(len(Relli)):
        c = cmap(0)
        ax[0].plot(Relli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].plot(Rdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].plot(Rbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, 2.5)
        ax[i].set_ylim(0, 2.5)
        ax[i].set_xlabel(r'R/R$_e$')
        ax[i].set_ylabel(r'X')
    plt.suptitle('Swing Amplification Parameter X')
    norm = mpl.colors.Normalize(vmin=0, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ticks=np.linspace(
        0.5, valmax+0.5, valmax+1), ax=ax.ravel().tolist(), boundaries=np.arange(0, valmax+1, 1))
    clb.ax.set_yticklabels(['E']+['S' + t for t in types[1:]])
    clb.set_label('Hubble Type')
    plt.savefig("figures/report/plot_X_total_type.pdf")
    plt.show()
    plt.close()


def plot_Q_sigZ_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].plot(sigZelli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].plot(sigZdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].plot(sigZbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_xlim(0, 3.5)
        ax[i].set_ylim(0, 3.5)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/plot_Q_sigZ_mass.pdf")
    plt.show()
    plt.close()


def plot_Q_sigZ_type():
    tBar = ta.get_type(gBar)
    tDisk = ta.get_type(gDisk)
    types = ['E', '0', 'a', 'b', 'c', 'd']
    cmap = plt.cm.get_cmap('jet_r')
    valmax = len(types)
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    for i in range(len(Relli)):
        c = cmap(0)
        ax[0].plot(sigZelli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].plot(sigZdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].plot(sigZbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_xlim(0, 3.5)
        ax[i].set_ylim(0, 3.5)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=0, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ticks=np.linspace(
        0.5, valmax+0.5, valmax+1), ax=ax.ravel().tolist(), boundaries=np.arange(0, valmax+1, 1))
    clb.ax.set_yticklabels(['E']+['S' + t for t in types[1:]])
    clb.set_label('Hubble Type')
    plt.savefig("figures/report/plot_Q_sigZ_type.pdf")
    plt.show()
    plt.close()


def plot_X_sigZ_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].plot(sigZelli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].plot(sigZdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].plot(sigZbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_xlim(0, 2.5)
        ax[i].set_ylim(0, 2.5)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'X')
    plt.suptitle('Swing Amplification Parameter X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/plot_X_sigZ_mass.pdf")
    plt.show()
    plt.close()


def plot_X_sigZ_type():
    tBar = ta.get_type(gBar)
    tDisk = ta.get_type(gDisk)
    types = ['E', '0', 'a', 'b', 'c', 'd']
    cmap = plt.cm.get_cmap('jet_r')
    valmax = len(types)
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    for i in range(len(Relli)):
        c = cmap(0)
        ax[0].plot(sigZelli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].plot(sigZdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].plot(sigZbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_xlim(0, 2.5)
        ax[i].set_ylim(0, 2.5)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'X')
    plt.suptitle('Swing Amplification Parameter X')
    norm = mpl.colors.Normalize(vmin=0, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ticks=np.linspace(
        0.5, valmax+0.5, valmax+1), ax=ax.ravel().tolist(), boundaries=np.arange(0, valmax+1, 1))
    clb.ax.set_yticklabels(['E']+['S' + t for t in types[1:]])
    clb.set_label('Hubble Type')
    plt.savefig("figures/report/plot_X_sigZ_type.pdf")
    plt.show()
    plt.close()


def stack_Q_profile_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    nbins = 5
    step = (valmax-valmin) / nbins
    mbins = [valmin + i * step for i in range(nbins+1)]
    RSbar, QSbar, = ta.mass_bin(Rbar, Qbar, Mbar, mbins)
    RSdisk, QSdisk, = ta.mass_bin(Rdisk, Qdisk, Mdisk, mbins)
    RSelli, QSelli, = ta.mass_bin(Relli, Qelli, Melli, mbins)
    step = 0.5
    for i in range(len(RSelli)):
        c = cmap((i*step+step/2)/(valmax-valmin))
        R, Q = ta.stack_curves(RSelli[i], QSelli[i], step, threshold=1)
        ax[0].plot(R, Q, color=c)
    for i in range(len(RSdisk)):
        c = cmap((i*step+step/2)/(valmax-valmin))
        R, Q = ta.stack_curves(RSdisk[i], QSdisk[i], step, threshold=1)
        ax[1].plot(R, Q, color=c)
    for i in range(len(RSbar)):
        c = cmap((i*step+step/2)/(valmax-valmin))
        R, Q = ta.stack_curves(RSbar[i], QSbar[i], step, threshold=1)
        ax[2].plot(R, Q, color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, 3.5)
        ax[i].set_ylim(0, 3.5)
        ax[i].set_xlabel(r'R/R$_e$')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/stack_Q_profile_mass.pdf")
    plt.show()
    plt.close()


def hist_Q_mean():
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    # ax[0].hist(Qelli_mean)
    # ax[0].set(title='Elliptical Galaxies', xlabel='Q')
    bins = np.linspace(0, 2.5, 11)
    ax.hist(Qdisk_mean, bins=bins, label='No Bar', alpha=0.7)
    ax.margins(0.05)
    ax.hist(Qbar_mean, bins=bins, label='Bar', alpha=0.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='Q')
    plt.savefig("figures/report/hist_Q_mean.pdf")
    plt.show()
    plt.close()


def hist_Q_median():
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    # ax[0].hist(Qelli_median)
    # ax[0].set(title='Elliptical Galaxies', xlabel='Q')
    bins = np.linspace(0, 2.5, 11)
    ax.hist(Qdisk_median, bins=bins, label='No Bar', alpha=0.7)
    ax.margins(0.05)
    ax.hist(Qbar_median, bins=bins, label='Bar', alpha=0.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='Q')
    plt.savefig("figures/report/hist_Q_median.pdf")
    plt.show()
    plt.close()


def hist_X_mean():
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    # ax[0].hist(Xelli_mean)
    # ax[0].set(title='Elliptical Galaxies', xlabel='X')
    ax.hist(Xdisk_median, label='No Bar', alpha=.7)
    ax.hist(Xbar_median, label='Bar', alpha=.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='X')
    plt.savefig("figures/report/hist_X_mean.pdf")
    plt.show()
    plt.close()


def hist_X_median():
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    ax.hist(Xdisk_median, label='No Bar', alpha=.7)
    ax.hist(Xbar_median, label='Bar', alpha=.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='X')
    plt.savefig("figures/report/hist_X_median.pdf")
    plt.show()
    plt.close()


def scatter_Q_mass_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(sigZelli_mean, np.append(
        sigZbar_mean, sigZdisk_mean)))
    valmin = min(np.append(sigZelli_mean, np.append(
        sigZbar_mean, sigZdisk_mean)))
    for i in range(len(Melli)):
        c = cmap((sigZelli_mean[i]-valmin)/(valmax-valmin))
        ax[0].scatter(Melli[i], Qelli_mean[i], color=c)
    for i in range(len(Mdisk)):
        c = cmap((sigZdisk_mean[i]-valmin)/(valmax-valmin))
        ax[1].scatter(Mdisk[i], Qdisk_mean[i], color=c)
    for i in range(len(Mbar)):
        c = cmap((sigZbar_mean[i]-valmin)/(valmax-valmin))
        ax[2].scatter(Mbar[i], Qbar_mean[i], color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_xlim(0, 3.5)
        # ax[i].set_ylim(0, 3.5)
        ax[i].set_xlabel(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q vs Stellar Mass')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'$\sigma_z$ [km/s]')
    plt.savefig("figures/report/scatter_Q_mass_mean.pdf")
    plt.show()
    plt.close()


def scatter_X_mass_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(sigZelli_mean, np.append(
        sigZbar_mean, sigZdisk_mean)))
    valmin = min(np.append(sigZelli_mean, np.append(
        sigZbar_mean, sigZdisk_mean)))
    for i in range(len(Melli)):
        c = cmap((sigZelli_mean[i]-valmin)/(valmax-valmin))
        ax[0].scatter(Melli[i], Xelli_mean[i], color=c)
    for i in range(len(Mdisk)):
        c = cmap((sigZdisk_mean[i]-valmin)/(valmax-valmin))
        ax[1].scatter(Mdisk[i], Xdisk_mean[i], color=c)
    for i in range(len(Mbar)):
        c = cmap((sigZbar_mean[i]-valmin)/(valmax-valmin))
        ax[2].scatter(Mbar[i], Xbar_mean[i], color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_xlim(0, 3.5)
        # ax[i].set_ylim(0, 3.5)
        ax[i].set_xlabel(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
        ax[i].set_ylabel(r'X')
    plt.suptitle('Swing Amplification Parameter X vs Stellar Mass')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'$\sigma_z$ [km/s]')
    plt.savefig("figures/report/scatter_X_mass_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_sigZ_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(gElli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].scatter(sigZelli_mean[i], Qelli_mean[i], color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].scatter(sigZdisk_mean[i], Qdisk_mean[i], color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].scatter(sigZbar_mean[i], Qbar_mean[i], color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/scatter_Q_sigZ_mean.pdf")
    plt.show()
    plt.close()


def scatter_X_sigZ_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(gElli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].scatter(sigZelli_mean[i], Xelli_mean[i], color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].scatter(sigZdisk_mean[i], Xdisk_mean[i], color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].scatter(sigZbar_mean[i], Xbar_mean[i], color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/scatter_X_sigZ_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_lambRe_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    lambReElli = ta.read_lambda_re(gElli)
    lambReDisk = ta.read_lambda_re(gDisk)
    lambReBar = ta.read_lambda_re(gBar)
    for i in range(len(gElli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].scatter(lambReElli[i], Qelli_mean[i], color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].scatter(lambReDisk[i], Qdisk_mean[i], color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].scatter(lambReBar[i], Qbar_mean[i], color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].set_xlabel(r'$\lambda_{Re}$')
        ax[i].set_ylabel(r'Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/scatter_Q_lambRe_mean.pdf")
    plt.show()
    plt.close()


def scatter_X_lambRe_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    lambReElli = ta.read_lambda_re(gElli)
    lambReDisk = ta.read_lambda_re(gDisk)
    lambReBar = ta.read_lambda_re(gBar)
    for i in range(len(gElli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].scatter(lambReElli[i], Xelli_mean[i], color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].scatter(lambReDisk[i], Xdisk_mean[i], color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].scatter(lambReBar[i], Xbar_mean[i], color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Witout Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].set_xlabel(r'$\lambda_{Re}$')
        ax[i].set_ylabel(r'X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/scatter_X_lambRe_mean.pdf")
    plt.show()
    plt.close()


def check_mass():
    MAbar = []
    for i in range(len(gBar)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gBar[i])
        MAbar.append(np.log10(G.get_stellar_mass()))
    plt.scatter(Mbar, MAbar)
    plt.xlabel('Stellar Mass')
    plt.ylabel('Dynamical Mass')
    plt.show()
    plt.close()


if __name__ == '__main__':
    # n = 'NGC6278'
    # plot_Q_profile(n)
    # plot_X_profile(n)
    # plot_Q_X_profile(n)
    # plot_sigR_profile(n)
    # plot_vcirc_profile(n)
    # plot_SMD_profile(n)
    # plot_Q_total_mass()
    # plot_X_total_mass()
    # plot_Q_total_type()
    # plot_X_total_type()
    # plot_Q_sigZ_mass()
    # plot_Q_sigZ_type()
    # plot_X_sigZ_mass()
    # plot_X_sigZ_type()
    # stack_Q_profile_mass()
    # hist_Q_mean()
    # hist_Q_median()
    # hist_X_mean()
    # hist_X_median()
    scatter_Q_sigZ_mean()
    scatter_X_sigZ_mean()
    scatter_Q_mass_mean()
    scatter_X_mass_mean()
    scatter_Q_lambRe_mean()
    scatter_X_lambRe_mean()
