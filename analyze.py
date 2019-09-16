import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import photometrics as ph
import galaxy as g
import gaussFilter as gf
import tool_analyze as ta

names = ta.read_filenames("data/massmge/", "_mass_mge.txt")
not_there = ['IC1755', 'IC1652', 'IC0540', 'IC0480', 'IC2101', 'IC2247', 'IC2487', 'IC5376', 'MCG-01-54-016', 'MCG-02-02-040', 'MCG-02-03-015', 'NGC0169', 'NGC0177', 'NGC0216', 'NGC0217', 'NGC0429', 'NGC0444', 'NGC0504', 'NGC0529', 'NGC0681', 'NGC0741', 'NGC0781', 'NGC0810', 'NGC0825', 'NGC1056', 'NGC1542', 'NGC2480', 'NGC2481', 'NGC2554', 'NGC3160', 'NGC3303', 'NGC4149', 'NGC4676A', 'NGC4676B', 'NGC4841A', 'NGC4874', 'NGC5216', 'NGC5218', 'NGC5614', 'NGC5797', 'NGC5908', 'NGC5930', 'NGC5934', 'NGC5953', 'NGC5966', 'NGC5987', 'NGC6081', 'NGC6168', 'NGC6310', 'NGC6338', 'NGC7025', 'NGC7436B', 'NGC7608',
             'NGC7625', 'NGC7684', 'NGC7711', 'NGC7783NED01', 'NGC7800', 'SDSS-C4-DR33247', 'UGC00148', 'UGC00335NED02', 'UGC00809', 'UGC00841', 'UGC01057', 'UGC03151', 'UGC03539', 'UGC03899', 'UGC03969', 'UGC04029', 'UGC04197', 'UGC04280', 'UGC04722', 'UGC05113', 'UGC05498NED01', 'UGC05598', 'UGC05990', 'UGC06036', 'UGC08107', 'UGC08778', 'UGC09537', 'UGC09665', 'UGC09873', 'UGC09892', 'UGC10123', 'UGC10205', 'UGC10257', 'UGC10297', 'UGC10331', 'UGC10380', 'UGC10384', 'UGC10650', 'UGC10710', 'UGC10972', 'UGC11680NED01', 'UGC11717', 'UGC12054', 'UGC12308', 'UGC12518', 'UGC12519', 'UGC12723', 'UGC12857', 'VV488NED02']

not_there += ['IC2402', 'LSBCF560-04',
              'MCG-02-02-086', 'NGC2484', 'IC2378', 'NGC0647']

# because inclination angle to low for deprojection for vcirc modelling
not_there += ['NGC0171', 'NGC3381']


for n in not_there:
    names.remove(n)


gBar = ta.get_barred_galaxies(names)
gDisk = ta.get_disk_galaxies(names)
gElli = ta.get_elliptical_galaxies(names)

# print(len(gBar+gDisk))

# Rbar, Qbar, ebar = ta.read_toomre("data/toomre/", gBar)
# Rdisk, Qdisk, edisk = ta.read_toomre("data/toomre/", gDisk)

Rbar, Qbar = ta.read_vcirc_toomre("data/toomre/", gBar)
Rdisk, Qdisk = ta.read_vcirc_toomre("data/toomre/", gDisk)
Relli, Qelli = ta.read_vcirc_toomre("data/toomre/", gElli)

Rbar, Qbar = ta.filter_dom_regions(gBar, Rbar, Qbar, filt='bar')
Rdisk, Qdisk = ta.filter_dom_regions(gDisk, Rdisk, Qdisk, filt='disk')

REbar = ta.get_effectiveR(gBar)
REdisk = ta.get_effectiveR(gDisk)
REelli = ta.get_effectiveR(gElli)

Rbar = np.array(Rbar) / np.array(REbar)
Rdisk = np.array(Rdisk) / np.array(REdisk)
Relli = np.array(Relli) / np.array(REelli)

Mbar = np.array(ta.read_mass(gBar))   # in solar masses
Mdisk = np.array(ta.read_mass(gDisk))


sig = 0.35

QbarF = []
QdiskF = []
for i in range(len(Qbar)):
    QbarF.append(gf.smooth_data(Rbar[i], Qbar[i], sig))
for i in range(len(Qdisk)):
    QdiskF.append(gf.smooth_data(Rdisk[i], Qdisk[i], sig))


Qbar_median = []
Qbar_mean = []
Qbar_std = []
# QbarVc = []
# QbarVc_median = []
# QbarVc_mean = []
# QbarVc_std = []
# R, QbarVc = ta.read_vcirc_toomre("data/toomre/", gBar)
# R, QbarVc = ta.filter_dom_regions(gBar, R, QbarVc, filt='bar')

for i in range(len(Qbar)):
    Qbar_median.append(np.median(Qbar[i]))
    Qbar_mean.append(np.mean(Qbar[i]))
    Qbar_std.append(np.std(Qbar[i], ddof=1))
    # QbarVc_median.append(np.median(QbarVc[i]))
    # QbarVc_mean.append(np.mean(QbarVc[i]))
    # QbarVc_std.append(np.std(QbarVc[i], ddof=1))

Qdisk_median = []
Qdisk_mean = []
Qdisk_std = []
# R, QdiskVc = ta.read_vcirc_toomre("data/toomre/", gDisk)
# R, QdiskVc = ta.filter_dom_regions(gDisk, R, QdiskVc, filt='disk')
# QdiskVc_median = []
# QdiskVc_mean = []
# QdiskVc_std = []
for i in range(len(Qdisk)):
    Qdisk_median.append(np.median(Qdisk[i]))
    Qdisk_mean.append(np.mean(Qdisk[i]))
    Qdisk_std.append(np.std(Qdisk[i], ddof=1))
    # QdiskVc_median.append(np.median(QdiskVc[i]))
    # QdiskVc_mean.append(np.mean(QdiskVc[i]))
    # QdiskVc_std.append(np.std(QdiskVc[i], ddof=1))


def plot_Q(vcirc=False, smooth=False):
    cmap = plt.cm.get_cmap('jet')
    maxval = 2
    for i in range(len(Rbar)):
        col = cmap(Qbar_median[i]/maxval)
        if smooth:
            plt.semilogx(Rbar[i], QbarF[i],  color=col)
        else:
            plt.semilogx(Rbar[i], Qbar[i],  color=col)
    plt.grid()
    plt.xlim(5e-1, 10)
    plt.ylim(0, 6)
    norm = mpl.colors.Normalize(vmin=0, vmax=maxval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ticks=np.linspace(0, maxval, 11),
                 boundaries=np.arange(-0.05, 2.1, .1))
    plt.show()


def plot_Q_mu0(Qdisk, Qbar):
    cmap = plt.cm.get_cmap('jet')
    mu = []
    for i in range(len(Rbar)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gBar[i])
        mu0 = G.mge_LM(0, 0)
        mu.append(mu0)
    for i in range(len(Rdisk)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gDisk[i])
        mu0 = G.mge_LM(0, 0)
        mu.append(mu0)
    mu = np.array(mu)
    maxval = np.max(mu)
    minval = np.min(mu)
    maxval = 12000
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for i in range(len(Rbar)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gBar[i])
        mu0 = G.mge_LM(0, 0)
        col = cmap(mu0/maxval)
        ax[1].semilogx(Rbar[i], Qbar[i], '-', ms=4, color=col)
    for i in range(len(Rdisk)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gDisk[i])
        mu0 = G.mge_LM(0, 0)
        col = cmap(mu0/maxval)
        ax[0].semilogx(Rdisk[i], Qdisk[i], '-', ms=4, color=col)
    ax[0].grid()
    ax[0].set_ylim(0, 4)
    ax[0].set_xlim(0, 8)
    ax[0].set_xlabel(r'R in $[R_e]$')
    ax[0].set_ylabel(r'Q')
    ax[0].set_title("Without Bar")
    ax[1].grid()
    ax[1].set_ylim(0, 4)
    ax[1].set_xlim(0, 8)
    ax[1].set_xlabel(r'R in $[R_e]$')
    ax[1].set_ylabel(r'Q')
    ax[1].set_title("With Bar")
    fig.suptitle(r'Q depending on $\mu_0$ in [L$_{sun}$]')
    norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ticks=np.linspace(minval, maxval, 11),
                 boundaries=np.arange(minval, maxval, .1),  ax=ax.ravel().tolist())
    plt.show()


def plot_type(Qdisk, Qbar):
    tBar = ta.get_type(gBar)
    tDisk = ta.get_type(gDisk)
    types = ['0', 'a', 'ab', 'b', 'bc', 'c', 'cd', 'd']
    cmap = plt.cm.get_cmap('jet')
    maxval = len(types)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for i in range(len(Rdisk)):
        t = tDisk[i][1:]
        try:
            col = cmap(types.index(t)/maxval)
            ax[0].semilogx(Rdisk[i], Qdisk[i], '-', ms=4, color=col)
        except:
            pass
    for i in range(len(Rbar)):
        t = tBar[i][2:]
        try:
            col = cmap(types.index(t)/maxval)
            ax[1].semilogx(Rbar[i], Qbar[i], '-', ms=4, color=col)
        except:
            pass
    ax[0].grid()
    ax[0].set_ylim(0, 3)
    ax[0].set_xlim(0, 8)
    ax[0].set_xlabel(r'R in $[R_e]$')
    ax[0].set_ylabel(r'Q')
    ax[0].set_title("Without Bar")
    ax[1].grid()
    ax[1].set_ylim(0, 3)
    ax[1].set_xlim(0, 8)
    ax[1].set_xlabel(r'R in $[R_e]$')
    ax[1].set_ylabel(r'Q')
    ax[1].set_title("With Bar")
    fig.suptitle(r'Q Of Late Type Galaxies S0 to Sd')
    norm = mpl.colors.Normalize(vmin=0, vmax=maxval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ticks=np.linspace(0, maxval, 9),
                       boundaries=np.arange(0, maxval, 1),  ax=ax.ravel().tolist())
    clb.ax.set_yticklabels(['S'+t for t in types])
    plt.show()


def plot_mass(Qdisk, Qbar):
    mBar = ta.read_stellar_mass(gBar)
    mDisk = ta.read_stellar_mass(gDisk)
    cmap = plt.cm.get_cmap('jet')
    maxval = max(np.append(mBar, mDisk))
    minval = min(np.append(mBar, mDisk))
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for i in range(len(Rdisk)):
        m = mDisk[i]
        try:
            col = cmap((m-minval)/(maxval-minval))
            ax[0].semilogx(Rdisk[i], Qdisk[i], 'o', ms=4, color=col)
        except Exception as e:
            print(e)
    for i in range(len(Rbar)):
        m = mBar[i]
        try:
            col = cmap((m-minval)/(maxval-minval))
            ax[1].semilogx(Rbar[i], Qbar[i], 'o', ms=4, color=col)
        except Exception as e:
            print(e)
    ax[0].grid()
    ax[0].set_ylim(0, 3)
    ax[0].set_xlim(0, 8)
    ax[0].set_xlabel(r'R in $[R_e]$')
    ax[0].set_ylabel(r'Q')
    ax[0].set_title("Without Bar")
    ax[1].grid()
    ax[1].set_ylim(0, 3)
    ax[1].set_xlim(0, 8)
    ax[1].set_xlabel(r'R in $[R_e]$')
    ax[1].set_ylabel(r'Q')
    ax[1].set_title("With Bar")
    fig.suptitle(r'Q With Galaxy Mass in $[M_{sun}*10^9]$')
    norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    plt.show()
    plt.close()


def plot_hist(typ=None):
    global Qbar_median
    global Qbar_mean
    global Qdisk_median
    global Qdisk_mean
    if typ != None:
        nQbar_med, nQbar_mea = [], []
        nQdisk_med, nQdisk_mea = [], []
        tBar = ta.get_type(gBar)
        tDisk = ta.get_type(gDisk)
        for i in range(len(tBar)):
            if tBar[i] in typ:
                print(len(Qbar_median))
                nQbar_med.append(Qbar_median[i])
                nQbar_mea.append(Qbar_mean[i])
        for i in range(len(tDisk)):
            if tDisk[i] in typ:
                nQdisk_med.append(Qdisk_median[i])
                nQdisk_mea.append(Qdisk_mean[i])
        Qbar_median = nQbar_med
        Qbar_mean = nQbar_mea
        Qdisk_median = nQdisk_med
        Qdisk_mean = nQdisk_mea

    bins = np.linspace(0, 3, 12)

    w1 = np.ones_like(Qbar_median)/len(Qbar_median)
    w2 = np.ones_like(Qdisk_median)/len(Qdisk_median)
    n, bins, patches = plt.hist(
        [Qbar_median, Qdisk_median], bins, label=['Bar', 'No Bar'], weights=[w1, w2])
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    plt.title("Distribution of the median of Q for Sc, Scd and Sd Type Galaxies")
    # plt.errorbar(bin_centers, n, yerr=Qbar_std, fmt='r.')
    plt.legend()
    plt.show()

    n, bins, patches = plt.hist(
        [Qbar_mean, Qdisk_mean], bins, label=['Bar', 'No Bar'], weights=[w1, w2])
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    # plt.errorbar(bin_centers, n, yerr=Qdisk_std, fmt='r.')
    plt.title("Distribution of the mean of Q for Sc, Scd and Sd Type Galaxies")
    plt.legend()
    plt.show()
    plt.close()


def plot_hist_mass(m_lower, m_upper):
    global Qbar_median
    global Qbar_mean
    global Qdisk_median
    global Qdisk_mean
    nQbar_med, nQbar_mea = [], []
    nQdisk_med, nQdisk_mea = [], []
    mBar = ta.read_mass(gBar)
    mDisk = ta.read_mass(gDisk)
    for i in range(len(gBar)):
        if mBar[i] >= m_lower and mBar[i] <= m_upper:
            nQbar_med.append(Qbar_median[i])
            nQbar_mea.append(Qbar_mean[i])
    for i in range(len(gDisk)):
        if mDisk[i] >= m_lower and mDisk[i] <= m_upper:
            nQdisk_med.append(Qdisk_median[i])
            nQdisk_mea.append(Qdisk_mean[i])
    bins = np.linspace(0, 3, 12)

    w1 = np.ones_like(nQbar_med)/len(nQbar_med)
    w2 = np.ones_like(nQdisk_med)/len(nQdisk_med)
    n, bins, patches = plt.hist(
        [nQbar_med, nQdisk_med], bins, label=['Bar', 'No Bar'], weights=[w1, w2])
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    plt.title(
        "Distribution of the median of Q for Galaxies with masses" +
        "{:2.1g} - {:2.1g}".format(m_lower, m_upper))
    # plt.errorbar(bin_centers, n, yerr=Qbar_std, fmt='r.')
    plt.legend()
    plt.show()

    n, bins, patches = plt.hist(
        [nQbar_mea, nQdisk_mea], bins, label=['Bar', 'No Bar'], weights=[w1, w2])
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    # plt.errorbar(bin_centers, n, yerr=Qdisk_std, fmt='r.')
    plt.title("Distribution of the mean of Q for Galaxies with masses" +
              "{:2.1g} - {:2.1g}".format(m_lower, m_upper))
    plt.legend()
    plt.show()
    plt.close()


def plot_error_hist():
    print(max(Qbar_std))
    print(max(Qdisk_std))
    bins = np.linspace(0, 6, 12)
    n, bins, patches = plt.hist(
        [Qbar_std, Qdisk_std], bins, label=['Bar', 'No Bar'])
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    plt.title("Distribution of the standard deviation of Q")
    plt.legend()
    plt.show()
    plt.close()


def Q_sigmaZ_scatter_profile(Qdisk, Qbar, Mbar, Mdisk):
    Mbar = Mbar * 1e-9
    Mdisk = Mdisk * 1e-9
    sigZBar = []
    sigZDisk = []
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    cmap = plt.cm.get_cmap('jet')
    maxval = max(np.append(Mbar, Mdisk))
    minval = min(np.append(Mbar, Mdisk))
    maxval = 1000
    for i in range(len(gBar)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gBar[i])
        sigmaZ = np.sqrt(G.data_cyl['vzz'])
        # discard the bulge dominated values
        sigZBar.append(sigmaZ[len(sigmaZ)-len(Rbar[i]):])
        col = cmap(Mbar[i]/maxval)
        ax[1].scatter(sigZBar[i], Qbar[i], s=4, color=col)
    for i in range(len(gDisk)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gDisk[i])
        sigmaZ = np.sqrt(G.data_cyl['vzz'])
        # discard the bulge dominated values
        sigZDisk.append(sigmaZ[len(sigmaZ)-len(Rdisk[i]):])
        col = cmap(Mbar[i]/maxval)
        ax[0].scatter(sigZDisk[i], Qdisk[i], s=4, color=col)
    plt.suptitle(r'Q vs $\sigma_Z$')
    ax[0].set_title("Without Bar")
    ax[1].set_title("With Bar")
    ax[0].set_ylabel("Q")
    ax[1].set_ylabel("Q")
    ax[0].set_xlabel(r'$\sigma_Z$ in [km/s]')
    ax[1].set_xlabel(r'$\sigma_Z$ in [km/s]')
    ax[0].set_xlim(10, 250)
    ax[1].set_xlim(10, 250)
    ax[0].set_ylim(0.06, 20)
    ax[1].set_ylim(0.06, 20)
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].grid()
    ax[1].grid()
    norm = mpl.colors.LogNorm(vmin=minval, vmax=maxval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    plt.show()
    plt.close()


def Q_sigmaZ_scatter(method, Qdisk, Qbar, Mbar, Mdisk):
    Mbar = Mbar * 1e-9
    Mdisk = Mdisk * 1e-9
    sigZbar = []
    sigZdisk = []
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    cmap = plt.cm.get_cmap('jet')
    maxval = max(np.append(Mbar, Mdisk))
    minval = min(np.append(Mbar, Mdisk))
    maxval = 1000
    for i in range(len(gBar)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gBar[i])
        sigmaZ = np.sqrt(G.data_cyl['vzz'])
        sigmaZ = sigmaZ[len(sigmaZ)-len(Rbar[i]):]
        m_sigZ = np.mean(sigmaZ)
        col = cmap(Mbar[i]/maxval)
        ax[1].scatter(m_sigZ, Qbar[i], s=5, color=col)
    for i in range(len(gDisk)):
        G = g.Galaxy("data/sRsz_data/", "data/massmge/", gDisk[i])
        sigmaZ = np.sqrt(G.data_cyl['vzz'])
        m_sigZ = np.mean(sigmaZ)
        col = cmap(Mdisk[i]/maxval)
        ax[0].scatter(m_sigZ, Qdisk_mean[i], s=5, color=col)
    plt.suptitle(r'Q vs the '+method[0].upper() +
                 method[1:] + r' of $\sigma_z$')
    ax[0].set_title("Without Bar")
    ax[1].set_title("With Bar")
    ax[0].grid()
    ax[1].grid()
    ax[0].set_xlabel(r'$\sigma_Z$ in [km/s]')
    ax[1].set_xlabel(r'$\sigma_Z$ in [km/s]')
    ax[0].set_ylabel("Q")
    ax[1].set_ylabel("Q")
    ax[0].set_xlim(20, 250)
    ax[1].set_xlim(20, 250)
    ax[0].set_ylim(0.1, 4)
    ax[1].set_ylim(0.1, 4)
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[0].set_yscale('log')
    norm = mpl.colors.LogNorm(vmin=minval, vmax=maxval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    plt.show()
    plt.close()


def Q_rbar_scatter():
    P = ph.Photometrics("data/photometric_decomposition.fits")
    rbar = np.zeros(len(gBar))
    for i in range(len(gBar)):
        r = P.get_param(gBar[i], "RBAR", "R")
        rbar[i] = r
    rbar = rbar / np.array(REbar)
    plt.scatter(rbar, Qbar_mean)
    plt.title("Q vs Bar Length In Effectiv Radius")
    plt.ylabel("Q")
    plt.xlabel(r'Bar Length in [$R_e$]')
    plt.grid()
    plt.show()
    plt.close()


def Q_mass_scatter(method, Qdisk, Qbar):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for i in range(len(gBar)):
        ax[1].scatter(Mbar[i], Qbar[i], s=5, color='red')
    for i in range(len(gDisk)):
        ax[0].scatter(Mdisk[i], Qdisk[i], s=5, color='blue')
    ax[0].set_title("Without Bar")
    ax[1].set_title("With Bar")
    ax[0].grid()
    ax[1].grid()
    ax[0].set_xlabel(r'Mass in [$M_{Sun}$]')
    ax[1].set_xlabel(r'Mass in [$M_{Sun}$]')
    ax[0].set_ylabel("Q")
    ax[1].set_ylabel("Q")
    # ax[0].set_xlim(20, 250)
    # ax[1].set_xlim(20, 250)
    ax[0].set_ylim(0.1, 10)
    ax[1].set_ylim(0.1, 10)
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[0].set_yscale('log')
    plt.show()
    plt.close()


if __name__ == "__main__":
    plot_mass(Qdisk, Qbar)
    # Q_sigmaZ_scatter_profile(Qdisk, Qbar, Mbar, Mdisk)
    # Q_sigmaZ_scatter('mean', Qdisk_mean, Qbar_mean, Mbar, Mdisk)
    # Q_sigmaZ_scatter('median', Qdisk, Qbar_median, Mbar, Mdisk)
