import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tool_analyze as ta
import galaxy as g
import photometrics as pm
from scipy import stats

# get the data
from prepare_data import *


def plot_error_hist_X():
    bins = np.linspace(0, 1, 10)
    n, bins, patches = plt.hist([Xbar_std, Xdisk_std, Xelli_std], bins, label=[
        'Bar', 'No Bar', 'Ellipticals'])
    plt.title("Distribution of the standard deviation of X")
    plt.legend()
    plt.show()
    plt.close()


def plot_error_hist_X():
    bins = np.linspace(0, 1, 10)
    n, bins, patches = plt.hist([Qbar_std, Qdisk_std, Qelli_std], bins, label=[
        'Bar', 'No Bar', 'Ellipticals'])
    plt.title("Distribution of the standard deviation of Q")
    plt.legend()
    plt.show()
    plt.close()


def plot_Q_profile(name):
    G = g.Galaxy(name)
    Q, eQ = G.get_toomre()
    R = G.data['R']
    plt.plot(R, Q, color='red')
    plt.errorbar(R, Q, yerr=eQ, fmt='.-', color='red', capsize=2)
    plt.grid()
    plt.title("Toomre Parameter Q " + name)
    plt.xlabel(r'Radius [arcsec]')
    plt.ylabel(r'Q')
    plt.savefig("figures/report/RadialProfile/Q_{}.pdf".format(name))
    plt.show()
    plt.close()


def plot_X_profile(name):
    G = g.Galaxy(name)
    R = G.data['R']
    X = G.swing_amplif_X()
    plt.plot(R, X)
    plt.title(r'Swing Amplification Parameter X '+name)
    plt.ylabel(r'X')
    plt.xlabel(r'R [arcsec]')
    plt.grid()
    plt.savefig("figures/report/RadialProfile/X_{}.pdf".format(name))
    plt.show()
    plt.close()


def plot_Q_X_profile(name):
    G = g.Galaxy(name)
    R = G.data['R']
    X = G.swing_amplif_X()
    Q, eQ = G.get_toomre()
    plt.plot(R, Q, color='red', label='Q')
    plt.errorbar(R, Q, yerr=eQ, fmt='.-', color='red', capsize=2)
    plt.plot(R, X, color='black', label='X')
    plt.title(r'Q and X '+name)
    plt.xlabel(r'R [arcsec]')
    plt.ylabel(r'Q and X')
    plt.grid()
    plt.legend()
    plt.savefig("figures/report/RadialProfile/QX_{}.pdf".format(name))
    plt.show()
    plt.close()


def plot_sigR_profile(name):
    G = g.Galaxy(name)
    R = G.data['R']
    sigR = np.sqrt(G.data['vRR'])
    e_sigR = 0.5 * 1/sigR * G.data['e_vRR']
    plt.plot(R, sigR, color='red')
    plt.errorbar(R, sigR, yerr=e_sigR, fmt='.-',
                 color='red', capsize=2)
    plt.grid()
    plt.title(r"Radial Velocity Dispersion $\sigma_R$ " + name)
    plt.xlabel(r'Radius [arcsec]')
    plt.ylabel(r'$\sigma_R$ [km/s]')
    plt.savefig("figures/report/RadialProfile/sigmaR_{}.pdf".format(name))
    plt.show()


def plot_vcirc_profile(name):
    G = g.Galaxy(name)
    R = G.data['R']
    vcirc = G.velocity_vcirc()
    plt.plot(R, vcirc, color='red')
    plt.grid()
    plt.title(r'Circular Velocity ' + name)
    plt.xlabel(r'Radius [arcsec]')
    plt.ylabel(r'$v_c$ [km/s]')
    plt.savefig("figures/report/RadialProfile/Vcirc_{}.pdf".format(name))
    plt.show()
    plt.close()


def plot_SMD_profile(name):
    G = g.Galaxy(name)
    R = G.data['R']
    smd_lm = G.stellar_surfMD(R)
    smd_dm = G.darkmatter_surfMD(R)
    smd = smd_lm + smd_dm
    plt.semilogy(R, smd_lm, '--', color='red', label='Stellar Mass')
    plt.semilogy(R, smd_dm, '--', color='black', label='Dark Matter')
    plt.semilogy(R, smd, color='blue', label='Total')
    plt.title(r'Surface Mass Density')
    plt.xlabel('R [arcsec]')
    plt.ylabel(r'$\Sigma$ [M$_{sun}$/pc$^2$]')
    plt.grid()
    plt.legend()
    plt.savefig("figures/report/RadialProfile/surfMD_{}.pdf".format(name))
    plt.show()
    plt.close()


def plot_Q_total_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(Relli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(Rdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(Rbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.8)
        ax[i].set_ylim(0.2, 120)
        ax[i].set_xlabel(r'R/R$_e$')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/Combined/plot_Q_total_mass.pdf")
    plt.show()
    plt.close()


def plot_X_total_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(Relli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(Rdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(Rbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.8)
        ax[i].set_ylim(.04, 150)
        ax[i].set_xlabel(r'R/R$_e$')
        ax[i].set_ylabel(r'X')
    plt.suptitle('Swing Amplification Parameter X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/Combined/plot_X_total_mass.pdf")
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
        ax[0].semilogy(Relli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].semilogy(Rdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].semilogy(Rbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.8)
        ax[i].set_ylim(.2, 120)
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
    plt.savefig("figures/report/Combined/plot_Q_total_type.pdf")
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
        ax[0].semilogy(Relli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].semilogy(Rdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].semilogy(Rbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.5)
        ax[i].set_ylim(0.04, 150)
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
    plt.savefig("figures/report/Combined/plot_X_total_type.pdf")
    plt.show()
    plt.close()


def plot_Q_sigZ_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(sigZelli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(sigZdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(sigZbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_ylim(0.2, 120)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/Combined/plot_Q_sigZ_mass.pdf")
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
        ax[0].semilogy(sigZelli[i], Qelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].semilogy(sigZdisk[i], Qdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].semilogy(sigZbar[i], Qbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_ylim(0.2, 120)
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
    plt.savefig("figures/report/Combined/plot_Q_sigZ_type.pdf")
    plt.show()
    plt.close()


def plot_X_sigZ_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(Relli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(sigZelli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(sigZdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(sigZbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_ylim(0.04, 150)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'X')
    plt.suptitle('Swing Amplification Parameter X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/Combined/plot_X_sigZ_mass.pdf")
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
        ax[0].semilogy(sigZelli[i], Xelli[i], 'x', color=c)
    for i in range(len(Rdisk)):
        if tDisk[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tDisk[i][1])/valmax)
        ax[1].semilogy(sigZdisk[i], Xdisk[i], 'x', color=c)
    for i in range(len(Rbar)):
        if tBar[i][1] in ['1', '2', '3', '4', '5', '6', '7']:
            c = cmap(0)
        else:
            c = cmap(types.index(tBar[i][1])/valmax)
        ax[2].semilogy(sigZbar[i], Xbar[i], 'x', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_ylim(0.04, 150)
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
    plt.savefig("figures/report/Combined/plot_X_sigZ_type.pdf")
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
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_xlim(0, 3.5)
        # ax[i].set_ylim(0, 3.5)
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
    res = stats.ks_2samp(Qbar_mean, Qdisk_mean)
    print("\n"+"*"*20)
    print("KS test for mean of Q:")
    print(res)
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    bins = np.linspace(0, 6, 20)
    ax.hist(Qdisk_mean, bins=bins, label='No Bar', alpha=0.7)
    ax.margins(0.05)
    ax.hist(Qbar_mean, bins=bins, label='Bar', alpha=0.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='Mean of Q')
    plt.savefig("figures/report/hist_Q_mean.pdf")
    plt.show()
    plt.close()


def hist_Q_median():
    res = stats.ks_2samp(Qbar_median, Qdisk_median)
    print("\n"+"*"*20)
    print("KS test for median of Q: ")
    print(res)
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    bins = np.linspace(0, 6, 20)
    ax.hist(Qdisk_median, bins=bins, label='No Bar', alpha=0.7)
    ax.margins(0.05)
    ax.hist(Qbar_median, bins=bins, label='Bar', alpha=0.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='Q')
    plt.savefig("figures/report/hist_Q_median.pdf")
    plt.show()
    plt.close()


def hist_X_mean():
    res = stats.ks_2samp(Xbar_mean, Xdisk_mean)
    print("\n"+"*"*20)
    print("KS test for mean of X: ")
    print(res)
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    bins = np.linspace(0, 2, 15)
    ax.hist(Xdisk_median, bins=bins, label='No Bar', alpha=.7)
    ax.hist(Xbar_median, bins=bins, label='Bar', alpha=.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='X')
    plt.savefig("figures/report/hist_X_mean.pdf")
    plt.show()
    plt.close()


def hist_X_median():
    res = stats.ks_2samp(Xbar_median, Xdisk_median)
    print("\n"+"*"*20)
    print("KS test for median of X: ")
    print(res)
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    bins = np.linspace(0, 2, 15)
    ax.hist(Xdisk_median, bins=bins, label='No Bar', alpha=.7)
    ax.hist(Xbar_median, bins=bins, label='Bar', alpha=.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='X')
    plt.savefig("figures/report/hist_X_median.pdf")
    plt.show()
    plt.close()


def hist_Q_types_mean(hTypes, title):
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    bins = np.linspace(0, 5, 15)
    x, Qbar_m = ta.filter_type(Qbar_mean, Qbar_mean, gBar, hTypes)
    x, Qdisk_m = ta.filter_type(Qdisk_mean, Qdisk_mean, gDisk, hTypes)
    res = stats.ks_2samp(Qbar_m, Qdisk_m)
    print("\n"+"*"*20)
    print("KS test for mean of Q: ", hTypes)
    print(res)
    ax.hist(Qdisk_m, bins=bins, label='No Bar', alpha=0.7)
    ax.margins(0.05)
    ax.hist(Qbar_m, bins=bins, label='Bar', alpha=0.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies of Type ' +
           title, xlabel='Mean of Q')
    title = title.replace(' ', '_')
    plt.savefig("figures/report/hist_Q_types_"+title+".pdf")
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
        ax[0].semilogy(Melli[i], Qelli_mean[i], 'o', color=c)
    for i in range(len(Mdisk)):
        c = cmap((sigZdisk_mean[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(Mdisk[i], Qdisk_mean[i], 'o', color=c)
    for i in range(len(Mbar)):
        c = cmap((sigZbar_mean[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(Mbar[i], Qbar_mean[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(9, 11.8)
        ax[i].set_ylim(0.3, 45)
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


def scatter_Q_mass_median():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(sigZelli_mean, np.append(
        sigZbar_mean, sigZdisk_mean)))
    valmin = min(np.append(sigZelli_mean, np.append(
        sigZbar_mean, sigZdisk_mean)))
    for i in range(len(Melli)):
        c = cmap((sigZelli_mean[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(Melli[i], Qelli_median[i], 'o', color=c)
    for i in range(len(Mdisk)):
        c = cmap((sigZdisk_mean[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(Mdisk[i], Qdisk_median[i], 'o', color=c)
    for i in range(len(Mbar)):
        c = cmap((sigZbar_mean[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(Mbar[i], Qbar_median[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(9, 11.8)
        ax[i].set_ylim(0.3, 45)
        ax[i].set_xlabel(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
        ax[i].set_ylabel(r'Q')
    plt.suptitle('Toomre Parameter Q vs Stellar Mass')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'$\sigma_z$ [km/s]')
    plt.savefig("figures/report/GalaxyProp/scatter_Q_mass_median.pdf")
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
        ax[0].semilogy(Melli[i], Xelli_mean[i], 'o', color=c)
    for i in range(len(Mdisk)):
        c = cmap((sigZdisk_mean[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(Mdisk[i], Xdisk_mean[i], 'o', color=c)
    for i in range(len(Mbar)):
        c = cmap((sigZbar_mean[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(Mbar[i], Xbar_mean[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(9, 11.8)
        ax[i].set_ylim(0.2, 1e2)
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
        ax[0].semilogy(sigZelli_mean[i], Qelli_mean[i], 'o', color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(sigZdisk_mean[i], Qdisk_mean[i], 'o', color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(sigZbar_mean[i], Qbar_mean[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].set_xlim(25, 280)
        ax[i].set_ylim(0.3, 45)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'Q')
        ax[i].grid()
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_Q_sigZ_mean.pdf")
    plt.show()
    plt.close()


def scatter_X_sigZ_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(gElli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(sigZelli_mean[i], Xelli_mean[i], 'o', color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(sigZdisk_mean[i], Xdisk_mean[i], 'o', color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(sigZbar_mean[i], Xbar_mean[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(25, 280)
        ax[i].set_ylim(0.2, 1e2)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
        ax[i].set_ylabel(r'X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_X_sigZ_mean.pdf")
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
        ax[0].semilogy(lambReElli[i], Qelli_mean[i], 'o', color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(lambReDisk[i], Qdisk_mean[i], 'o', color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(lambReBar[i], Qbar_mean[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, .9)
        ax[i].set_ylim(0.3, 45)
        ax[i].set_xlabel(r'$\lambda_{Re}$')
        ax[i].set_ylabel(r'Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_Q_lambRe_mean.pdf")
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
        ax[0].semilogy(lambReElli[i], Xelli_mean[i], 'o', color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(lambReDisk[i], Xdisk_mean[i], 'o', color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(lambReBar[i], Xbar_mean[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, .9)
        ax[i].set_ylim(0.2, 1e2)
        ax[i].set_xlabel(r'$\lambda_{Re}$')
        ax[i].set_ylabel(r'X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_X_lambRe_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_RBar_mean():
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
    plt.savefig("figures/report/BarProp/scatter_Q_Rbar_mean.pdf")
    plt.show()
    plt.close()


def scatter_X_RBar_mean():
    rbar = np.zeros(len(gBar))
    for i in range(len(gBar)):
        r = P.get_param(gBar[i], "RBAR", "R")
        rbar[i] = r
    rbar = rbar / np.array(REbar)
    plt.scatter(rbar, Xbar_mean)
    plt.title("X vs Bar Length In Effectiv Radius")
    plt.ylabel("X")
    plt.xlabel(r'Bar Length in [$R_e$]')
    plt.grid()
    plt.savefig("figures/report/BarProp/scatter_X_Rbar_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_e_mean():
    ebar = np.zeros(len(gBar))
    for i in range(len(gBar)):
        q = P.get_param(gBar[i], 'BABAR', 'R')
        e = 1-q
        ebar[i] = e
    plt.scatter(ebar, Qbar_mean)
    plt.title("Q vs Bar Ellipticity")
    plt.ylabel('Q')
    plt.xlabel(r'$\epsilon$')
    plt.ylim(0, 5)
    plt.grid()
    plt.savefig("figures/report/BarProp/scatter_Q_e_mean.pdf")
    plt.show()
    plt.close()


def scatter_X_e_mean():
    ebar = np.zeros(len(gBar))
    for i in range(len(gBar)):
        q = P.get_param(gBar[i], 'BABAR', 'R')
        e = 1-q
        ebar[i] = e
    plt.scatter(ebar, Xbar_mean)
    plt.title("X vs Bar Ellipticity")
    plt.ylabel('X')
    plt.xlabel(r'$\epsilon$')
    plt.ylim(0, 5)
    plt.grid()
    plt.savefig("figures/report/BarProp/scatter_X_e_mean.pdf")
    plt.show()
    plt.close()


def check_mass():
    MAbar = []
    for i in range(len(gBar)):
        G = g.Galaxy(gBar[i])
        MAbar.append(np.log10(G.get_stellar_mass()))
    plt.scatter(Mbar, MAbar)
    plt.xlabel('Stellar Mass')
    plt.ylabel('Dynamical Mass')
    plt.show()
    plt.close()


if __name__ == '__main__':
    n = 'NGC6278'
    # plot_error_hist()
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
    # hist_Q_types_mean(['Sa', 'Sab', 'Sb'], 'Sa Sab Sb')
    # hist_Q_types_mean(['Sc', 'Scd', 'Sd'], 'Sc Scd Sd')
    # scatter_Q_mass_mean()
    # scatter_X_mass_mean()
    # scatter_Q_sigZ_mean()
    # scatter_X_sigZ_mean()
    # scatter_Q_lambRe_mean()
    # scatter_X_lambRe_mean()
    # scatter_Q_RBar_mean()
    # scatter_X_RBar_mean()
    scatter_Q_e_mean()
    scatter_X_e_mean()
    # scatter_Q_mass_median()
    # check_mass()
