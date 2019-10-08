import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tool_analyze as ta
import galaxy as g
import gaussFilter as gf
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
    X = G.swing_amplif_X(R)
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
    X = G.swing_amplif_X(R)
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
    vcirc = G.velocity_vcirc(R)
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.8)
        ax[i].set_xlabel(r'R/R$_e$')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.8)
        ax[i].set_xlabel(r'R/R$_e$')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.8)
        ax[i].set_xlabel(r'R/R$_e$')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(-.1, 3.5)
        ax[i].set_xlabel(r'R/R$_e$')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(5, 380)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    nbins = 5
    step = (valmax-valmin) / nbins
    mbins = [valmin, 10.25, 10.75, 12]
    RSbar, QSbar, = ta.mass_bin(iRbar, iQbar, Mbar, mbins)
    RSdisk, QSdisk, = ta.mass_bin(iRdisk, iQdisk, Mdisk, mbins)
    RSelli, QSelli, = ta.mass_bin(iRelli, iQelli, Melli, mbins)
    step = 0.05
    width = 0.05
    for i in range(len(gElli)):
        c2 = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(Relli[i], Qelli[i], 'x', color=c2, alpha=0.8)
    for i in range(len(RSelli)):
        c = cmap(((mbins[i]+mbins[i+1])/2-valmin)/(valmax-valmin))
        R, Q = ta.stack_curves(RSelli[i], QSelli[i], step, threshold=4)
        R, Q = gf.smooth_data(R, Q, width)
        ax[0].semilogy(R, Q, lw=2, color=c)
    for i in range(len(gDisk)):
        c2 = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(Rdisk[i], Qdisk[i], 'x', color=c2, alpha=0.8)
    for i in range(len(RSdisk)):
        c = cmap(((mbins[i]+mbins[i+1])/2-valmin)/(valmax-valmin))
        R, Q = ta.stack_curves(RSdisk[i], QSdisk[i], step, threshold=4)
        R, Q = gf.smooth_data(R, Q, width)
        ax[1].semilogy(R, Q, lw=2, color=c)
    for i in range(len(Rbar)):
        c2 = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(Rbar[i], Qbar[i], 'x', color=c2, alpha=0.8)
    for i in range(len(RSbar)):
        if len(QSbar[i]) == 0:
            print('Empty bin: {}'.format(i))
        c = cmap(((mbins[i]+mbins[i+1])/2-valmin)/(valmax-valmin))
        R, Q = ta.stack_curves(RSbar[i], QSbar[i], step, threshold=1)
        R, Q = gf.smooth_data(R, Q, width)
        ax[2].semilogy(R, Q, lw=2, color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, 3.45)
        ax[i].set_xlabel(r'R/R$_e$')
    plt.suptitle('Toomre Parameter Q')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/stack_Q_profile_mass.pdf")
    plt.show()
    plt.close()


def stack_X_profile_mass():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    nbins = 5
    step = (valmax-valmin) / nbins
    mbins = [valmin, 10.25, 10.75, 12]
    RSbar, XSbar, = ta.mass_bin(iRbar, iXbar, Mbar, mbins)
    RSdisk, XSdisk, = ta.mass_bin(iRdisk, iXdisk, Mdisk, mbins)
    RSelli, XSelli, = ta.mass_bin(iRelli, iXelli, Melli, mbins)
    step = 0.05
    for i in range(len(gElli)):
        c2 = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(Relli[i], Xelli[i], 'x', color=c2, alpha=0.8)
    for i in range(len(RSelli)):
        c = cmap(((mbins[i]+mbins[i+1])/2-valmin)/(valmax-valmin))
        R, X = ta.stack_curves(RSelli[i], XSelli[i], step, threshold=4)
        ax[0].semilogy(R, X, lw=2, color=c)
    for i in range(len(gDisk)):
        c2 = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].semilogy(Rdisk[i], Xdisk[i], 'x', color=c2, alpha=0.8)
    for i in range(len(RSdisk)):
        c = cmap(((mbins[i]+mbins[i+1])/2-valmin)/(valmax-valmin))
        R, X = ta.stack_curves(RSdisk[i], XSdisk[i], step, threshold=4)
        ax[1].semilogy(R, X, lw=2, color=c)
    for i in range(len(Rbar)):
        c2 = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].semilogy(Rbar[i], Xbar[i], 'x', color=c2, alpha=0.8)
    for i in range(len(RSbar)):
        c = cmap(((mbins[i]+mbins[i+1])/2-valmin)/(valmax-valmin))
        R, X = ta.stack_curves(RSbar[i], XSbar[i], step, threshold=4)
        ax[2].semilogy(R, X, lw=2, color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, 3.45)
        ax[i].set_xlabel(r'R/R$_e$')
    plt.suptitle('Swing Amplification Parameter X')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/stack_X_profile_mass.pdf")
    plt.show()
    plt.close()


def hist_Q_mean():
    # Qbar = ta.clearly_classified(gBar, Qbar_mean, 'SB', barredness)
    # Qdisk = ta.clearly_classified(gDisk, Qdisk_mean, 'SA', barredness)
    Qbar = Qbar_mean
    Qdisk = Qdisk_mean
    res = stats.ks_2samp(Qbar, Qdisk)
    print("\n"+"*"*20)
    print("KS test for mean of Q:")
    print(res)
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    bins = np.linspace(0, 6, 20)
    ax.hist(Qdisk, bins=bins, label='No Bar', alpha=0.7)
    ax.margins(0.05)
    ax.hist(Qbar, bins=bins, label='Bar', alpha=0.7)
    ax.legend()
    ax.set(title='Barred and Non-Barred Galaxies', xlabel='Mean of Q')
    plt.savefig("figures/report/Integrated/hist_Q_mean.pdf")
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
    plt.savefig("figures/report/Integrated/hist_Q_median.pdf")
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
    plt.savefig("figures/report/Integrated/hist_X_mean.pdf")
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
    plt.savefig("figures/report/Integrated/hist_X_median.pdf")
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
    plt.savefig("figures/report/Integrated/hist_Q_types_"+title+".pdf")
    plt.show()
    plt.close()


def scatter_Q_mass_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(9, 11.8)
        ax[i].set_ylim(0.3, 45)
        ax[i].set_xlabel(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.suptitle('Toomre Parameter Q vs Stellar Mass')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'$\sigma_z$ [km/s]')
    plt.savefig("figures/report/GalaxyProp/scatter_Q_mass_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_mass_median():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(9, 11.8)
        ax[i].set_ylim(0.3, 45)
        ax[i].set_xlabel(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(9, 11.8)
        ax[i].set_ylim(0.2, 1e2)
        ax[i].set_xlabel(r'Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.suptitle('Swing Amplification Parameter X vs Stellar Mass')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label(r'$\sigma_z$ [km/s]')
    plt.savefig("figures/report/GalaxyProp/scatter_X_mass_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_sigZ_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(gElli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].semilogy(sigZelli_mean[i], Qelli_mean[i], 'o', color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        symb = 'o'
        ms = 4
        if barredness[gDisk[i]] == 'SA':
            symb = '*'
            ms = 8
        ax[1].semilogy(sigZdisk_mean[i], Qdisk_mean[i],
                       symb, ms=ms, color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        symb = 'o'
        ms = 4
        if barredness[gBar[i]] == 'SB':
            symb = '*'
            ms = 8
        ax[2].semilogy(sigZbar_mean[i], Qbar_mean[i], symb, ms=ms, color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].set_xlim(25, 280)
        ax[i].set_ylim(0.3, 45)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
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
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(25, 280)
        ax[i].set_ylim(0.2, 1e2)
        ax[i].set_xlabel(r'$\sigma_z$ [km/s]')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_X_sigZ_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_lambRe_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, .9)
        ax[i].set_ylim(0.3, 45)
        ax[i].set_xlabel(r'$\lambda_{Re}$')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_Q_lambRe_mean.pdf")
    plt.show()
    plt.close()


def scatter_X_lambRe_mean():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
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
    ax[0].set_ylabel(r'X')
    for i in [0, 1, 2]:
        ax[i].grid()
        ax[i].set_xlim(0, .9)
        ax[i].set_ylim(0.2, 1e2)
        ax[i].set_xlabel(r'$\lambda_{Re}$')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_X_lambRe_mean.pdf")
    plt.show()
    plt.close()


def scatter_Q_frac():
    fig, ax = plt.subplots(1, 3, figsize=(18, 5), sharey=True, sharex=True)
    fig.subplots_adjust(wspace=0)
    cmap = plt.cm.get_cmap('jet')
    valmax = max(np.append(Melli, np.append(Mbar, Mdisk)))
    valmin = min(np.append(Melli, np.append(Mbar, Mdisk)))
    for i in range(len(gElli)):
        c = cmap((Melli[i]-valmin)/(valmax-valmin))
        ax[0].loglog(fracElli[i], Qelli_mean[i], 'o', color=c)
    for i in range(len(gDisk)):
        c = cmap((Mdisk[i]-valmin)/(valmax-valmin))
        ax[1].loglog(fracDisk[i], Qdisk_mean[i], 'o', color=c)
    for i in range(len(gBar)):
        c = cmap((Mbar[i]-valmin)/(valmax-valmin))
        ax[2].loglog(fracBar[i], Qbar_mean[i], 'o', color=c)
    ax[0].set_title('Elliptical Galaxies')
    ax[1].set_title('Late Type Galaxies Without Bar')
    ax[2].set_title('Late Type Galaxies With Bar')
    ax[0].set_ylabel(r'Q')
    for i in [0, 1, 2]:
        ax[i].grid()
        # ax[i].set_ylim(0.2, 1e2)
        ax[i].set_xlabel(r'M$_{stellar}$/M$_{dark}$')
    norm = mpl.colors.Normalize(vmin=valmin, vmax=valmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    clb = fig.colorbar(sm, ax=ax.ravel().tolist())
    clb.set_label('Stellar Mass [log$_{10}($M/M$_{sun}$)]')
    plt.savefig("figures/report/GalaxyProp/scatter_Q_frac_mean.pdf")
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


def plot_error_mass():
    plt.plot(Mbar, Qbar_std, 'o', color='black', label='Bar')
    plt.plot(Mdisk, Qdisk_std, 'o', color='red', label='No Bar')
    plt.plot(Melli, Qelli_std, 'o', color='orange', label='Elliptical')
    plt.title('Standard deviation of Q')
    plt.ylabel('Q')
    plt.xlabel(r'Stellar Mass in [log(M/M$_{sun}$)]')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    n = 'NGC6278'
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
    # stack_X_profile_mass()
    # hist_Q_mean()
    # hist_Q_median()
    # hist_X_mean()
    # hist_X_median()
    # hist_Q_types_mean(['Sa', 'Sab', 'Sb'], 'Sa Sab Sb')
    # hist_Q_types_mean(['Sc', 'Scd', 'Sd'], 'Sc Scd Sd')
    # scatter_Q_mass_mean()
    # scatter_X_mass_mean()
    scatter_Q_frac()
    scatter_Q_sigZ_mean()
    scatter_X_sigZ_mean()
    scatter_Q_lambRe_mean()
    scatter_X_lambRe_mean()
    scatter_Q_RBar_mean()
    scatter_X_RBar_mean()
    scatter_Q_e_mean()
    scatter_X_e_mean()
    scatter_Q_mass_median()
    check_mass()
    plot_error_mass()
