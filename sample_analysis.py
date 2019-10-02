import numpy as np
import matplotlib.pyplot as plt
import tool_analyze as ta
from astropy.io import fits

# get data
from prepare_data import *


def bin_type(names):
    """Count the number of galaxies of a Hubble type in bins for
    E, S0, Sa, Sb, Sc, Sd
    Args:
        names (list): names of the galaxies.
    Returns:
        np.array: bins of types.
    """
    bins = np.zeros(6)
    hTypes = ta.get_type(names)
    for i in range(len(names)):
        if hTypes[i][0] == 'E':
            bins[0] += 1
        elif hTypes[i][1] == '0':
            bins[1] += 1
        elif hTypes[i][1] == 'a':
            bins[2] += 1
        elif hTypes[i][1] == 'b':
            bins[3] += 1
        elif hTypes[i][1] == 'c':
            bins[4] += 1
        elif hTypes[i][1] == 'd':
            bins[5] += 1
    return bins


def bin_mass(names):
    """Count the number of galaxies in a stellar mass range.
    9, 9.5, 10, 10.5, 11, 11.5
    Args:
        names (list)
    Returns:
        np.array: bins of mass"""
    M = ta.read_stellar_mass(names)
    bins = np.zeros(6)
    for i in range(len(names)):
        if M[i] < 9.25:
            bins[0] += 1
        elif M[i] < 9.75:
            bins[1] += 1
        elif M[i] < 10.25:
            bins[2] += 1
        elif M[i] < 10.75:
            bins[3] += 1
        elif M[i] < 11.25:
            bins[4] += 1
        else:
            bins[5] += 1
    return bins


def bin_redshift(names):
    """Bin the redshift"""
    z = ta.get_redshift(names)*1e2
    bins = np.zeros(6)
    for i in range(len(names)):
        if z[i] < 0.75:
            bins[0] += 1
        elif z[i] < 1.25:
            bins[1] += 1
        elif z[i] < 1.75:
            bins[2] += 1
        elif z[i] < 2.25:
            bins[3] += 1
        elif z[i] < 2.75:
            bins[4] += 1
        else:
            bins[5] += 1
    return bins


def MSbin_type():
    hdul = fits.open("data/sample/CALIFA_DR3_sample.fits")
    ids = hdul[1].data['califaid']
    hdulMS = fits.open("data/sample/CALIFA_2_MS_class.fits")
    MSids = hdulMS[1].data['CALIFAID']
    MStypes = hdulMS[1].data['hubtyp']+hdulMS[1].data['hubsubtyp']
    MStypes = dict(zip(MSids, MStypes))
    hdulES = fits.open("data/sample/CALIFA_2_ES_class.fits")
    ESids = hdulES[1].data['califaid']
    EStypes = hdulES[1].data['hubtyp']+hdulES[1].data['hubsubtyp']
    EStypes = dict(zip(ESids, EStypes))
    bins = np.zeros(6)
    for j in range(len(ids)):
        t = None
        if ids[j] not in MSids:
            if ids[j] not in ESids:
                print(
                    "ID {} not in mother and extension-sample".format(ids[j]))
            else:
                t = EStypes[ids[j]]
        else:
            t = MStypes[ids[j]]
        if t != None:
            if t[0] == 'E':
                bins[0] += 1
            elif t[1] == '0':
                bins[1] += 1
            elif t[1] == 'a':
                bins[2] += 1
            elif t[1] == 'b':
                bins[3] += 1
            elif t[1] == 'c':
                bins[4] += 1
            elif t[1] == 'd':
                bins[5] += 1
    return bins


def MSbin_mass():
    hdul = fits.open("data/sample/CALIFA_DR3_sample.fits")
    ids = hdul[1].data['califaid']
    hdulMS = fits.open("data/sample/CALIFA_6_MS_GC_Mstar_optical.fits")
    MSids = hdulMS[1].data['CALIFAID']
    MSmass = hdulMS[1].data['mstar']
    MSmass = dict(zip(MSids, MSmass))
    hdulES = fits.open("data/sample/CALIFA_9_ES_SDSS_Mstar_optical.fits")
    ESids = hdulES[1].data['CALIFAID']
    ESmass = hdulES[1].data['mstar']
    ESmass = dict(zip(ESids, ESmass))
    bins = np.zeros(6)
    for i in range(len(ids)):
        m = None
        if ids[i] not in MSids:
            if ids[i] not in ESids:
                print("ID {} not in mother or extension sample".format(ids[i]))
            else:
                m = ESmass[ids[i]]
        else:
            m = MSmass[ids[i]]
        if m < 9.25:
            bins[0] += 1
        elif m < 9.75:
            bins[1] += 1
        elif m < 10.25:
            bins[2] += 1
        elif m < 10.75:
            bins[3] += 1
        elif m < 11.25:
            bins[4] += 1
        else:
            bins[5] += 1
    return bins


def MSbin_redshift():
    hdul = fits.open("data/sample/CALIFA_DR3_sample.fits")
    ids = hdul[1].data['califaid']
    hdulMS = fits.open("data/sample/CALIFA_1_MS_basic.fits")
    MSids = hdulMS[1].data['CALIFAID']
    MSz = hdulMS[1].data['redshift'] * 1e2
    MSz = dict(zip(MSids, MSz))
    hdulES = fits.open("data/sample/CALIFA_1_ES_basic.fits")
    ESids = hdulES[1].data['califaid']
    ESz = hdulES[1].data['redshift'] * 1e2
    ESz = dict(zip(ESids, ESz))
    bins = np.zeros(6)
    for i in range(len(ids)):
        z = None
        if ids[i] not in MSids:
            if ids[i] not in ESids:
                print("ID {} not in mother or extension sample".format(ids[i]))
            else:
                z = ESz[ids[i]]
        else:
            z = MSz[ids[i]]
        if z < 0.75:
            bins[0] += 1
        elif z < 1.25:
            bins[1] += 1
        elif z < 1.75:
            bins[2] += 1
        elif z < 2.25:
            bins[3] += 1
        elif z < 2.75:
            bins[4] += 1
        else:
            bins[5] += 1
    return bins


def plot_overview():
    sTypes = ['E', 'S0', 'Sa', 'Sb', 'Sc', 'Sd']
    nTypes = bin_type(names)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
    ax[0].bar(sTypes, nTypes)
    ax[0].set_ylabel('Number')
    ax[0].set_title('Hubble Type')
    sMass = ['9', '9.5', '10', '10.5', '11', '11.5']
    nMass = bin_mass(names)
    ax[1].bar(sMass, nMass)
    ax[1].set_title(r'log(M/M$_\odot$)')
    sRedshift = ['0.5', '1', '1.5', '2', '2.5', '3']
    nRedshift = bin_redshift(names)
    ax[2].bar(sRedshift, nRedshift)
    ax[2].set_title(r'Redshift (x10$^{-2}$)')
    plt.savefig("figures/report/sample_overview.pdf")
    plt.show()
    plt.close()


def plot_bar_sample():
    sTypes = ['E', 'S0', 'Sa', 'Sb', 'Sc', 'Sd']
    nTypes = bin_type(gBar)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
    ax[0].bar(sTypes, nTypes)
    ax[0].set_ylabel('Number')
    ax[0].set_title('Hubble Type')
    sMass = ['9', '9.5', '10', '10.5', '11', '11.5']
    nMass = bin_mass(gBar)
    ax[1].bar(sMass, nMass)
    ax[1].set_title(r'log(M/M$_\odot$)')
    sRedshift = ['0.5', '1', '1.5', '2', '2.5', '3']
    nRedshift = bin_redshift(gBar)
    ax[2].bar(sRedshift, nRedshift)
    ax[2].set_title(r'Redshift (x10$^{-2}$)')
    plt.savefig("figures/report/sample_bar.pdf")
    plt.show()
    plt.close()


def plot_disk_sample():
    sTypes = ['E', 'S0', 'Sa', 'Sb', 'Sc', 'Sd']
    nTypes = bin_type(gDisk)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
    ax[0].bar(sTypes, nTypes)
    ax[0].set_ylabel('Number')
    ax[0].set_title('Hubble Type')
    sMass = ['9', '9.5', '10', '10.5', '11', '11.5']
    nMass = bin_mass(gDisk)
    ax[1].bar(sMass, nMass)
    ax[1].set_title(r'log(M/M$_\odot$)')
    sRedshift = ['0.5', '1', '1.5', '2', '2.5', '3']
    nRedshift = bin_redshift(gDisk)
    ax[2].bar(sRedshift, nRedshift)
    ax[2].set_title(r'Redshift (x10$^{-2}$)')
    plt.savefig("figures/report/sample_disk.pdf")
    plt.show()
    plt.close()


def plot_joint_sample():
    sTypes = ['E', 'S0', 'Sa', 'Sb', 'Sc', 'Sd']
    nDTypes = bin_type(gDisk)
    nBTypes = bin_type(gBar)
    nETypes = bin_type(gElli)
    MStype = MSbin_type()
    fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.subplots_adjust(wspace=0)
    ax[0].bar(sTypes, nBTypes)
    ax[0].bar(sTypes, nDTypes, bottom=nBTypes)
    ax[0].bar(sTypes, nETypes, bottom=nDTypes+nBTypes)
    ax[0].bar(sTypes, MStype, fill=False, linestyle='--')
    ax[0].set_ylabel('Number')
    ax[0].set_title('Hubble Type')
    sMass = ['9', '9.5', '10', '10.5', '11', '11.5']
    nDMass = bin_mass(gDisk)
    nBMass = bin_mass(gBar)
    nEMass = bin_mass(gElli)
    MSmass = MSbin_mass()
    ax[1].bar(sMass, nBMass)
    ax[1].bar(sMass, nDMass, bottom=nBMass)
    ax[1].bar(sMass, nEMass, bottom=nBMass+nDMass)
    ax[1].bar(sMass, MSmass, fill=False, linestyle='--')
    ax[1].set_title(r'log(M/M$_\odot$)')
    sRedshift = ['0.5', '1', '1.5', '2', '2.5', '3']
    nDRedshift = bin_redshift(gDisk)
    nBRedshift = bin_redshift(gBar)
    nERedshift = bin_redshift(gElli)
    MSredshift = MSbin_redshift()
    ax[2].bar(sRedshift, nBRedshift)
    ax[2].bar(sRedshift, nDRedshift, bottom=nBRedshift)
    ax[2].bar(sRedshift, nERedshift, bottom=nBRedshift+nDRedshift)
    ax[2].bar(sRedshift, MSredshift, fill=False, linestyle='--')
    ax[2].set_title(r'Redshift (x10$^{-2}$)')
    ax[2].legend(['Barred', 'Non-Barred', 'Ellipticals', 'CALIFA-Sample'])
    plt.savefig("figures/report/sample_stacked.pdf")
    plt.show()
    plt.close()


if __name__ == "__main__":
    # plot_overview()
    # plot_bar_sample()
    # plot_disk_sample()
    plot_joint_sample()
    # MSbin_type()
    # MSbin_mass()
