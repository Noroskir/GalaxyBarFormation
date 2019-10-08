import numpy as np
import tool_analyze as ta
import matplotlib.pyplot as plt
import matplotlib.image as img

# get data
from prepare_data import *


def filter_mass(names, bound=(10, 10.5)):
    M = ta.read_stellar_mass(names)
    g = []
    for i in range(len(names)):
        if M[i] >= bound[0] and M[i] <= bound[1]:
            g.append(names[i])
    return g


def filter_Q(names, Q, bound=(1, 2)):
    g = []
    for i in range(len(names)):
        if Q[i] >= bound[0] and Q[i] <= bound[1]:
            g.append(names[i])
    return g


def filter_sigZ(names, bound=(50, 100)):
    sigZ, e_sigZ = ta.get_sigma_Z(names)
    sigZmedian, sigZmean, sigZstd = ta.get_median_mean_std(sigZ)
    g = []
    for i in range(len(names)):
        if sigZmean[i] >= bound[0] and sigZmean[i] <= bound[1]:
            g.append(names[i])
    return g


def plot_QX(nBar, nDisk):
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)
    for i in range(len(nDisk)):
        ind = gDisk.index(nDisk[i])
        ax[0].plot(Rdisk[ind], Qdisk[ind], color='red')
        ax[0].plot(Rdisk[ind], Xdisk[ind], color='black')
    for i in range(len(nBar)):
        ind = gBar.index(nBar[i])
        ax[1].plot(Rbar[ind], Qbar[ind], color='red')
        ax[1].plot(Rbar[ind], Xbar[ind], color='black')
    ax[0].grid()
    ax[0].set_title('Disks')
    ax[0].set_ylabel("Q, X")
    ax[1].grid()
    ax[1].set_title('Bars')
    plt.show()
    plt.close()


def plot_images(nBar, nDisk):
    """"""
    fig, ax = plt.subplots(
        2, max([len(nBar), len(nDisk)]), figsize=(15, 10), sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    ax[0][-1].remove()
    # ax[0][-2].remove()
    for i in range(len(nBar)):
        im = img.imread("figures/report/images/bars/{}.jpeg".format(nBar[i]))
        ax[0][i].imshow(im)
        ax[0][i].get_xaxis().set_visible(False)
        ax[0][i].get_yaxis().set_visible(False)
    for i in range(len(nDisk)):
        im = img.imread("figures/report/images/disks/{}.jpeg".format(nDisk[i]))
        ax[1][i].imshow(im)
        ax[1][i].get_xaxis().set_visible(False)
        ax[1][i].get_yaxis().set_visible(False)
    plt.show()
    plt.close()


if __name__ == "__main__":
    # TODO: remove after investigation
    # BEGIN INVESTIGATION
    simBar = filter_Q(gBar, Qbar_mean, bound=(10, 100))
    simBar = filter_mass(simBar, bound=(1, 100))
    simBar = filter_sigZ(simBar, bound=(1, 1000))
    print('Mean Q over 10 (BAR): ')
    print(simBar)
    simDisk = filter_Q(gDisk, Qdisk_mean, bound=(10, 100))
    simDisk = filter_mass(simDisk, bound=(1, 100))
    simDisk = filter_sigZ(simDisk, bound=(1, 1000))
    print('Mean Q over 10 (Disk):')
    print(simDisk)
    # END INVESTIGATION

    # simBar = filter_Q(gBar, Qbar_mean, bound=(1.5, 2.5))
    # #simBar = filter_Q(gBar, Qbar_mean, bound=(.5, 1))
    # simBar = filter_mass(simBar, bound=(10, 10.5))
    # simBar = filter_sigZ(simBar, bound=(50, 100))
    # simDisk = filter_Q(gDisk, Qdisk_mean, bound=(1.5, 2.5))
    # #simDisk = filter_Q(gDisk, Qdisk_mean, bound=(.5, 1))
    # simDisk = filter_mass(simDisk, bound=(10, 10.5))
    # simDisk = filter_sigZ(simDisk, bound=(50, 100))
    # print("Bars: ", simBar)
    # xBar = []
    # for i in range(len(simBar)):
    #     xBar.append(Xbar_mean[gBar.index(simBar[i])])
    # xDisk = []
    # for i in range(len(simDisk)):
    #     xDisk.append(Xdisk_mean[gDisk.index(simDisk[i])])
    # print("X:", xBar)
    # print("Disks: ", simDisk)
    # print("X:", xDisk)
    # # simDisk.remove('NGC4644')
    # plot_images(simBar, simDisk)
