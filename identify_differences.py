import numpy as np
import tool_analyze as ta

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


if __name__ == "__main__":
    simBar = filter_Q(gBar, Qbar_mean, bound=(0.3, 1))
    simBar = filter_mass(simBar, bound=(10.5, 11))
    simBar = filter_sigZ(simBar)
    simDisk = filter_Q(gDisk, Qdisk_mean, bound=(0.3, 1))
    simDisk = filter_mass(simDisk, bound=(10.5, 11))
    simDisk = filter_sigZ(simDisk)
    print("Bars: ", simBar)
    xBar = []
    for i in range(len(simBar)):
        xBar.append(Xbar_mean[gBar.index(simBar[i])])
    xDisk = []
    for i in range(len(simDisk)):
        xDisk.append(Xdisk_mean[gDisk.index(simDisk[i])])
    print("X:", xBar)
    print("Disks: ", simDisk)
    print("X:", xDisk)
