import numpy as np
import galaxy as g
import tool_analyze as ta
from prepare_data import *


f = open("data/table.txt", 'w')
f.write("{:15s}{:8s}{:7s}{:8s}{:9s}".format('Name', 'HubTyp', 'Bar', 'Q', 'X'))
f.write("{:10s}{:10s}{:10s}".format('sigma_z', 'M_stellar', 'M_dark'))
f.write("{:10}{:10s}".format('M_total', 'lambda_re'))
f.write("\n")
f.write("{:15s}{:8s}{:7s}{:8s}{:9s}".format('', '', '', '', ''))
f.write("{:10s}{:<10s}{:10s}".format('[km/s]', '[', 'log(M/M_sun)'))
f.write("{:>6s}{:10s}".format(']', ''))
f.write("\n")
hType = ta.get_type(gBar)
lambReBar = ta.read_lambda_re(gBar)
for i in range(len(gBar)):
    G = g.Galaxy(gBar[i])
    f.write("{:15s}".format(G.name))
    f.write("{:8s}".format(hType[i]))
    f.write("{:6s}".format(barredness[gBar[i]]))
    f.write("{:6.3f}".format(Qbar_mean[i]))
    f.write("{:8.3f}".format(Xbar_mean[i]))
    f.write("{:10.3f}".format(sigZbar_mean[i]))
    f.write("{:10.3f}".format(np.log10(G.get_stellar_mass())))
    f.write("{:10.3f}".format(np.log10(G.get_dark_matter())))
    f.write("{:10.3f}".format(
        np.log10(G.get_stellar_mass()+G.get_dark_matter())))
    f.write("{:10.3f}".format(lambReBar[i]))
    f.write("\n")

hType = ta.get_type(gDisk)
lambReDisk = ta.read_lambda_re(gDisk)
for i in range(len(gDisk)):
    G = g.Galaxy(gDisk[i])
    f.write("{:15s}".format(G.name))
    f.write("{:8s}".format(hType[i]))
    f.write("{:6s}".format(barredness[gDisk[i]]))
    f.write("{:6.3f}".format(Qdisk_mean[i]))
    f.write("{:8.3f}".format(Xdisk_mean[i]))
    f.write("{:10.3f}".format(sigZdisk_mean[i]))
    f.write("{:10.3f}".format(np.log10(G.get_stellar_mass())))
    f.write("{:10.3f}".format(np.log10(G.get_dark_matter())))
    f.write("{:10.3f}".format(
        np.log10(G.get_stellar_mass()+G.get_dark_matter())))
    f.write("{:10.3f}".format(lambReDisk[i]))
    f.write("\n")

hType = ta.get_type(gElli)
lambReElli = ta.read_lambda_re(gElli)
for i in range(len(gElli)):
    G = g.Galaxy(gElli[i])
    f.write("{:15s}".format(G.name))
    f.write("{:8s}".format(hType[i]))
    f.write("{:6s}".format(barredness[gElli[i]]))
    f.write("{:6.3f}".format(Qelli_mean[i]))
    f.write("{:8.3f}".format(Xelli_mean[i]))
    f.write("{:10.3f}".format(sigZelli_mean[i]))
    f.write("{:10.3f}".format(np.log10(G.get_stellar_mass())))
    f.write("{:10.3f}".format(np.log10(G.get_dark_matter())))
    f.write("{:10.3f}".format(
        np.log10(G.get_stellar_mass()+G.get_dark_matter())))
    f.write("{:10.3f}".format(lambReElli[i]))
    f.write("\n")
