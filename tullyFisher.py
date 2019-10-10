import numpy as np
import matplotlib.pyplot as plt
import tool_analyze as ta
import galaxy as g
from prepare_data import *

vcBar = np.zeros(len(gBar))
vcDisk = np.zeros(len(gDisk))
vcElli = np.zeros(len(gElli))
for i in range(len(gBar)):
    G = g.Galaxy(gBar[i])
    R = np.linspace(G.data['R'][0], G.data['R'][-1], 300)
    vc = G.velocity_vcirc(R)
    vcBar[i] = np.max(vc)

for i in range(len(gDisk)):
    G = g.Galaxy(gDisk[i])
    R = np.linspace(G.data['R'][0], G.data['R'][-1], 300)
    vc = G.velocity_vcirc(R)
    vcDisk[i] = np.max(vc)

for i in range(len(gElli)):
    G = g.Galaxy(gDisk[i])
    R = np.linspace(G.data['R'][0], G.data['R'][-1], 300)
    vc = G.velocity_vcirc(R)
    vcElli[i] = np.max(vc)


mBar, dBar = ta.read_magnitude_dist(gBar)
mDisk, dDisk = ta.read_magnitude_dist(gDisk)
mElli, dElli = ta.read_magnitude_dist(gElli)
dBar = dBar * 1e6  # in pc
dDisk = dDisk * 1e6
dElli = dElli * 1e6


magBar = mBar + 5 - 5*np.log10(dBar)
magDisk = mDisk + 5 - 5*np.log10(dDisk)
magElli = mElli + 5 - 5*np.log10(dElli)


def bestFit(x):
    m = -7.5
    c = -4
    return (x - c)/m


xx = np.linspace(-24, -19)

plt.scatter(magBar, np.log10(vcBar), label='Barred')
plt.scatter(magDisk, np.log10(vcDisk), label='Non-Barred')
plt.scatter(magElli, np.log10(vcElli), label='Ellipticals')
plt.plot(xx, bestFit(xx), label=r'CALFIA $v_{circ}$', color='red')
plt.title("Tully-Fisher")
plt.ylabel(r'log$_{10}(v_{c,max})$')
plt.xlabel(r'Magnitude')
plt.xlim(-18.5, -24.5)
plt.legend()
plt.grid()
plt.show()
