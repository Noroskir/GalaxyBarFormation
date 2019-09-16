import unittest

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
from galaxy import Galaxy
import mge_radial_density as mrd

g = Galaxy('data/sRsz_data/', 'data/massmge/', 'NGC4210')


class TestGalaxy(unittest.TestCase):
    def test_kappa_sq(self):
        a = 10
        R = np.linspace(1, 40)
        vc = a * R**2
        kappa_sq = 6 * a**2 * R**2
        result = g.kappa_sq(R, vc)
        self.assertAlmostEqual(result[19], kappa_sq[19])
        self.assertAlmostEqual(result[2], kappa_sq[2])
        self.assertAlmostEqual(result[35], kappa_sq[35])

    def test_toomre_param(self):
        g.distance = 1
        self.assertEqual(g.distance, 1)
        R = np.linspace(1, 40)
        a = 10
        b = 200
        L = 4 * 500
        sig = 50
        surfMD = L * np.exp(-1/(2*sig**2)*R**2)
        kappa = np.sqrt(6) * a * R
        vc = a * R**2
        sigma_R = b*R
        result = g.toomre_param(R, vc, sigma_R, surfMD)
        # convert units from M_sun to kg and pc^2 to m^2
        surfMD = surfMD * const.M_sun.value / (const.pc.value)**2
        # in km^2/(m arcsec)
        analytical_result = sigma_R * \
            (kappa)/(3.36*const.G.value*surfMD)
        # km to m and arcsec to km
        analytical_result *= 1e3 / (np.pi / (180*3600) *
                                    const.pc.value * g.distance * 1e3)
        for i in [3, 5, 10, 14, 17, 20, 24, 26, 30, 34, 35]:
            self.assertAlmostEqual(result[i], analytical_result[i])


def test_toomre():
    a = 10
    b = 200
    L = 500*4
    sigma = 50
    R = g.data_cyl['R']
    omega = a * R
    Domega = np.gradient(omega, R)
    sigma_R = 80 + -1/10 * R

    def mge(x, y):
        return L * np.exp(-1/(2*sigma**2)*x**2)

    surfMD = g.massLight * L * \
        np.exp(-1/(2*sigma**2)*R**2) * const.M_sun.value / const.pc.value**2

    kappa_sq = 4 * omega**2 + 2 * R * omega*Domega
    kappa = np.sqrt(kappa_sq) / (np.pi / (180*3600) *
                                 const.pc.value * g.distance * 1e3)

    sigma_crit = 3.36 * const.G.value * surfMD / kappa * 1e-3

    Q = sigma_R / sigma_crit

    print(np.sqrt(g.data_cyl['vRR']))

    g.data_cyl['vRR'] = sigma_R**2
    g.mge_LM = mge
    g.mge_DM = lambda x, y:  0*x
    g.data_cyl['vp'] = a * R**2

    plt.plot(R, Q, '-x')
    plt.grid()
    g.plot_total_cyl()


def test_vcric():
    R = g.data_cyl['R']
    vcric = g.velocity_vcirc()
    omega = vcric / R
    Domega = np.gradient(omega, R)
    kappa_sq = 4*omega**2 + 2 * R * omega * Domega
    plt.plot(R, kappa_sq)
    plt.grid()
    plt.show()


def test_interpolation(name):
    g = Galaxy('data/sRsz_data/', "data/massmge/", name)
    R = g.data_int['R']
    Rnorm = g.data_cyl['R']
    Vphi = np.abs(g.data_int['vp'])
    sigR = np.sqrt(g.data_int['vRR'])
    sigPhi = np.sqrt(g.data_int['vpp'])
    fig, ax = plt.subplots(1, 4, figsize=(18, 5))
    ax[0].plot(R, Vphi, label=r'$V_{\phi}$')
    ax[0].plot(Rnorm, np.sqrt(g.data_cyl['vpp']), 'o')
    ax[0].plot(R, sigR, label=r'$\sigma_R$')
    ax[0].plot(Rnorm, np.sqrt(g.data_cyl['vRR']), 'o')
    ax[0].plot(R, sigPhi, label=r'$\sigma_{\Phi}$')
    ax[0].plot(Rnorm, g.data_cyl['vp'], 'o')
    ax[0].set_xlim(0, 35)
    ax[0].legend()
    ax[0].set_title("Quadratic interpolated data")
    ax[0].set_ylabel('V [km/s]')
    ax[0].grid()
    nu = mrd.mge_radial_density(
        np.array(g.paramGaussLM[0]), np.array(g.paramGaussLM[1]), np.array(g.paramGaussLM[2]), g.inclAngle, R, g.distance)
    sum1 = np.gradient(np.log(nu*sigR**2), np.log(R))
    sum2 = (sigPhi/sigR)**2 - 1
    vc_sq = Vphi**2 + sigR**2 * (sum1 + sum2)
    ax[2].plot(R, sum1, label='Derivative Term')
    ax[2].legend()
    ax[2].grid()
    ax[3].plot(R, np.log(nu*sigR**2))
    ax[3].grid()
    ax[3].set_title('nu')
    R_vc = []
    vc = []
    # for i in range(len(vc_sq)):
    #     if vc_sq[i] > 0:
    #         vc.append(np.sqrt(vc_sq[i]))
    #         R_vc.append(R[i])
    # vc = np.array(vc)
    # R = np.array(R_vc)
    vc = np.sqrt(vc_sq)

    ax[1].plot(R, vc, label='ADC vc')
    vp = g.data_cyl['vp']
    vcirc = g.velocity_vcirc()
    ax[1].plot(Rnorm, vp, label='vp')
    ax[1].plot(Rnorm, vcirc, label='vcirc')
    Radc, vcadc = g.vc_ADC()
    ax[1].set_title('Velocities')
    ax[1].plot(Radc, vcadc, label='ADC no interpolation')
    ax[1].grid()
    ax[1].set_xlim(0, 35)
    ax[1].set_ylabel('V [km/s]')
    ax[1].legend()
    plt.suptitle(name)
    plt.savefig("figures/interpolation/"+name+"_interpol.pdf")
    plt.show()
    plt.close()


if __name__ == "__main__":
    test_interpolation('NGC6278')
    test_interpolation('NGC4210')
    test_interpolation('IC2378')
    test_interpolation('IC0674')
    test_interpolation('NGC5630')
    # unittest.main()
