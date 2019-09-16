import numpy as np
import matplotlib.pylab as plt
import tqdm
import gaussianExpansion as mge
import montecarlo as mc
from astropy.io import fits
import astropy.constants as const
from scipy.interpolate import interp1d
import mge_vcirc
import mge_radial_density


class Galaxy:
    def __init__(self, betaz, massmge, name):
        """Read the data from the files: <name>_anisotropy.dat,
                                         <name>_betaz.dat,
                                         <name>_mass_mge.txt"""
        self.name = name
        self.distance = 0.0
        self.mge_Re_sph = 0.0
        self.mge_Re_cyl = 0.0
        self.data_cyl = self.read_cylindrical(betaz + name + '_betaz.dat')
        self.data_cyl['vp'] = np.abs(self.data_cyl['vp'])
        self.massLight = 0.0
        self.inclAngle = 0.0
        self.paramGaussLM = list()
        self.paramGaussDM = list()
        self.read_mass(massmge + name + '_mass_mge.txt')
        self.mge_LM = mge.multi_gauss_exp(self.paramGaussLM)  # in [L_sun]
        self.mge_DM = mge.multi_gauss_exp(self.paramGaussDM)  # in [M_sun]
        self.mcSigmaR = self.mc_sigma_R(N=1000)
        self.mcKappaSq = 0
        self.mass = self.get_total_mass()
        self.bKine = True
        try:
            self.kine = fits.open(
                "data/kinemetry/"+name+".V1200.rscube_INDOUSv2_SN20_" +
                "stellar_kin_elliptical_kinemetry.fits")[1].data
        except:
            self.bKine = False
        # interpolate data
        self.data_int = self.get_interpolated_data()

    def read_cylindrical(self, filename):
        """Read the columns of the data-file in the format:
        Radius, betaZ, e_betaZ, vzz, error, vRR, error, vpp, error, vp, error
        Args:
            filename (str): name of the file
        Returns:
            tuple of lists
        """
        data = {'R': [], 'betaZ': [], 'e_betaZ': [], 'vzz': [], 'e_vzz': [],
                'vRR': [], 'e_vRR': [],  'vpp': [], 'e_vpp': [], 'vp': [], 'e_vp': []}
        f = open(filename, 'r')
        lines = f.readlines()
        line = lines[0].split()
        self.distance = float(line[0])
        line = lines[0].split()
        self.mge_Re_cyl = float(line[0])
        for line in lines:
            line = line.split()
            if len(line) == 11:
                try:
                    data['R'].append(float(line[0]))
                    data['betaZ'].append(float(line[1]))
                    data['e_betaZ'].append(float(line[2]))
                    data['vzz'].append(float(line[3]))
                    data['e_vzz'].append(float(line[4]))
                    data['vRR'].append(float(line[5]))
                    data['e_vRR'].append(float(line[6]))
                    data['vpp'].append(float(line[7]))
                    data['e_vpp'].append(float(line[8]))
                    data['vp'].append(float(line[9]))
                    data['e_vp'].append(float(line[10]))
                except:
                    print("Error: Converting string into float:")
                    print(line)

        # exclude data points with values 0 at n > 0
        n = 0
        for i in range(len(data['R'])):
            # if i > 0 and data['R'][i] == 0:
            #     n = i
            #     break
            if i > 0 and np.isnan(data['e_betaZ'][i]):
                n = i
                break
        skipFirst = np.isnan(data['e_betaZ'][0])
        for i in data:
            data[i] = np.array(data[i][:n])
            if skipFirst:
                data[i] = np.array(data[i][1:])
        return data

    def read_mass(self, filename):
        """Read the galaxy_mass_mge.txt file.
        Args:
            filename (str): name of the file
        Returns:
            None."""
        f = open(filename, 'r')
        lines = f.readlines()
        line = lines[0].split()
        self.massLight = float(line[0])
        line = lines[1].split()
        self.inclAngle = float(line[0])
        line = lines[2].split()
        nGaussLM = int(line[0])
        line = lines[3].split()
        nGaussDM = int(line[0])
        i = 4
        n = 0
        while n < nGaussLM:
            l = lines[i].split()
            try:
                self.paramGaussLM.append(
                    (float(l[0]), float(l[1]), float(l[2])))
                i += 1
                n += 1
            except:
                i += 1
                pass
        n = 0
        while n < nGaussDM:
            l = lines[i].split()
            try:
                self.paramGaussDM.append(
                    (float(l[0]), float(l[1]), float(l[2])))
                i += 1
                n += 1
            except:
                i += 1
                pass

    def interpolate_quad(self, x, y, xnew):
        """Quadratic interpolation of the given data.
        Returns the interpolated y-values at the xnew x-values.
        Args:
            x (np.array): array of the x-values
            y (np.array): array of the y-values
        Returns:
            np.array: ynew - the interpolated data.
        """
        f = interp1d(x, y, kind='quadratic')
        ynew = f(xnew)
        return ynew

    def get_interpolated_data(self):
        """Interpolate vRR, vpp, vzz and vp.
        Args:
            None
        Returns:
            dict: the interpolated data and the radius.
        """
        data = {}
        rmax = 35
        if rmax > self.data_cyl['R'][-1]:
            rmax = self.data_cyl['R'][-1]
        Rint = np.linspace(self.data_cyl['R'][0], rmax, 5000)
        R = self.data_cyl['R']
        data['R'] = Rint
        data['vpp'] = self.interpolate_quad(R, self.data_cyl['vpp'], Rint)
        data['vRR'] = self.interpolate_quad(R, self.data_cyl['vRR'], Rint)
        data['vzz'] = self.interpolate_quad(R, self.data_cyl['vzz'], Rint)
        data['vp'] = self.interpolate_quad(R, self.data_cyl['vp'], Rint)
        return data

    def kappa_sq(self, R, vc):
        """Calculate the epicyclic frequency kappa^2.
        Args:
            R  (np.array): radii values
            vc (np.array): velocity profile at the radii valus
        Returns:
            np.array: kappa^2"""
        omega = vc / R
        Domega = np.gradient(omega, R)
        kappa_squared = 4*omega**2 + 2*R*omega*Domega
        # TODO: investigate
        # kappa_alt = 2*omega/R * np.gradient(R**2*omega, R)
        return kappa_squared

    def mc_kappa_sq_curves(self, R, vc, e_vc, N=1000):
        """Simulate variations of the kappa^2 profile based on normal distributed
        errors.
        Args:
            R (np.array): radius values
            vc (np.array): velocity profile at radius values
            e_vc (np.array): error of velocity profile
            N (int): number of randomly chosen points per data-point
        Returns:
            np.array: array of arrays with different curves."""
        l_omega = []
        for i in range(len(R)):
            mc_vc = np.random.normal(vc[i], e_vc[i], N)
            omega = mc_vc / R[i]
            l_omega.append(omega)
        omega = np.transpose(np.array(l_omega))
        trash, Domega = np.gradient(omega, 1, R)
        kappa_sq = 4 * omega**2 + 2 * R * omega*Domega
        return kappa_sq

    def mc_sigma_R(self, N, coord='cyl'):
        """Generate N curves of the velocity dispersion sigma_R with the data points
        normal distributed within the error.
        Args:
            N (int): number of curves
        Returns:
            np.array: 2D array with the curves."""
        if coord == "sph":
            R = self.data_sph['R']
            vrr = self.data_sph['vrr']
            e_vrr = self.data_sph['e_vrr']
        elif coord == "cyl":
            R = self.data_cyl['R']
            vrr = self.data_cyl['vRR']
            e_vrr = self.data_cyl['e_vRR']
        else:
            raise Exception("coord must be a string with value 'sph' or 'cyl'.\n\
            The value was: {}, type {}".format(coord, type(coord)))
        mc_vrr = mc.simulate_errors_fast(R, vrr, e_vrr, N=N)
        while True:
            mc_vrr = mc_vrr[mc_vrr.min(axis=1) >= 0, :]
            if len(mc_vrr) >= N:
                mc_vrr = mc_vrr[:N]
                break
            else:
                ar = mc.simulate_errors_fast(
                    R, vrr, e_vrr, N=2*(N-len(mc_vrr)))
                mc_vrr = np.append(mc_vrr, ar, axis=0)
        mc_sig_R = np.sqrt(mc_vrr)
        return mc_sig_R

    def mc_toomre(self, N=1000):
        """Generate N curves of the toomre parameter Q by the variation of the
        velocity dispersion sigma_R and the epicyclic frequency kappa^2.
        Note that curves with negative kappa^2 are excluded.
        Args:
            N (int): number of curves.
        Returns:
            np.array: the toomre parameter at radii"""
        R = self.data_cyl['R']
        vc = self.data_cyl['vp']
        e_vc = self.data_cyl['e_vp']
        kappa_sq = self.kappa_sq(R, vc)
        e_kappa_sq = self.propagate_e_kappa_sq(R, vc/R, e_vc)

        l_kappa_sq = []
        # make sure no negative values for kappa^2 appear
        with tqdm.tqdm(total=N, desc='Kappa Squared') as pbar:
            barProgress = 0
            while len(l_kappa_sq) < N:
                barProgress = len(l_kappa_sq)
                # pbar.write('%i' % len(l_kappa_sq))
                k_sq = self.mc_kappa_sq_curves(R, vc, e_vc, N)
                for i in range(len(kappa_sq)):
                    # exclude values with 3-sigma below 0
                    # set them to 0
                    if kappa_sq[i] + 3*e_kappa_sq[i] < 0:
                        # pbar.write("Skipped Kappa^2 value: %f" %
                        #            kappa_sq[i], end='')
                        # pbar.write(" error: %f" % e_kappa_sq[i])
                        k_sq[:, i] = 0.0
                    if len(k_sq[i][k_sq[i] >= 0]) == len(R):
                        l_kappa_sq.append(k_sq[i])
                pbar.update(len(l_kappa_sq)-barProgress)
        kappa_sq = np.array(l_kappa_sq)[:N]
        kappa = np.sqrt(kappa_sq)
        mcSigR = self.mc_sigma_R(N, coord='cyl')
        surfMD = self.massLight * \
            self.mge_LM(R, np.zeros(len(R))) + self.mge_DM(R, np.zeros(len(R)))
        # unit conversions
        # Sigma in kg / m^2
        surfMD = surfMD * const.M_sun.value / const.pc.value**2
        # arcsec -> km
        kappa = kappa / (self.distance * const.kpc.value * np.pi / (3600*180))
        Q = mcSigR * kappa * 1e3 / (3.36 * const.G.value * surfMD)
        return Q

    def plot_luminosity_LM(self):
        x = np.linspace(0, 50, 500)
        y = np.linspace(0, 50, 500)
        plt.semilogy(x, self.mge_LM(x, np.zeros(len(x))),
                     label='major axis', color='red')
        plt.semilogy(y, self.mge_LM(np.zeros(len(y)), y),
                     '--', label='minor axis', color='red')
        plt.ylim(1, 1e4)
        plt.xlim(0, 50)
        plt.legend()
        plt.grid()
        plt.title('Surface Brightness ' + self.name)
        plt.xlabel('R [arcsec]')
        plt.ylabel(r'Flux [L$_{sun}$/pc$^2$]')
        plt.show()

    def toomre_param(self, R, vc, sigma_R, surfMD):
        """Calculate the Toomre paramter from the data.
        Args:
            R (np.array): radius values in [arcsec]
            vc (np.array): velocity profile in [km/s]
            sigma_R (np.array): velocity dispersion in R direction in [km/s]
            surfMD (np.array): surface mass density in [M_sun/pc^2]
        Returns:
            np.array: array of Toomre parameter"""
        kappa = np.sqrt(self.kappa_sq(R, vc))
        # unit conversions
        # Sigma in kg / m^2
        surfMD = surfMD * const.M_sun.value / (const.pc.value)**2
        # arcsec -> km
        kappa = kappa / (np.pi / (180*3600) *
                         const.pc.value * self.distance * 1e3)
        sigma_crit = 3.36 * const.G.value * surfMD / kappa * 1e-3
        Q = sigma_R / sigma_crit
        return Q

    def plot_mass_profiles(self):
        def mass_LM(x, y): return self.massLight * self.mge_LM(x, y)
        R = np.linspace(0, 50, 500)
        y = np.zeros(len(R))
        plt.semilogy(R, mass_LM(R, y), '--',
                     color='blue', label='Stellar Mass')
        plt.semilogy(R, self.mge_DM(R, y), '--',
                     color='black', label='Dark Matter')
        plt.semilogy(R, self.mge_DM(R, y)+mass_LM(R, y),
                     color='red', label='Sum')
        plt.title(r'Surface Mass Density $\Sigma$ '+self.name)
        plt.xlabel('R [arcsec]')
        plt.ylabel(r'$\Sigma$ [M$_{sun}/$pc$^2$]')
        plt.grid()
        plt.xlim(0, 50)
        plt.legend()
        plt.show()

    def plot_total_cyl(self, show=False, save=True):
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        R = self.data_cyl['R']

        # velocity profile
        vp, e_vp = self.data_cyl['vp'], self.data_cyl['e_vp']
        R_ADC, vcADC = self.vc_ADC()
        axes[0][1].plot(R_ADC, vcADC, label=r'$V_c$ with ADC')
        bVcirc = False
        try:
            vcirc = self.velocity_vcirc()
            axes[0][1].plot(R, vcirc, label=r'vcirc', color='green')
            bVcirc = True
        except:
            pass
        if self.bKine:
            axes[0][1].plot(self.kine["RADIUS"][0],
                            self.kine["VPROF"][0], label='Kinemetry VP')
        axes[0][1].plot(R, np.abs(vp), '-x', label=r'V$_c$', color='blue')
        axes[0][1].errorbar(R, np.abs(vp), yerr=e_vp, fmt='none',
                            capsize=2, color='blue')
        axes[0][1].grid()
        axes[0][1].legend()
        axes[0][1].set_xlabel(r'R [arcsec]')
        axes[0][1].set_ylabel(r'|V$_c$| [km/s]')
        axes[0][1].set_title(r'Velocity Profile V$_c$')

        kappa_sq = self.kappa_sq(R, vp)
        e_kappa_sq = self.propagate_e_kappa_sq(R, vp/R, e_vp)
        mc_kappa_sq = self.mc_kappa_sq_curves(R, vp, e_vp, N=1000)
        mc_kappa_sq, e_mc_kappa_sq = mc.mean_2d_array(mc_kappa_sq)

        axes[0][0].plot(R, kappa_sq, '-x', label=r'$\kappa^2$')
        axes[0][0].errorbar(R, kappa_sq, yerr=e_kappa_sq,
                            fmt='none', capsize=2)
        axes[0][0].plot(R, mc_kappa_sq, label=r'MC $\kappa^2$', color='blue')
        axes[0][0].errorbar(R, mc_kappa_sq, yerr=e_mc_kappa_sq,
                            fmt='none', capsize=2, color='blue')
        if bVcirc:
            kappa_vcirc_sq = self.kappa_vcirc()**2
            axes[0][0].plot(R, kappa_vcirc_sq,
                            label=r'$\kappa^2$ with vcirc', color='green')
        axes[0][0].legend()
        axes[0][0].grid()
        axes[0][0].set_xlabel(r'R [arcsec]')
        axes[0][0].set_ylabel(r'$\kappa^2$ [(km s$^{-1}$arcsec$^{-1}$)$^2$]')
        axes[0][0].set_title(r'Epicyclic Frequency $\kappa^2$')

        # surface mass density
        surfMD = self.massLight * \
            self.mge_LM(R, np.zeros(len(R))) + self.mge_DM(R, np.zeros(len(R)))
        axes[1][0].semilogy(R, surfMD, color='red', label=r'$\Sigma$')
        axes[1][0].semilogy(R, self.mge_DM(R, np.zeros(len(R))),
                            '--', color='black', label=r'$\Sigma$ Dark Matter')
        axes[1][0].semilogy(R, self.massLight * self.mge_LM(R, np.zeros(len(R))),
                            '--', color='blue', label=r'$\Sigma$ Stellar Mass')
        axes[1][0].legend()
        axes[1][0].grid()
        axes[1][0].set_xlabel('R [arcsec]')
        axes[1][0].set_ylabel(r'$\Sigma$ [M$_{Sun}$ pc$^{-2}$]')
        axes[1][0].set_title(r'Surface Mass Density $\Sigma$')

        # toomre parameter
        sigma_R = np.sqrt(self.data_cyl['vRR'])
        # Q = self.toomre_param(R, vp, sigma_R, surfMD)
        mcQ = self.mc_toomre(1000)
        m_mcQ, e_m_mcQ = mc.mean_2d_array(mcQ)
        # axes[1][1].plot(R, Q, 'x', label=r'$Q_{star}$')
        axes[1][1].plot(R, m_mcQ, label=r'MC Q')
        axes[1][1].errorbar(R, m_mcQ, yerr=e_m_mcQ, fmt='none', capsize=2)
        if bVcirc:
            Q_vcirc = self.toomre_param(R, vcirc, sigma_R, surfMD)
            axes[1][1].plot(R, Q_vcirc, label=r'Q vcirc', color='green')
        axes[1][1].set_title(r'Toomre Parameter Q')
        axes[1][1].set_xlabel(r'R [arcsec]')
        axes[1][1].grid()
        axes[1][1].legend()
        if save:
            plt.savefig('figures/toomre/'+self.name+'_toomre.pdf')
            self.write_plot(mc_kappa_sq, e_mc_kappa_sq, m_mcQ, e_m_mcQ)
        if show:
            plt.show()
        plt.close()

    def propagate_e_kappa_sq(self, R, omega, e_vc):
        Domega = np.gradient(omega, R)
        return np.abs(8 * omega + 2 * R * Domega) / R * e_vc

    def get_toomre_vcirc(self):
        """Calculate Q from the vcirc modelling and the error.
        Args:
            None
        Returns:
            np.array: array of Q.
        """
        R = self.data_cyl['R']
        vcirc = self.velocity_vcirc()
        sigma_R = np.sqrt(self.data_cyl['vRR'])
        e_sigma_R = 0.5 / sigma_R * self.data_cyl['e_vRR']
        surfMD = self.massLight * \
            self.mge_LM(R, np.zeros(len(R))) + self.mge_DM(R, np.zeros(len(R)))
        Q_vcirc = self.toomre_param(R, vcirc, sigma_R, surfMD)
        e_Q_vcirc = Q_vcirc / sigma_R * e_sigma_R
        return Q_vcirc, e_Q_vcirc

    def velocity_vcirc(self):
        """Calculate the velocity profile from the MGE parameters.
        Args:
            None
        Returns:
            np.array: velocity profile at the radii"""
        R = self.data_cyl['R']
        mbh = 0  # solar masses
        surfPot = [self.paramGaussLM[i][0] *
                   self.massLight for i in range(len(self.paramGaussLM))]
        surfPot = surfPot + [self.paramGaussDM[i][0]
                             for i in range(len(self.paramGaussDM))]
        surfPot = np.array(surfPot)
        sigPot = [self.paramGaussLM[i][1]
                  for i in range(len(self.paramGaussLM))]
        sigPot = sigPot + [self.paramGaussDM[i][1]
                           for i in range(len(self.paramGaussDM))]
        sigPot = np.array(sigPot)
        qPot = [self.paramGaussLM[i][2] for i in range(len(self.paramGaussLM))]
        qPot = qPot + [self.paramGaussDM[i][2]
                       for i in range(len(self.paramGaussDM))]
        qPot = np.array(qPot)
        vcirc = mge_vcirc.mge_vcirc(
            surfPot, sigPot, qPot, self.inclAngle, mbh, self.distance, R)
        return vcirc

    def kappa_vcirc(self):
        """Calculate the epicyclic frequency from the vcirc modelled velocity profile.
        Args:
            None:
        Returns:
            np.array: array of kappa
        """
        R = np.linspace(self.data_cyl['R'][0]/2,
                        self.data_cyl['R'][-1]+1, 1000)
        R = np.append(R, self.data_cyl['R'])
        R = np.sort(R)
        mbh = 0  # solar masses
        surfPot = [self.paramGaussLM[i][0] *
                   self.massLight for i in range(len(self.paramGaussLM))]
        surfPot = surfPot + [self.paramGaussDM[i][0]
                             for i in range(len(self.paramGaussDM))]
        surfPot = np.array(surfPot)
        sigPot = [self.paramGaussLM[i][1]
                  for i in range(len(self.paramGaussLM))]
        sigPot = sigPot + [self.paramGaussDM[i][1]
                           for i in range(len(self.paramGaussDM))]
        sigPot = np.array(sigPot)
        qPot = [self.paramGaussLM[i][2] for i in range(len(self.paramGaussLM))]
        qPot = qPot + [self.paramGaussDM[i][2]
                       for i in range(len(self.paramGaussDM))]
        qPot = np.array(qPot)
        vcirc = mge_vcirc.mge_vcirc(
            surfPot, sigPot, qPot, self.inclAngle, mbh, self.distance, R)

        omega = vcirc/R
        Domega = np.gradient(omega, R)
        kappa = np.sqrt(4*omega**2 + 2*R*omega*Domega)
        l_kappa = []

        for i in range(len(self.data_cyl['R'])):
            ind = np.where(R == self.data_cyl['R'][i])
            l_kappa.append(kappa[ind][0])
        return np.array(l_kappa)

    def get_stellar_mass(self):
        """Calculate the stellar mass from the MGE expansion.
        Args:
            None
        Returns:
            float: stellar mass in solar masses.
        """
        stellM = []
        for i in range(len(self.paramGaussLM)):
            M0 = self.massLight * self.paramGaussLM[i][0]  # M_sun/pc^2
            sig = self.paramGaussLM[i][1] * \
                (self.distance * np.pi / (3600*180) * 1e6)  # pc
            q = self.paramGaussLM[i][2]
            M = 2 * np.pi * M0 * sig**2 * q
            stellM.append(M)
        return sum(stellM)

    def get_total_mass(self):
        """Calculate the total mass from the MGE expansion.
        Args:
            None.
        Returns:
            float: total mass in solar masses
        """
        totM = []
        for i in range(len(self.paramGaussLM)):
            M0 = self.massLight * self.paramGaussLM[i][0]  # M_sun/pc^2
            sig = self.paramGaussLM[i][1] * \
                (self.distance * np.pi / (3600*180) * 1e6)  # pc
            q = self.paramGaussLM[i][2]
            M = 2 * np.pi * M0 * sig**2 * q
            totM.append(M)
        for i in range(len(self.paramGaussDM)):
            M0 = self.paramGaussDM[i][0]
            sig = self.paramGaussDM[i][1] * \
                (self.distance * np.pi / (3600*180) * 1e6)
            q = self.paramGaussDM[i][2]
            M = 2 * np.pi * M0 * sig**2 * q
            totM.append(M)
        return sum(totM)

    def get_surfaceMD(self, matter='all'):
        """Return the surface mass density at the radial values
        in [M_sun/pc^2].
        Args:
            None
        Returns:
            np.array: the values.
        """
        R = self.data_cyl['R']
        y = np.zeros(len(R))
        if matter == 'all':
            smd = self.massLight * self.mge_LM(R, y) + self.mge_DM(R, y)
        elif matter == 'lm':
            smd = self.massLight * self.mge_LM(R, y)
        elif matter == 'dm':
            smd = self.mge_DM(R, y)
        return smd

    def test_total_mass(self):
        from scipy.integrate import nquad
        from scipi.special import erf
        m = []
        ml = []
        md = []
        R = np.linspace(0, 40, 10)
        for i in range(len(R)):
            MassLm = nquad(self.mge_LM, [[-R[i], R[i]],  [-R[i], R[i]]])
            MassLm = MassLm[0] * self.massLight * \
                (self.distance * np.pi / (3600*180) * 1e6)**2
            MassDm = nquad(self.mge_DM, [[-R[i], R[i]],  [-R[i], R[i]]])[0]
            MassDm *= (self.distance * np.pi / (3600*180) * 1e6)**2
            ml.append(MassLm)
            md.append(MassDm)
            print("cululative mass at radius {:} arcs: {:2.5g}".format(
                R[i], MassLm))
            m.append(MassLm + MassDm)
        tot_ml = nquad(self.mge_DM, [[-np.inf, np.inf], [-np.inf, np.inf]])
        tot_ml = tot_ml[0] * self.massLight * \
            (self.distance * np.pi / (3600*180)*1e6)**2
        tot_md = nquad(self.mge_DM, [[-np.inf, np.inf], [-np.inf, np.inf]])
        tot_md = tot_md[0] * (self.distance * np.pi / (3600*180) * 1e6)**2
        print('Total integrated Mass: {:2.5g}'.format(tot_md+tot_ml))
        print("Total calculated Mass: {:2.5g}".format(self.get_total_mass()))
        plt.plot(R, m, label='Total')
        plt.plot(R, ml, label='Light Matter')
        plt.plot(R, md, label='Dark Matter')
        plt.grid()
        plt.legend()
        plt.show()
        plt.close()

    def write_plot(self, kappa_sq, e_kappa_sq, Q, e_Q):
        """Plot Toomre Parameter and write data to file.
        Args:
            kappa_sq (np.array): array of the kappa_sq values.
        Returns:
            None."""
        R = self.data_cyl['R']  # in arcsec
        R_pc = self.data_cyl['R'] * \
            (np.pi / (180*3600) * self.distance * 1e3)  # in kpc
        sigma_R = np.sqrt(self.data_cyl['vRR'])
        e_sigma_R = 0.5 / sigma_R * self.data_cyl['e_vRR']
        vc = self.data_cyl['vp']
        e_vc = self.data_cyl['e_vp']
        bVcirc = False
        try:
            vcirc = self.velocity_vcirc()
            kappa_vcirc = self.kappa_vcirc()
            bVcirc = True
        except:
            pass

        # Q = self.
        Sigma = self.massLight * \
            self.mge_LM(R, np.zeros(len(R))) + self.mge_DM(R, np.zeros(len(R)))

        f = open('data/toomre/'+self.name+'_data.txt', 'w')
        f.write(self.name+'\t#Name\n')
        f.write(str(self.distance)+'\t#Distance [Mpc]\n')
        f.write(str(self.mass)+'\t#Total Mass [M_sun]\n')
        f.write(str(self.inclAngle)+'\t#Inclination Angle [deg]\n')
        f.write(
            '# Radius[arcs]\tRadius[kpc]\tQ\te_Q\tKappa^2[(km/s/arcs)^2]\t')
        f.write('e_Kappa^2\tSigma[M_sun/pc ^ 2]\tsigma_R[km/s]\t')
        f.write('e_sigma_R\tvc[km/s]\te_vc\tvcirc\tkappa_vcirc\t\n')
        for i in range(len(R)):
            f.write(str(R[i])+'\t')
            f.write(str(R_pc[i])+'\t')
            f.write(str(Q[i])+'\t')
            f.write(str(e_Q[i])+'\t')
            f.write(str(kappa_sq[i])+'\t')
            f.write(str(e_kappa_sq[i])+'\t')
            f.write(str(Sigma[i])+'\t')
            f.write(str(sigma_R[i])+'\t')
            f.write(str(e_sigma_R[i])+'\t')
            f.write(str(vc[i])+'\t')
            f.write(str(e_vc[i])+'\t')
            if bVcirc:
                f.write(str(vcirc[i])+'\t')
                f.write(str(kappa_vcirc[i])+'\t')
            else:
                f.write('nan'+'\t')
                f.write('nan'+'\t')
            f.write('\n')

    def write_swing_ampli(self):
        """Write swing amplification parameter X to file 'data/swing/<name>_X_data.txt.
        Args:
            None
        Returns:
            None.
        """
        R = self.data_cyl['R']
        X = self.swing_amplif_X()
        f = open("data/swing/"+self.name+"_X_data.txt", 'w')
        f.write(self.name+'\t#Name\n')
        f.write('# Radius[arcsec]\t\tX\n')
        for i in range(len(R)):
            f.write(str(R[i])+'\t\t')
            f.write(str(X[i]) + '\t\t')
            f.write('\n')

    def vc_ADC(self):
        """Calculate the velocity with the asymetric drift correction ADC.
        The radial density is computed by the sum of stellar mass and dark matter.
        Exclude all negative values which are created by the derivative.
        Args:
             None
        Returns:
             R, Vc: Radius and the corrected circular velocity.
        """
        R = self.data_cyl["R"]
        Vphi = np.abs(self.data_cyl["vp"])
        sig_R = np.sqrt(self.data_cyl["vRR"])
        sig_Phi = np.sqrt(self.data_cyl["vpp"])
        try:
            nu = mge_radial_density.mge_radial_density(
                np.array(self.paramGaussLM[0]), np.array(
                    self.paramGaussLM[1]), np.array(self.paramGaussLM[2]),
                self.inclAngle, R, self.distance)
            # 1/pc^3 -> 1/km^3
            # nu = nu / (const.pc.value * 1e-3)**3

        except:
            print("inclination Angle:", self.inclAngle, self.name)
        sum1 = np.gradient(np.log(nu*sig_R**2), np.log(R))
        sum2 = (sig_Phi/sig_R)**2 - 1
        vc_sq = Vphi**2 + sig_R**2 * (sum1 + sum2)
        vc = []
        R_vc = []

        # fig, ax = plt.subplots(1, 4, figsize=(15, 5))
        # ax[0].plot(R, nu)
        # ax[0].set_title("Intrinsic luminosity density")
        # ax[0].grid()
        # ax[1].plot(R, sum1)
        # ax[1].grid()
        # ax[1].set_title("Derivative")
        # ax[2].plot(R, sig_R)
        # ax[2].set_title("sig_R")
        # ax[2].grid()
        # plt.show()
        # plt.close()

        for i in range(len(vc_sq)):
            if vc_sq[i] > 0:
                vc.append(np.sqrt(vc_sq[i]))
                R_vc.append(R[i])
        vc = np.array(vc)
        R_vc = np.array(R_vc)
        return R_vc, vc

    def read_data(self, key):
        """Read data written to '_data.txt' files and return selected columns.
        Args:
            key (str or list): the columns to read.
        Returns:
            np.array: a 1D or 2D array of the column values
        """
        f = open("data/toomre/"+self.name+"_data.txt", 'r')
        lines = f.readlines()
        if type(key) == list:
            data = [[] for i in len(key)]
        elif type(key) == str:
            data = []
        for l in lines:
            l = l.split()
            if len(l) >= 6 and l[0] != '#':
                if type(key) == str:
                    if key == 'kappa_sq':
                        data.append(float(l[4]))
                else:
                    print("Not implemented!")
        return np.array(data)

    def swing_amplif_X(self):
        """Calculate the swing amplification parameter X profile.
        Using the modelled vcirc.
        Args:
            None:
        Returns:
            np.array: array of the X parameter.
        """
        R = self.data_cyl['R']
        kappa_sq = self.kappa_vcirc()**2
        # 1/arcsec -> 1/pc
        kappa_sq = kappa_sq / (np.pi / (180*3600) * self.distance * 1e6)
        Sigma = self.massLight * \
            self.mge_LM(R, np.zeros(len(R))) + self.mge_DM(R, np.zeros(len(R)))
        # pc -> m
        Sigma = Sigma * const.M_sun.value / const.pc.value
        X = R * kappa_sq * 1e6 / (4 * np.pi * const.G.value * Sigma)
        return X

    def plot_Vc(self):
        """Doc modelling"""
        R = self.data_cyl['R']
        vp, e_vp = self.data_cyl['vp'], self.data_cyl['e_vp']
        R_ADC, vcADC = self.vc_ADC()
        try:
            vcirc = self.velocity_vcirc()
            plt.plot(R, vcirc, label=r'vcirc MGE modelled')
        except:
            pass
        if self.bKine:
            plt.plot(self.kine['RADIUS'][0], self.kine['VPROF']
                     [0], label=r'Kinemetry $v_p$')
        plt.plot(R, vp, label=r'v_p', color='blue')
        plt.errorbar(R, vp, yerr=e_vp, fmt='none',
                     capsize=2, color='blue')
        plt.plot(R_ADC, vcADC, label=r'$v_c$ with ADC')
        plt.grid()
        plt.legend()
        plt.xlabel('R [arcsec]')
        plt.ylabel('Velocities [km/s]')
        plt.title(r'Velocity profiles ' + self.name)
        plt.savefig("figures/vcProfiles/"+self.name+'_vc.pdf')
        plt.show()
        plt.close()


if __name__ == "__main__":
    g = Galaxy('data/sRsz_data/', 'data/massmge/', 'NGC6278')
    # g.vc_ADC()
    g.swing_amplif_X()
    # g.velocity_vcirc()
    # g.kappa_sq_vcirc()
    # g.plot_total_cyl(show=True, save=False)
    # print("Total mass NGC: {:2.5g}".format(g.get_total_mass()))
