import numpy as np
import matplotlib.pylab as plt
import gaussianExpansion as mge
import montecarlo as mc
from astropy.io import fits
from scipy.interpolate import interp1d
import astropy.constants as const
import mge_vcirc
import mge_radial_density
import tool_analyze as ta


class Galaxy:
    def __init__(self, name):
        """Read the data from the files: <name>_anisotropy.dat,
                                         <name>_betaz.dat,
                                         <name>_mass_mge.txt"""
        self.name = name
        self.distance = 0.0
        self.mge_Re_cyl = 0.0
        self.data = self.read_data("data/sRsz_data/" + name + '_betaz.dat')
        self.data['vp'] = np.abs(self.data['vp'])
        self.massLight = 0.0
        self.inclAngle = 0.0
        self.pGaussLM, self.pGaussDM = self.read_mass(
            "data/massmge/" + name + '_mass_mge.txt')
        self.mge_LM = mge.multi_gauss_exp(self.pGaussLM)  # in [L_sun]
        self.mge_DM = mge.multi_gauss_exp(self.pGaussDM)  # in [M_sun]

    def read_data(self, filename):
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
            paramLM (list), paramDM (list).
        """
        paramLM = []
        paramDM = []
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
                paramLM.append(
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
                paramDM.append(
                    (float(l[0]), float(l[1]), float(l[2])))
                i += 1
                n += 1
            except:
                i += 1
                pass
        return paramLM, paramDM

    def get_stellar_mass(self):
        """Calculate the stellar mass from the MGE expansion.
        Args:
            None
        Returns:
            float: stellar mass in solar masses.
        """
        stellM = []
        for i in range(len(self.pGaussLM)):
            M0 = self.massLight * self.pGaussLM[i][0]  # M_sun/pc^2
            sig = self.pGaussLM[i][1] * \
                (self.distance * np.pi / (3600*180) * 1e6)  # pc
            q = self.pGaussLM[i][2]
            M = 2 * np.pi * M0 * sig**2 * q

            stellM.append(M)
        return sum(stellM)

    def interpolate_data(self, x, y, xnew, kind='cubic'):
        """Interpolate data at the new given x values.
        Args:
            x (np.array): array of x-values
            y (np.array): array of y-values
            xnew (np.array): new x-values
        Returns:
            ynew (np.array).
        """
        f = interp1d(x, y, kind=kind)
        ynew = f(xnew)
        return ynew

    def kappa_sq(self, R, vc):
        """Calculate the epicyclic frequency kappa^2.
        Args:
            R  (np.array): radii values
            vc (np.array): velocity profile at the radii valus
        Returns:
            np.array: kappa^2 in [km^2/s^2]"""
        omega = vc / R
        Domega = np.gradient(omega, R)
        kappa_squared = 4*omega**2 + 2*R*omega*Domega
        # TODO: investigate
        # kappa_alt = 2*omega/R * np.gradient(R**2*omega, R)
        return kappa_squared

    def stellar_surfMD(self, R):
        """Calculate the stellar surface mass density at radius R.
        Args:
            R (np.array): the radial values
        Returns:
            np.array: stellar surface mass density in [M_sun/pc^2].
        """
        smd = self.massLight * self.mge_LM(R, np.zeros(len(R)))
        return smd

    def darkmatter_surfMD(self, R):
        """Calculate the stellar surface mass density at radius R.
        Args:
            R (np.array): the radial values
        Returns:
            np.array: stellar surface mass density in [M_sun/pc^2].
        """
        smd = self.mge_DM(R, np.zeros(len(R)))
        return smd

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

    def get_toomre(self):
        """Return the Toomre parameter and error.
        Args:
            None
        Returns:
            Q, eQ: np.arrays.
        """
        R = self.data['R']
        sigma_R = np.sqrt(self.data['vRR'])
        e_sigma_R = 0.5 / sigma_R * self.data['e_vRR']
        vc = self.velocity_vcirc(R)
        surfMD = self.stellar_surfMD(R)
        Q = self.toomre_param(R, vc, sigma_R, surfMD)
        eQ = Q / sigma_R * e_sigma_R
        return Q, eQ

    def swing_amplif_X(self, R):
        """Calculate the swing amplification parameter X profile.
        Using the modelled vcirc.
        Args:
            R (np.array): the radial vaules. 
        Returns:
            np.array: array of the X parameter.
        """
        vc = self.velocity_vcirc(R)
        kappa_sq = self.kappa_sq(R, vc)
        # 1/arcsec -> 1/pc
        kappa_sq = kappa_sq / (np.pi / (180*3600) * self.distance * 1e6)
        Sigma = self.massLight * self.mge_LM(R, np.zeros(len(R)))
        # pc -> m
        Sigma = Sigma * const.M_sun.value / const.pc.value
        X = R * kappa_sq * 1e6 / (4 * np.pi * const.G.value * Sigma)
        return X

    def velocity_vcirc(self, R):
        """Calculate the velocity profile from the MGE parameters.
        Args:
            R (np.array): the radial values
        Returns:
            np.array: velocity profile at the radii.
        """
        mbh = 0  # solar masses
        surfPot = [self.pGaussLM[i][0] *
                   self.massLight for i in range(len(self.pGaussLM))]
        surfPot = surfPot + [self.pGaussDM[i][0]
                             for i in range(len(self.pGaussDM))]
        surfPot = np.array(surfPot)
        sigPot = [self.pGaussLM[i][1]
                  for i in range(len(self.pGaussLM))]
        sigPot = sigPot + [self.pGaussDM[i][1]
                           for i in range(len(self.pGaussDM))]
        sigPot = np.array(sigPot)
        qPot = [self.pGaussLM[i][2] for i in range(len(self.pGaussLM))]
        qPot = qPot + [self.pGaussDM[i][2]
                       for i in range(len(self.pGaussDM))]
        qPot = np.array(qPot)
        vcirc = mge_vcirc.mge_vcirc(
            surfPot, sigPot, qPot, self.inclAngle, mbh, self.distance, R)
        return vcirc

    def write_toomre(self):
        """Write the Toomre parameter to file data/toomre/<name>_Q_data.txt.
        Args:
            None
        Returns:
            None.
        """
        R = self.data['R']
        sigma_R = np.sqrt(self.data['vRR'])
        vc = self.velocity_vcirc(R)
        surfMD = self.stellar_surfMD(R)
        Q = self.toomre_param(R, vc, sigma_R, surfMD)
        f = open("data/toomre/"+self.name+"_Q_data.txt", 'w')
        f.write('# Name: {}\n'.format(self.name))
        f.write(
            '# Radius[arcsec]\t\tQ\t\tsigmaR[km/s]\t\tvc[km/s]\t\tsurfMD[M_sun/pc^2]\n')
        for i in range(len(R)):
            f.write(str(R[i]) + '\t\t')
            f.write(str(Q[i]) + '\t\t')
            f.write(str(sigma_R[i])+'\t\t')
            f.write(str(vc[i])+'\t\t')
            f.write(str(surfMD[i])+'\t\t')
            f.write('\n')

    def write_toomre_interpl(self):
        """Write the toomre parameter Q to the file 'data/toomre/<name>_Qintpol_data.txt.
        Args:
            None
        Returns:
            None.
        """
        R = np.linspace(np.min(self.data['R']), np.max(self.data['R']), 300)
        vRR = self.interpolate_data(self.data['R'], self.data['vRR'], R)
        sigma_z = np.sqrt(self.interpolate_data(
            self.data['R'], self.data['vzz'], R))
        sigma_R = np.sqrt(vRR)
        vc = self.velocity_vcirc(R)
        surfMD = self.stellar_surfMD(R)
        Q = self.toomre_param(R, vc, sigma_R, surfMD)
        f = open("data/toomre/"+self.name+"_Qintpol_data.txt", 'w')
        f.write('# Name: {}\n'.format(self.name))
        f.write(
            '# Radius[arcsec]\t\tQ\t\tsigmaR[km/s]\t\tvc[km/s]\t\tsurfMD[M_sun/pc^2]\t\t')
        f.write('sigma_z[km/s]\n')
        for i in range(len(R)):
            f.write(str(R[i]) + '\t\t')
            f.write(str(Q[i]) + '\t\t')
            f.write(str(sigma_R[i])+'\t\t')
            f.write(str(vc[i])+'\t\t')
            f.write(str(surfMD[i])+'\t\t')
            f.write(str(sigma_z[i])+'\t\t')
            f.write('\n')

    def write_swing_ampli(self):
        """Write swing amplification parameter X to file 'data/swing/<name>_X_data.txt.
        Args:
            None
        Returns:
            None.
        """
        R = self.data['R']
        X = self.swing_amplif_X(R)
        f = open("data/swing/"+self.name+"_X_data.txt", 'w')
        f.write('# Name: {}\n'.format(self.name))
        f.write('# Radius[arcsec]\t\tX\n')
        for i in range(len(R)):
            f.write(str(R[i])+'\t\t')
            f.write(str(X[i]) + '\t\t')
            f.write('\n')

    def write_swing_ampli_interpl(self):
        """Write swing amplification parameter X to file 'data/swing/<name>_X_data.txt.
        Args:
            None
        Returns:
            None.
        """
        R = np.linspace(np.min(self.data['R']), np.max(self.data['R']), 300)
        X = self.swing_amplif_X(R)
        f = open("data/swing/"+self.name+"_Xintpol_data.txt", 'w')
        f.write('# Name: {}\n'.format(self.name))
        f.write('# Radius[arcsec]\t\tX\n')
        for i in range(len(R)):
            f.write(str(R[i])+'\t\t')
            f.write(str(X[i]) + '\t\t')
            f.write('\n')


if __name__ == "__main__":
    names = ta.read_filenames("data/massmge/", "_mass_mge.txt")
    for n in names:
        try:
            g = Galaxy(n)
            # g.write_toomre()
            g.write_toomre_interpl()
            # g.write_swing_ampli()
            # g.write_swing_ampli_interpl()
        except Exception as e:
            print(e)
            print(n)
