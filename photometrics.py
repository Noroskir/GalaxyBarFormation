import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import fsolve


class Photometrics:
    def __init__(self, filename):
        """Read fits file.
        Args:
            filename (str): name of the file
        Returns:
            HDUList."""
        self.hdul = fits.open(filename)
        self.data = self.hdul[1].data
        l = len(self.data['GALAXY'][0])
        self.index = {self.data['GALAXY'][0][i]: i for i in range(l)}

    def get_param(self, name, param, band='R'):
        """Get a parameter for a specific galaxy.
        Args:
            name (str): the name of the galaxy
            param (str): the parameter
        Returns:
            float: value of the parameter.
        """
        p = self.data[param+'_'+band][0][self.index[name]]
        return p

    def is_barred(self, names):
        """Check if galaxies are barred or not.
        Args:
            names (list or str): name(s) of galaxies
        Returns:
            list: list of True and False
            """
        if type(names) == list:
            barred = []
            for n in names:
                barred.append(self.data['BAR'][0][self.index[n]])
            return barred
        elif type(names) == str:
            return self.data['BAR'][0][self.index[names]]

    def is_disked(self, names):
        """Check if galaxies have a disk or not.
        Args:
            names (list or str): name(s) of galaxies
        Returns:
            list: list of True and False
            """
        if type(names) == list:
            barred = []
            for n in names:
                barred.append(self.data['DISK'][0][self.index[n]])
            return barred
        elif type(names) == str:
            return self.data['DISK'][0][self.index[names]]

    def is_bulged(self, names):
        """Check if galaxies have bulge component or not.
        Args:
            names (list or str): name(s) of galaxies
        Returns:
            list: list of True and False.
        """
        if type(names) == list:
            bulged = []
            for n in names:
                bulged.append(self.data['BULGE'][0][self.index[n]])
            return bulged
        elif type(names) == str:
            return self.data['BULGE'][0][self.index[names]]

    def bar_dominating(self, name, band='R'):
        """Calculate the radius at which the bar dominates
        over the bulge.
        Args:
            name (str): name of the galaxy
        Returns:
            float: radius in arcsec
        """
        if self.data["BULGE"][0][self.index[name]] == False:
            return 0
        sersic = self.get_sersic(name, band)
        ferres = self.get_ferres(name, band)
        try:
            rbar = self.data['RBAR_'+band][0][self.index[name]]
            rcrit = self.estimate_crossing(
                sersic, ferres, steps=0.01, end=rbar)
            if rcrit == -1:
                return self.disk_dominating(name)
            r = fsolve(lambda x: sersic(x) - ferres(x), rcrit)
            # if r[0] == 0 or r[0] == rcrit:
            #     return self.disk_dominating(name)
            return r[0]
        except:
            print("Could not find Bar dominating radius for {:}".format(name))
            return 0.0

    def disk_dominating(self, name, band='R'):
        """Calculate the radius at which the luminosity by the disk dominates.
        Evaluated by comparing the contribution of the Sersic (Bar) profile with
        the profile of the disk luminosity.
        Args:
            name (str): name of the galaxy
            band (str): SDSS band
        Returns:
            float: radius"""
        if self.data["BULGE"][0][self.index[name]] == False:
            return 0
        sersic = self.get_sersic(name, band)
        disk = self.get_disk(name, band)
        if sersic(0) < disk(0):
            return 0
        try:
            rcrit = self.estimate_crossing(sersic, disk, steps=0.1, end=40)
            if rcrit == -1:
                print(name)
            r = fsolve(lambda x: sersic(x) - disk(x), rcrit)
            if r[0] == 0:
                print(
                    "Could not calculate critical radius! disk dominating for {:}".format(name))
            return r[0]
        except:
            print("Could not find Disk dominating radius for {:}".format(name))
            return 0.0

    def get_sersic(self, name, band):
        """Get the Sersic function of the bulge in [L0].
        Args:
            name (str): name of galaxy
            band (str): SDSS band
        Returns:
            function: the Sersic function."""
        Mser = self.data['MUE_' +
                         band][0][self.index[name]]  # effective surface brightness
        Iser = np.power(10, -0.4 * Mser)  # in L0 luminosity
        Reser = self.data['RE_'+band][0][self.index[name]]  # effective radius
        nser = self.data['N_'+band][0][self.index[name]]
        bser = 0.868 * nser - 0.142
        return lambda x: Iser * np.power(10, (-bser*(np.power((x/Reser), (1/nser)) - 1)))

    def get_ferres(self, name, band):
        """Get the Ferres function of the Bar in [L0]
        Args:
            name (str): name of galaxy
            band (str): SDSS band
        Returns:
            function: the Ferres function."""
        if not self.is_barred([name])[0]:
            raise Exception("Galaxy has no bar!")
        Mbar = self.data['MU0B_'+band][0][self.index[name]]
        Ibar = np.power(10, - 0.4 * Mbar)
        nbar = self.data['NBAR_'+band][0][self.index[name]]
        rbar = self.data['RBAR_'+band][0][self.index[name]]
        return lambda x: Ibar * np.power((1 - (x / rbar)**2), (nbar+0.5))

    def get_disk(self, name, band):
        """Get the luminosity disk profile in [L0].
        Args:
            name (str): name of galaxy
            band (str): SDSS band
        Returns:
            function: the disk luminosity function."""
        M = self.data['MU0_'+band][0][self.index[name]]
        I = np.power(10, -0.4 * M)
        h = self.data['HI_'+band][0][self.index[name]]
        ho = self.data['HO_'+band][0][self.index[name]]
        rbreak = self.data['RBREAK_'+band][0][self.index[name]]
        if self.data["DBREAK"][0][self.index[name]] == False:
            rbreak = 999

        def f1(x): return np.exp(-x/h)*np.heaviside(-x+rbreak, 0)
        def f2(x): return np.exp(-rbreak*(ho-h)/(h*ho)) * np.exp(-x/ho)
        return lambda x: I * (f1(x) + f2(x)*(1 - np.heaviside(-x+rbreak, 0)))

    def get_redshift(self, names):
        """Get the redshift of a list of galaxies.
        Args:
            name (list): names of galaxies
        Returns:
            list: list of the redshifts.
        """
        redshift = []
        for n in names:
            z = self.data["REDSHIFT"][0][self.index[n]]
            redshift.append(z)
        return redshift

    def estimate_crossing(self, f1, f2, steps=1, end=50):
        """Calculate the estimated value of x at which the functions f1 and f2 are equal.
        Args:
            f1 (function): function 1 of one argument x
            f2 (function): function 2 of one argument x
            steps (float): step size
            end (float)  : max value for x
        Returns:
            float
        """
        step = 0.0001
        sign = 0
        if f1(step) - f2(step) > 0:
            sign = +1
        else:
            sign = -1
        while step < end:
            val = f1(step) - f2(step)
            if val / abs(val) != sign:
                break
            step += steps
        if step >= end:
            # print(
            #     "Could not find Crossing of the functions in interval [{:},{:}]".format(steps, end))
            return -1
        else:
            return step

    def plot_profile(self, name, band='R'):
        xx = np.linspace(0, 50)
        disk = self.get_disk(name, band)
        if self.is_barred([name])[0]:
            ferres = self.get_ferres(name, band)
            plt.plot(xx, -2.5*np.log10(ferres(xx)), label='bar')
            r0 = self.bar_dominating(name)
            plt.plot(r0, -2.5*np.log10(ferres(r0)), 'o')
            print('Cossing at:', r0)
            if self.is_bulged(name):
                sersic = self.get_sersic(name, band)
                plt.plot(xx, -2.5*np.log10(sersic(xx)), label='bulge')
        elif self.data["BULGE"][0][self.index[name]]:
            sersic = self.get_sersic(name, band)
            plt.plot(xx, -2.5*np.log10(sersic(xx)), label='bulge')
            r0 = self.disk_dominating(name)
            plt.plot(r0, -2.5*np.log10(sersic(r0)), 'o')
            plt.plot(r0, -2.5*np.log10(disk(r0)), 'o')
        plt.plot(xx, -2.5*np.log10(disk(xx)), label='disk')
        plt.legend()
        plt.ylim(25, 16)
        plt.grid()
        plt.show()
        plt.close()

    def get_doubles(self):
        """Get galaxies with in ambiguous photometric decomposition.
        Args:
            None.
        Returns:
            names (list).
        """
        galax = self.data['GALAXY'][0]
        g = []
        abundance = []
        doubles = []
        for i in range(len(galax)):
            if galax[i] in g:
                abundance[g.index(galax[i])] += 1
                doubles.append(galax[i])
            else:
                g.append(galax[i])
                abundance.append(1)
        return doubles


if __name__ == "__main__":
    p = Photometrics("data/photometric_decomposition.fits")
    # print(p.is_barred(["NGC4210", "NGC6278"]))
    # print(p.is_disked(["NGC4210", "NGC6278"]))
    p.plot_profile("NGC2880")
    # g = 'NGC4210'
    # p.plot_profile(g)
    # print(p.bar_dominating(g))
    # if p.is_barred([g])[0]:
    #     print(p.bar_dominating(g))
    # else:
    #     print(p.disk_dominating(g))
