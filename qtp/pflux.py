import numpy as np
from scipy import integrate

class PFlux(object):
    """Particle flux (current) through a 1D potential energy barrier

    Parameters
    ----------
    trans_coeff: TransCoeff object
        Transmission coefficient for the potential energy barrier
    """

    def __init__(self, trans_coeff):
        self.trans_coeff = trans_coeff
        self.peb = trans_coeff.peb
        self.pmass_me = trans_coeff.pmass_me
        self.Umax = self.peb.get_max_energy()

    def __call__(self, beta_au):
        """Return the classical and quantum components of the flux for a
        given reverse temperature

        Parameters
        ----------
        beta_au: float
            Reverse temperature in atomic units

        Return
        ------
        2-tuple
            Classical and quantum flux components, respectively. 
            The total flux is given by the sum of these two components. 
        """
        def integrand(E):
            """Integration function"""
            return np.exp(self.trans_coeff(E)) * np.exp(-beta_au*E)

        E1 = 0
        E2 = self.Umax
        integral = integrate.quad(integrand, E1, E2, limit=100)
        pflux_qm = np.sqrt(beta_au/(2*np.pi*self.pmass_me)) * integral[0]

        pflux_class = 1/np.sqrt(2*np.pi*self.pmass_me*beta_au) \
                      * np.exp(-beta_au*self.Umax)

        return pflux_class, pflux_qm

