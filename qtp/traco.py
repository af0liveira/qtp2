import numpy as np
from scipy import integrate

class TransCoeff(object):
    """Transmission coefficient of a 1D potential energy barrier

    Parameters
    ----------
    peb: PotentialEnergyBarrier instance
        Describes the potential energy barrier in function of the 
        position coordinates
    
    mass_me: float
        Mass of the particle crossing the potential energy barrier, in 
        electron-mass units (not atomic mass units!)
    """

    def __init__(self, peb, pmass_me): 
        self.peb = peb
        self.pmass_me = pmass_me

    def __call__(self, E_hartree):
        """Return the natural logarithm of the transmission coefficient 
        T(E) for a particle with energy E

        Parameters
        ----------
        E_hartree: float
            Particle energy in hartree
        
        Returns
        -------
        float
            Value of ln(T(E))
        """
        def integrand(z):
            """Integration function."""
            U_hartree = self.peb(z)
            if U_hartree < E_hartree:
                return 0
            return np.sqrt(U_hartree - E_hartree)

        Umax = self.peb.get_max_energy()
        zvals_for_E = self.peb.get_coords_from_energy(E_hartree)

        z1 = np.min(zvals_for_E) 
        z2 = np.max(zvals_for_E)
        integral = integrate.quad(integrand, z1, z2, limit=100)

        return -2*np.sqrt(2*self.pmass_me) * integral[0]

