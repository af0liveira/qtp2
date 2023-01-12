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


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from peb import PotentialEnergyBarrier
    from scipy import constants

    dalton2me = constants.atomic_mass/constants.m_e

    zvals = np.linspace(-5.0, 5.0, 30, endpoint=True)
    evals = 0.05 * np.exp(-(zvals**2)/4.0)

    PEB = PotentialEnergyBarrier(zvals, evals)
    TraCo = TransCoeff(PEB, dalton2me)

    zs = np.linspace(-6, 6, 100, endpoint=True)
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)

    ax1.plot(zvals, evals, 'o-', label="input")
    ax1.plot(zs, PEB(zs), label='cubic splines')
    ax1.legend()

    ax2.plot(evals, [np.exp(TraCo(e)) for e in evals])

    plt.show()