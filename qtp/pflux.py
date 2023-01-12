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


if __name__ == '__main__': 
    from matplotlib import pyplot as plt
    from peb import PotentialEnergyBarrier
    from traco import TransCoeff
    from scipy import constants

    dalton2me = constants.atomic_mass/constants.m_e
    kelvin2au = 1 / (constants.value('Hartree energy')/constants.k)

    zvals = np.linspace(-5.0, 5.0, 30, endpoint=True)
    evals = 0.05 * np.exp(-(zvals**2)/4.0)

    PEB = PotentialEnergyBarrier(zvals, evals)
    TraCo = TransCoeff(PEB, dalton2me)
    pflux = PFlux(TraCo)

    for temperature in [100, 200, 300]:
        beta_au = 1/(temperature*kelvin2au)
        print(temperature, beta_au, pflux(beta_au))

    zs = np.linspace(-6, 6, 100, endpoint=True)
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)

    ax1.plot(zvals, evals, 'o-', label="input")
    ax1.plot(zs, PEB(zs), label='cubic splines')
    ax1.legend()

    ax2.plot(evals, [np.exp(TraCo(e)) for e in evals])

    plt.show()