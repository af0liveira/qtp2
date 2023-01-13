import numpy as np
from scipy.interpolate import CubicSpline


def get_arrhenius(betavals_au, lnkvals):
    """Calculate Arrhenius activation energy and coefficient

    Parameters
    ----------
    betavals_au: array_like
        Array of reciprocal temperatures in atomic units.

    lnkvals: array_like
        Natural logarithms of the equilibrium constants.

    Returns
    -------
    e_act: float
        Arrhenius activation energy.
    coeff: float
        Arrhenius coefficient. 
    """
    # Construct a cubic spline representing the Arrhenius equation
    # ln(k) = ln(A) - beta*E_a 
    # with ln(A) := 0
    x, y = zip(*sorted(zip(betavals_au, lnkvals)))  # sort points by `x`
    arrhenius_cs = CubicSpline(x, y, extrapolate=False)

    # Obtain the activation energies from the derivative wrt. beta
    # d[ln(k)]/d(beta) = -E_a
    activation_energies = -arrhenius_cs(betavals_au, nu=1)

    # Obtain the pre-exponential factors A as
    # A = k * exp(beta*E_a)
    kvals = np.exp(arrhenius_cs(betavals_au))
    prefactors = kvals * np.exp(betavals_au*activation_energies)

    return activation_energies, prefactors
