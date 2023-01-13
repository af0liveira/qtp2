import numpy as np
from scipy.interpolate import CubicSpline


class PotentialEnergyBarrier(object):
    """Potential energy barrier in function of z coordinate
    
    Parameters
    ----------
    zvalues_bohr: array-like
        Sequence of 1D coordinate values in bohr
    energies_hartree: array-like
        Sequence of energy values in hartree corresponding to the values
        in `zvalues_bohr`
    """

    def __init__(self, zvalues_bohr: list, energies_hartree: list): 
        self._zvalues_bohr = zvalues_bohr
        self._energies_hartree = energies_hartree 
        self._cubic_spline = CubicSpline(self._zvalues_bohr,
                                         self._energies_hartree,
                                         extrapolate=False)

    def __call__(self, z_bohr, der=0):
        """Returns the energy values associated to input coordinates

        Parameters
        ----------
        z_bohr: float or array-like
            Coordinate values in bohr
        der: int
            Order of the derivative

        Return
        ------
        float or array-like
            Energy values corresponding to the `z_bohr`

        Notes
        -----
        It is assumed that the potential energy function is converged to
        zero on the right and left sides. Thus, it is assumed that 
        the energy -- and its derivatives -- are zero outside the bounds
        of the cubic spline representing the potential energy function.
        """
        ret = self._cubic_spline(z_bohr, nu=der, extrapolate=False)
        return np.nan_to_num(ret)

    def get_coords_from_energy(self, energy_hartree: float) -> np.array: 
        """Return the coordinates corresponding to a given energy value

        Parameters
        ----------
        energy_hartree: float
            Potential energy value in hartree

        Return
        ------
        array-like
            List of all coordinates (in bohr) associated with the input 
            energy value
        """
        assert energy_hartree > -1e-6, \
                "Potential energy barrier must be positive!"
        return self._cubic_spline.solve(y=energy_hartree, extrapolate=True)

    def get_max_energy(self) -> float: 
        """Return the maximum potential energy value."""
        critical_pts = self._cubic_spline.derivative(nu=1).roots()
        return np.max(self.__call__(critical_pts))

