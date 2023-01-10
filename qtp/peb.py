import numpy as np

class PotentialEnergyBarrier(object):
    """Onedimensional potential energy barrier in function of spatial coordinate z""" 

    def __init__(self, zvalues_bohr: list, energies_hartree: list): 
        raise NotImplementedError("Class not implemented yet!")

    def __call__(self, z_bohr: float, der=0: int) -> float:
        """Return the potential energy or its derivative for a given coordinate. 

        Parameters
        ----------
        z_bohr: float
            Coordinate in bohrs
        der: int
            Order of the derivative

        Return
        ------
        float
            Value of the potential energy at position z_bohr
        """
        raise NotImplementedError("Evaluation of PotentialEnergyBarrier not implemented yet!") 

    def get_coords_from_energy(self, energy_hartree: float) -> np.array: 
        """Return coordinates for a given potential energy value

        Parameters
        ----------
        energy_hartree: float
            Potential energy value in hartree

        Return
        ------
        array-like
            List of all coordinates (in bohr) associated with the input energy value
        """
        raise NotImplementedError("Finding coordinates for given potential energy is not implemented yet!")

    def get_max_energy(self) -> float: 
        """Return the maximum potential energy value"""
        raise NotImplementedError("Finding maximum energy is not yet possible.")