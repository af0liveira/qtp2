import os
import sys

import datetime

import numpy as np
import pandas as pd
from scipy import constants

sys.path.insert(0, os.path.dirname(__file__))
from cli import cli_parser
from peb import PotentialEnergyBarrier
from traco import TransCoeff
from pflux import PFlux
from arrhenius import get_arrhenius
sys.path.pop(0)

TITLE = 'Quantum Transport Properties 2 (v. 230112)'

dalton2me = constants.atomic_mass/constants.m_e
kelvin2au = 1 / (constants.value('Hartree energy')/constants.k)
angst2bohr = constants.angstrom/constants.value('Bohr radius')

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)


def qtp2():
    """QTP2 driver -- i.e., the main routine. ;)"""

    # Print output header
    # -------------------

    ruler = len(TITLE)*"="
    print(ruler, TITLE, ruler, sep='\n')
    print(datetime.datetime.now(), end='\n\n', flush=True)


    # Parse CLI options
    # -----------------

    cli_args = cli_parser.parse_args()
    
    pmass_Da = getattr(cli_args, 'pmass_Da')
    pmass_me = pmass_Da*dalton2me

    t0 = getattr(cli_args, 'tstart_K')
    tf = getattr(cli_args, 'tstop_K')
    if tf is None:
        tf = t0
    tstep = getattr(cli_args, 'tstep_K', 0)

    temperatures_K = np.arange(t0, tf+0.5*tstep, tstep)

    # Parse potential energy barrier
    # ------------------------------

    datafile = getattr(cli_args, 'datafile')
    with open(datafile, 'r') as fp:
        values = [float(v) for s in fp.readlines() if len(s)>0 and s[0]!='#'
                    for v in s.strip().split()]

    zvalues_bohr = np.array(values[0::2])
    energies_hartree = np.array(values[1::2]) \
                       + getattr(cli_args, 'eshift_hartree')
    peb = PotentialEnergyBarrier(zvalues_bohr, energies_hartree)

    # Build table of transmission coefficients
    # ----------------------------------------
    traco_calculator = TransCoeff(peb=peb, pmass_me=pmass_me)

    transco_df = pd.DataFrame(
            {
                r'z/bohr': zvalues_bohr, 
                r'U(z)/hartree': [peb(z) for z in zvalues_bohr], 
                r'ln[T(E)]': [traco_calculator(peb(z)) for z in zvalues_bohr],
            }
        )

    transco_df.insert(2, r'T(E)', np.exp(transco_df[r'ln[T(E)]']))
    transco_df.insert(0, r'z/angstrom', transco_df[r'z/bohr']/angst2bohr)
    transco_df.insert(0, r'pmass/m_e', pmass_me)
    transco_df.insert(0, r'pmass/Da', pmass_Da)

    print(transco_df, end='\n\n', flush=True)

    # Build table of particle fluxes
    # ------------------------------

    pflux_calculator = PFlux(traco_calculator)

    betas_au =  1/(temperatures_K*kelvin2au)
    jclass_jtunnel = [pflux_calculator(b) for b in betas_au]
    j_class = np.array([pair[0] for pair in jclass_jtunnel])
    j_tunnel = np.array([pair[1] for pair in jclass_jtunnel])

    j2k = np.sqrt(2*np.pi*pmass_me*betas_au)
    k_class = j2k*j_class
    k_tunnel = j2k*j_tunnel
    k_total = k_class + k_tunnel

    pflux_df = pd.DataFrame(
            {
                r'Temp/K': temperatures_K, 
                r'beta/(1/K)': 1/temperatures_K, 
                r'beta/au': betas_au,
                r'k_class': k_class,
                r'k_tunnel': k_tunnel,
                r'k_total': k_total,
                r'flux_class': j_class, 
                r'flux_tunnel': j_tunnel,
                r'flux_total': j_class+j_tunnel, 
            }
        )
    pflux_df.insert(3, r'pmass/m_e', pmass_me)
    pflux_df.insert(3, r'pmass/Da', pmass_Da)

    print(pflux_df, end='\n\n', flush=True)

    # Get Arrhenius parameters
    # ------------------------

    try:
        eactvals_class, coeffs_class = get_arrhenius(betas_au, np.log(k_class))
        eactvals_tunnel, coeffs_tunnel = get_arrhenius(betas_au, np.log(k_tunnel))
        eactvals_total, coeffs_total = get_arrhenius(betas_au, np.log(k_total))
    except ValueError:
        print("Skipping calculation of Arrhenius properties."
              " Not enough points for interpolation.", end='\n\n')
    else:
        arrhenius_df = pd.DataFrame(
                {
                    r'Temp/K': temperatures_K, 
                    r'beta/(1/K)': 1/temperatures_K, 
                    r'beta/au': betas_au,
                    r'Eact_class/au': eactvals_class,
                    r'Eact_tunnel/au': eactvals_tunnel,
                    r'Eact_total/au': eactvals_total,
                    r'Afactor_class': coeffs_class,
                    r'Afactor_tunnel': coeffs_tunnel,
                    r'Afactor_total': coeffs_total,
                }
            )
        arrhenius_df.insert(3, r'pmass/m_e', pmass_me)
        arrhenius_df.insert(3, r'pmass/Da', pmass_Da)

        print(arrhenius_df, end='\n\n', flush=True)

    # The end
    # -------
    print(datetime.datetime.now(), flush=True)
    print("=== THE END ===")


if __name__ == '__main__': 
    qtp2()
