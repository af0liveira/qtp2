import os
import sys

import datetime

import numpy as np
from scipy import constants

sys.path.insert(0, os.path.dirname(__file__))
from cli import cli_parser
from peb import PotentialEnergyBarrier
from traco import TransCoeff
from pflux import PFlux
sys.path.pop(0)

TITLE = 'Quantum Transport Properties 2 (v. 230112)'

dalton2me = constants.atomic_mass/constants.m_e
kelvin2au = 1 / (constants.value('Hartree energy')/constants.k)
angst2bohr = constants.angstrom/constants.value('Bohr radius')


def main():

    # Print the output header
    ruler = len(TITLE)*"="
    print(ruler, TITLE, ruler, sep='\n')
    print(datetime.datetime.now(), end='\n\n', flush=True)

    # Parse CLI options
    cli_args = cli_parser.parse_args()
    
    pmass_Da = getattr(cli_args, 'pmass_Da')
    pmass_me = pmass_Da*dalton2me
    print(f"Particle's mass: {pmass_Da:.2f} g/mol = {pmass_me:.12f} m_e")

    t0 = getattr(cli_args, 'tstart_K')
    tf = getattr(cli_args, 'tstop_K')
    if tf is None:
        tf = t0
    tstep = getattr(cli_args, 'tstep_K', 0)
    temperatures_K = np.arange(t0, tf+0.5*tstep, tstep)
    print(f"Temperature/K:", 
          ', '.join([f"{t:.2f}" for t in temperatures_K]))
    
    # Parse potential energy barrier
    datafile = getattr(cli_args, 'datafile')
    with open(datafile, 'r') as fp:
        values = [float(v) for s in fp.readlines() if len(s)>0 and s[0]!='#'
                    for v in s.strip().split()]

    zvalues_bohr = np.array(values[0::2])
    energies_hartree = np.array(values[1::2]) \
                       + getattr(cli_args, 'eshift_hartree')
    #energies_hartree -= min(energies_hartree)

    # Calculate transmission coefficients
    peb = PotentialEnergyBarrier(zvalues_bohr, energies_hartree)
    transcoeff = TransCoeff(peb=peb, pmass_me=pmass_me)

    header1 = '  '.join(
            [
                '{:^12s}'.format('z/angstrom'),
                '{:^12s}'.format('z/bohr'),
                '{:^12s}'.format('U(z)/hartree'),
                '{:^12s}'.format('T(E)'),
                '{:^12s}'.format('ln[T(E)]'),
            ]
    )
    ruler = len(header1)*'-'
    print(ruler, header1, ruler, sep='\n', flush=True)

    #zvals_to_show = np.linspace(min(zvalues_bohr), max(zvalues_bohr), num=50,
    #                            endpoint=True)
    #for z in zvals_to_show:
    for z in zvalues_bohr:
        E = peb(z)
        ln_traco = transcoeff(E)
        traco = np.exp(ln_traco)
        tab1_row = '  '.join(
            [
                '{:>12.6f}'.format(z/angst2bohr),
                '{:>12.6f}'.format(z),
                '{:>12.6f}'.format(E),
                '{:>12.6e}'.format(traco),
                '{:>12.6f}'.format(ln_traco),
            ]
        )
        print(tab1_row, flush=True) 

    print(ruler)
    print(datetime.datetime.now(), end='\n\n', flush=True)

    # Calculate particle flux
    pflux = PFlux(transcoeff)

    header2 = '  '.join(
        [
            '{:^10s}'.format("Temp/K"),
            '{:^10s}'.format("beta/(1/K)"),
            '{:^10s}'.format("beta/au"),
            '{:^14s}'.format("k_classic"),
            '{:^14s}'.format("k_tunnel"),
            '{:^14s}'.format("k_total"),
            '{:^14s}'.format("flux_classic"),
            '{:^14s}'.format("flux_tunnel"),
            '{:^14s}'.format("flux_total"),
        ]
    )
    ruler = len(header2)*'-'
    print(ruler, header2, ruler, sep='\n', flush=True)

    for t in temperatures_K:
        beta_au = 1/(t*kelvin2au)

        j_class, j_tunnel = pflux(beta_au)
        j_total = j_class + j_tunnel

        j2k = np.sqrt(2*np.pi*pmass_me*beta_au)
        k_class = j_class*j2k
        k_tunnel = j_tunnel*j2k
        k_total = j_total*j2k

        tab2_row = '  '.join(
            [
                '{:>10.2f}'.format(t),
                '{:>10.2e}'.format(1/t),
                '{:>10.2f}'.format(beta_au),
                '{:>14.6e}'.format(k_class),
                '{:>14.6e}'.format(k_tunnel),
                '{:>14.6e}'.format(k_total),
                '{:>14.6e}'.format(j_class),
                '{:>14.6e}'.format(j_tunnel),
                '{:>14.6e}'.format(j_total),
            ]
        )
        print(tab2_row, flush=True)

    print(ruler)
    print(datetime.datetime.now(), end='\n\n', flush=True)

if __name__ == '__main__': 
    main()