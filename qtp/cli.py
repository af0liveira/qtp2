import argparse

cli_parser = argparse.ArgumentParser(
        prog='qtp2', 
        fromfile_prefix_chars='@',
        )

cli_parser.add_argument('datafile', metavar='<file>', type=str,
        help=("XY file describing the potential energy barrier U(z)."
              " The coordinate values are given in the X column in bohrs,"
              " whereas the corresponding U values are specified in the"
              " Y column in hartrees."))

cli_parser.add_argument('-m', dest='pmass_Da', metavar='<mass/Da>', type=float,
        default=1.0, 
        help="particle mass in daltons (1 Da = 1 g/mol). (Default: 1 Da)")

cli_parser.add_argument('-tstart', dest='tstart_K', metavar='<temperature/K>',
        type=float, default=300.0,
        help="Starting temperature in kelvins. (Default = 300 K")

cli_parser.add_argument('-tstop', dest='tstop_K', metavar='<temperature/K>', 
        type=float, required=False, 
        help=("Final temperature in kelvins. This is an optional argument, in"
             " case a range of temperatures is required for the calculation."
             " The `-tstep` option can be used to control the interval between"
             " consecutive temperature values."))

cli_parser.add_argument('-tstep', dest='tstep_K', metavar='<step/K>',
                        type=float, default=10,
                        help=("Inteval between consecutive temperatures when"
                              " requesting calculations over a range of"
                              " temperatures. (Default: 10 K)"))

cli_parser.add_argument('-eshift', dest='eshift_hartree',
        metavar='<energy/Eh>', type=float, default=0.0, 
        help=("Apply energy shift to the potential energy barrier."
              " The energy barrier should have its baseline at zero." 
              " (Default: 0)"))
