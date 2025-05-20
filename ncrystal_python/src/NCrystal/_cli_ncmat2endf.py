
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

from ._cliimpl import ( create_ArgumentParser,
                        cli_entry_point )

def _parseArgs( progname, arglist, return_parser=False ):
    from argparse import RawTextHelpFormatter
    #FIXME: check before final merge
    #NOTE: Keep the description below synchronised with doc-string of the
    #ncmat2endf python API function:
    descr="""

Script for creating a set of ENDF-6 thermal scattering files from a .ncmat
file. Basic paramters are supported as arguments, but additional parameters
for the ENDF-6 can be set by using the Python API and pasing a custom
EndfMetaData object.

The script allows to handle multiple temperatures in one ENDF-6 file, but this
is not recommended, because NCrystal computes an optimal (alpha, beta) grid for
each material and temperature, while the ENDF format imposes the same grid on
all temperatures.

Ths script uses the endf-parserpy package from IAEA to format and check the
syntax of the ENDF-6 file.

G. Schnabel, D. L. Aldama, R. Capote, "How to explain ENDF-6 to computers:
A formal ENDF format description language", arXiv:2312.08249,
DOI:10.48550/arXiv.2312.08249

"""
    parser = create_ArgumentParser(prog = progname,
                                   description=descr,
                                   formatter_class=RawTextHelpFormatter)

    parser = create_ArgumentParser(
        description=descr,
        formatter_class=RawTextHelpFormatter)
    from .ncmat2endf import ( available_elastic_modes,
                              default_smin_value,
                              default_emax_value )
    import json

    parser.add_argument('input',
                        help='NCMAT cfg string to convert')
    parser.add_argument('-m', '--material_name',
                        help=('name of the material to be processed. '
                              'ENDF files will be named '
                              'tsl_element_in_name.endf for compounds '
                              'or tsl_element.endf for elements. E.g. '
                              'tsl_H_in_CH2.endf or tsl_Cu.endf'))
    parser.add_argument('-v', '--verbosity',
                        help='controls how verbose should be the output',
                        type=int, default=1)
    parser.add_argument('-f', '--force',action='store_true',
                        help=('overwrite existing file '
                              'if it already exists'))
    parser.add_argument('-e', '--elastic_mode',
                        help=('approximation used for the elastic component.',
                              ' An explanation of these modes is given in:',
                              'https://doi.org/10.1016/j.nima.2021.166227'),
                        type=str, choices=available_elastic_modes,
                        default='scaled')
    parser.add_argument('--metadata',
                        help=('JSON dictionary containing ENDF-6 metadata'),
                        type=json.loads)
    parser.add_argument('--set_date_to_now',action='store_true',
                        help=('Set ENDF6 fields EDATE, DDATE and RDATE'
                              ' to current month and year.'))

    parser.add_argument('-t', '--temperatures',
                        nargs='+',
                        type=float,
                        help=' (EXPERT ONLY)'
                             ' additional temperatures to process.'
                             ' If additional temperatures are provided,'
                             ' ncmat2endf will interpolate the S(a,b) into'
                             ' a grid for the first temperature. This should'
                             ' be done with care. It is preferred to run'
                             ' each temperature independently using the temp='
                             ' keyword in the cfg string.'   )
    parser.add_argument('--smin',
                        help=' (EXPERT ONLY)'
                             ' set the minimum value of S(alpha, beta) stored '
                             'in MF7/MT4', type=float,
                        default=default_smin_value)
    parser.add_argument('--emax',
                        help=' (EXPERT ONLY)'
                             ' maximum energy for the scatterig kernel',
                        type=float,
                        default=default_emax_value)
    parser.add_argument('--asymmetric_sab',action='store_true',
                        help=(' (EXPERT ONLY)'
                             ' store S(a,b) without the detailed balance'
                             ' factor exp(beta/2).'))
    parser.add_argument('--total_sab',action='store_true',
                        help=(' (EXPERT ONLY)'
                             ' store S(a,b) branches for positive and '
                              ' negative beta'))

    if return_parser:
        return parser
    args=parser.parse_args(arglist)
    return args

def create_argparser_for_sphinx( progname ):
    return _parseArgs(progname,[],return_parser=True)

@cli_entry_point
def main( progname, arglist ):
    from .ncmat2endf import EndfMetaData, ncmat2endf
    args = _parseArgs( progname, arglist )
    metadata = EndfMetaData()
    if args.set_date_to_now:
        metadata.set_all_dates_as_now()
    lasym = 0
    if args.total_sab:
        lasym = 1
    if args.asymmetric_sab:
        lasym += 2

    _ = ncmat2endf(args.input,
                   material_name=args.material_name,
                   endf_metadata=metadata,
                   temperatures=args.temperatures,
                   elastic_mode=args.elastic_mode,
                   force_save=args.force,
                   smin=args.smin,
                   emax=args.emax,
                   lasym=lasym,
                   verbosity=args.verbosity)

