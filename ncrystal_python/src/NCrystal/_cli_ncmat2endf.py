
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
    from .ncmat2endf import available_elastic_modes
    from .ncmat2endf import EndfMetaData
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
                        help=("overwrite existing file "
                              "if it already exists"))
    parser.add_argument('-e', '--elastic_mode',
                        help='approximation used for the elastic component',#fixme: consider docs, perhaps refer to njoy-ncrystal paper?
                        type=str, choices=available_elastic_modes,
                        default='scaled')

    #fixme: remove these, but add --metadata which takes json dict:
    parser.add_argument('--set_date_to_now',action='store_true',
                        help=('Set ENDF6 fields EDATE, DDATE and RDATE'
                              ' to current month and year.'))
    parser.add_argument('-m', '--mats',
                        help=('JSON dictionary containing material number '
                              'assignment for each element, e.g. \'{"C":37, '
                              '"H": 38}\''),
                        type=json.loads)

    #fixme: these should be marked as expert-only (also, add options controlling
    #emax, lasym)
    parser.add_argument('-t', '--temperatures',
                        nargs='+',
                        type=float,
                        help='''additional temperatures to process (fixme:
                        mention bad idea?)''')
    parser.add_argument('--smin',
                        help='set the minimum value of S(alpha, beta) stored '
                             'in MF7/MT4', type=float,
                        default=EndfMetaData().smin)

    #fixme: test the cli

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
    params = EndfMetaData()
    params.set_smin(args.smin)

    _ = ncmat2endf(args.input,
                   material_name=args.material_name,
                   endf_metadata=params,
                   temperatures=args.temperatures,
                   mat_numbers=args.mats,
                   elastic_mode=args.elastic_mode,
                   force_save=args.force,
                   set_date_to_now=args.set_date_to_now,
                   verbosity=args.verbosity)

