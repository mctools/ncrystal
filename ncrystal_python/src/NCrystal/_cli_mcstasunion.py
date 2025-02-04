
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
                        cli_entry_point,
                        print )

def parseArgs( progname, args, return_parser=False ):

    descr="""Output McStas-Union code in McStas instrument (.instr) format to
    define material with NAME based on specified NCrystal CFGSTR.  Providing the
    --split option will split NCrystal processes into physics types at the
    McStas-Union level, rather than internally in NCrystal."""

    parser = create_ArgumentParser(prog = progname,
                                   description=descr)
    parser.add_argument('union_mat_name', metavar='NAME', type=str,
                        help="Resulting NAME of McStas-Union material.")
    parser.add_argument('cfgstr', metavar='CFGSTR', type=str,
                        help="NCrystal cfg-string of the material.")
    parser.add_argument('--split','-s',action='store_true',
                        help=("Splits NCrystal processes"
                              "into physics types at the McStas-Union level,"
                              "rather than internally in NCrystal."))
    parser.add_argument("--output",'-o',default='',type=str,
                        help=("Output file name (defaults to stdout)."))
    parser.add_argument('--force','-f',action='store_true',
                        help=("Override output file if it already exists."))
    if return_parser:
        return parser
    args = parser.parse_args(args)
    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)

@cli_entry_point
def main( progname, args ):
    args = parseArgs( progname, args )
    from .mcstasutils import cfgstr_2_union_instrument_code

    code = cfgstr_2_union_instrument_code( cfgstr = args.cfgstr,
                                           name = args.union_mat_name,
                                           split_by_physics = args.split )
    if not args.output or args.output=='stdout':
        print( code )
    else:
        import pathlib
        p = pathlib.Path( args.output )
        if p.is_file() and not args.force:
            raise SystemExit(f'File already exists (use --force to overwrite): {p}')
        if not p.parent.is_dir():
            raise SystemExit(f'Output directory does not exist: {p.parent}')
        from ._common import write_text
        write_text(p,code)
        print(f"Wrote: {p}")
