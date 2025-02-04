
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

def _parseArgs( default_debye_temp, progname, arglist, return_parser=False ):
    from argparse import RawTextHelpFormatter

    fdt = default_debye_temp
    #NOTE: Keep the description below synchronised with doc-string of the
    #hfg2ncmat python API function:
    descr=f"""

Script which can be used to generate NCMAT files for hydrogen-rich amorphous
materials, in which the hydrogen atoms are bound to certain standard functional
groups (e.g. carbohydrates, polyimides, polymers, ...). Based on the material's
density, (empirical) chemical formula, and the specification of hydrogen
bindings in terms of standard functional groups, an NCMAT file is generated. In
this NCMAT file, non-hydrogen atoms are treated with a simplistic model
(idealised Debye model of phonon vibrations, assuming a Debye temperature of
{fdt}K for all atoms unless --debyetemp is specified), and the hydrogen atoms
are treated with a proper phonon density of state (VDOS) curve, which is
constructed based on the provided binding specifications. This is done using an
idea (and VDOS curves from) the following publication:

  "Thermal neutron cross sections of amino acids from average contributions
  of functional groups", G. Romanelli, et. al., J. Phys.: Condens. Matter,
  (2021). doi:10.1088/1361-648X/abfc13

Example of valid spec strings for the --spec parameter. Note only the hydrogen
bindings are important here:

     "1xCHali+2xCH2+8xCHaro+2xCH3+1xOH" (20 hydrogen atoms)
     "1xNH+3xSH+2xCH2" (8 hydrogen atoms).

List of valid bindings which can be used in --spec:

      CHali (1 hydrogen atom, aliphatic binding)
      CHaro (1 hydrogen atom, aromatic binding, e.g. on phenyl group).
      CH2 (2 hydrogen atoms)
      CH3 (3 hydrogen atoms)
      NH (1 hydrogen atom)
      NH2 (2 hydrogen atoms)
      NH3 (3 hydrogen atoms)
      OH (1 hydrogen atom)
      SH (1 hydrogen atom)

Note that as only hydrogen atoms will have a realistic VDOS curve, and as the
incoherent approximation is employed, the realism of the resulting material
modelling will be higher for materials with more hydrogen atoms. As a metric of
this, the script prints out how much of the total scattering cross section (sum
of sigma_coh+sigma_incoh for all atoms), is due to the sigma_incoh contribution
from hydrogen atoms.

"""
    parser = create_ArgumentParser(prog = progname,
                                   description=descr,
                                   formatter_class=RawTextHelpFormatter)
    parser.add_argument("--output",'-o',default='autogen.ncmat',type=str,
                        help=("Output file name (defaults to autogen.ncmat)."
                              " Can be stdout."))
    parser.add_argument('--force',action='store_true',
                        help=("Will overwrite existing file "
                              "if it already exists."))
    parser.add_argument("--spec",'-s',metavar='SPEC',type=str,required=True,
                        help="Hydrogen binding specification (see above).")
    parser.add_argument("--formula",'-f',metavar='FORMULA',
                        type=str,
                        required=True,
                        help=("Chemical formula such as C20H18O3 (only"
                              " the relative ratios matter)."))
    parser.add_argument("--density",'-d',metavar='DENSITY',
                        type=float,
                        required=True,
                        help="Material density in g/cm3.")
    parser.add_argument("--debyetemp",metavar='VALUE',
                        type=float,
                        default=default_debye_temp,
                        help=("Debye temperature (kelvin) of non-hydrogen"
                              f" atoms (default is {default_debye_temp})."))
    parser.add_argument("--title",type=str,required=True,
                        help=('Title string of material (will be placed as'
                              ' comment near top of output file). Use \\n '
                              'for line-breaks.'))
    parser.add_argument('--notrim',action='store_true',
                        help="No trimming of resulting VDOS curve.")
    if return_parser:
        return parser
    args=parser.parse_args(arglist)
    return args

def create_argparser_for_sphinx( progname ):
    from .hfg2ncmat import _default_debye_temp
    return _parseArgs(_default_debye_temp(),progname,[],return_parser=True)

@cli_entry_point
def main( progname, arglist ):
    import pathlib
    from .hfg2ncmat import _default_debye_temp, hfg2ncmat
    from ._common import write_text as nc_write_text
    args = _parseArgs( _default_debye_temp(), progname, arglist )
    do_stdout = args.output=='stdout'
    try:
        ncmat =  hfg2ncmat( spec = args.spec,
                            formula = args.formula,
                            density = args.density,
                            title = args.title,
                            debyetemp = args.debyetemp,
                            verbose = not do_stdout,
                            notrim = args.notrim )
    except RuntimeError as e:
        raise SystemExit('Error: %s'%str(e))
    if do_stdout:
        print(ncmat)
        return
    outfile = pathlib.Path(args.output)
    if outfile.exists() and not args.force:
        raise SystemExit('Error: output file already exists'
                         ' (run with --force to overwrite)')
    if not outfile.parent.is_dir():
        raise SystemExit('Error: output directory does not exist:'
                         f' { outfile.parent }')
    nc_write_text(outfile,ncmat)
    print(f"Wrote: {outfile}")
