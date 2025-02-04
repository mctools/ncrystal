
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

def parseArgs( progname, arglist, return_parser=False ):
    import textwrap
    helpw = 60
    descrw = helpw + 22
    descr=( textwrap.fill(
        """Script which can be used to generate NCMAT files from crystal
structures in local or remote CIF files. The script can optionally also compare
existing NCMAT files with online structures indicated in their comments.

It is essentially a command-line front-end for the functionality provided by the
NCrystal.cifutils Python module. Behind the scenes, CIF data is loaded and
processed with the gemmi and spglib modules, which must be available on the
system. They can be either installed via PyPi ("python3 -mpip install gemmi
spglib") or conda ("conda install -c conda-forge gemmi spglib").""",descrw)+'\n\n'
            + textwrap.fill(
                """IMPORTANT NOTE: It is important to realise that while most
CIF files contain crystal structure information, most of them lack appropriate
information about atomic dynamics (i.e. atomic displacements, Uiso, Biso, mean
squared displacements, Debye temperature, phonon DOS, ...). Even when available
as a Uiso (or Biso) value in the CIF data, the CIF data usually does not specify
which temperature the Uiso value is associated with, making it meaningless. The
--uisotemp parameter can be used to manually specify the temperature at which
any Uiso in the file is valid (it might of course require the user to perform a
bit of detective work to figure out the correct value). In other cases, this
script will simply assign each atom a Debye temperature value of 300K and use
that for dynamics (the --debyetemp flag can be used to modify this value if
desired). For more detailed options, one can either manually edit the generated
NCMAT data, or use the Python API (specifically the CIFAnalyser class from the
NCrystal.cifutils module) which offers more options and flexibility than the
present script.  """,descrw)+f"""\n\nExample invocations:

          $> {progname} myfile.cif
          $> {progname} myfile.cif --output=myfile.ncmat
          $> {progname} myfile.cif --output=stdout
          $> {progname} myfile.cif --uisotemp=25
          $> {progname} myfile.cif --debyetemp=500
          $> {progname} mpid::127
          $> {progname} codid::9008460 --debyetemp=412

""")

    from argparse import RawTextHelpFormatter

    parser = create_ArgumentParser( prog = progname,
                                    description=descr.strip()+'\n',
                                    usage=(f'{progname} CIFFILE [-o OUTPUT] [--debyetemp VAL] [--uisotemp VAL]\n'
                                           +' '*(len(progname)+8)+'[<<other advanced options described below>>]'),
                                    formatter_class=RawTextHelpFormatter )
    def wrap(t):
        return textwrap.fill(t,width=helpw)

    #CIFFILE arg use nargs='?' if --valplots detected, to avoid spurious error:
    parser.add_argument('CIFFILE', type=str,nargs=('?' if any( a.startswith('--val') for a in arglist ) else 1 ),
                        help=wrap('CIF file to use as source. This can either be a path '
                                  'to a physical file, or an online database ID, specified '
                                  'using a syntax like "codid::1000055" or "mpid::127". Do not specify with --valplot.'))
    parser.add_argument('--output','-o',default=None,
                        help=wrap('Name of output NCMAT file (specify "stdout" to print to stdout). If not'
                                  ' supplied, a name will be automatically generated. Alternatively, this'
                                  ' can be used with --valplot to force output to a pdf file'))
    value_fallback_debye_temp = 300.0
    parser.add_argument("--debyetemp",default=value_fallback_debye_temp, type=float,metavar='VALUE',
                        help=wrap( 'Fall-back Debye temperature value (kelvin) to use for all'
                                   f' atoms (default: {value_fallback_debye_temp}K)' ))
    parser.add_argument("--uisotemp",default=None, type=float,metavar='VALUE',
                        help=wrap( 'Temperature value (kelvin) associated with any atomic displacement (Uiso/Biso) values'
                                   ' in the input. Without this parameter, Uiso/Biso values in the input will be ignored.' ))
    parser.add_argument('--quiet','-q',default=False,action='store_true',
                        help=wrap('Silence non-error output (automatic if --output=stdout).'))
    parser.add_argument('--noval',default=False,action='store_true',
                        help=wrap('Do not validate generated NCMAT data by trying to load it.'))
    parser.add_argument("--remap",metavar='COMP',default=[],type=str,action='append',
                        help=wrap('Redefine the composition of elements in the input CIF data, using'
                                  ' a syntax like the one in @ATOMDB sections of NCMAT files.'
                                  ' Examples: "H is D" and "B is 0.9 B10 0.1 B11". Colons can'
                                  ' be used in place of spaces if desired ("H:is:D").'))
    parser.add_argument("--atomdata",default=[],type=str,action='append',
                        help=wrap('Provide or override atom data of particular elements and isotopes, using'
                                  ' a syntax like the one in @ATOMDB sections of NCMAT files. '
                                  'Example: "Si 28.09u 4.1491fm 0.004b 0.171b". Colons can'
                                  ' be used in place of spaces if desired.'))
    parser.add_argument('--dynamics','-d',default=None,metavar='FILE',
                        help=wrap('Data file from which to copy over dynamics information (matching up atoms'
                                  ' based on Z-values). If this argument is not provided, or if a '
                                  'required element is not present in the indicated file, '
                                  'dynamics will be created with a fall-back Debye temperature (cf. --debyetemp).'))

    parser.add_argument('--nomerge',default=False,action='store_true',
                        help=wrap('Do not use same label for different non-equivalent positions of the same atom species.'))
    parser.add_argument('--showcif',default=False,action='store_true',
                        help=wrap('Simply show the input CIF data (mostly useful with mpid:: or codid:: sources).'))
    parser.add_argument("--valplot",nargs='*',metavar="FILE",
                        help=wrap('List of existing NCMAT files which should be compared'
                                  ' to crystal structures directly taken from online DBs (thus hopefully verifying'
                                  ' that online DB references mentioned in comments in those files are'
                                  ' actually compatible). This will launch interactive matplotlib plots, unless'
                                  ' the --output option is used to request pdf output (either as -opdf or '
                                  '-osomename.pdf, where the former will autogenerate a suitable name for the output PDF file'))

    parser.add_argument("--spacegroup",default=None,type=str,
                        help=wrap('Override the space group, disregarding the spacegroup embedded in the file. '))
    parser.add_argument('--via-ase',default=False,action='store_true',
                        help=wrap('Load input file via ASE (this can be useful to support non-CIF input files).'))


    if return_parser:
        return parser

    args = parser.parse_args(arglist)
    ll=[]
    for c in args.remap:
        p=c.replace(':',' ').split()
        if not len(p)>=3 or p[1]!='is':
            parser.error('invalid --remap syntax in "%s"'%c)
        ll.append( (p[0],' '.join(p[2:]) ) )
    args.remap = ll

    ll=[]
    for c in args.atomdata:
        p = c.replace(':',' ').split()
        if not len(p)==5 or not p[1].endswith('u') or not p[2].endswith('fm') or not p[3].endswith('b') or not p[4].endswith('b'):
            parser.error('invalid --atomdata syntax in "%s"'%c)
        ll.append( (p[0],' '.join(p[1:]) ) )
    args.atomdata = ll

    if (args.via_ase or args.showcif) and args.valplot:
        parser.error('Conflicting options')
    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)


@cli_entry_point
def main( progname, arglist ):
    args = parseArgs( progname, arglist )
    do_quiet = ( args.quiet or args.output == 'stdout' or args.showcif )

    if do_quiet:
        from ._common import modify_ncrystal_print_fct_ctxmgr
        with modify_ncrystal_print_fct_ctxmgr('block'):
            _main_impl(args,do_quiet)
    else:
        _main_impl(args,do_quiet)

def _main_impl( args, do_quiet ):

    from . import cifutils as nc_cifutils
    from . import _ncmatimpl as nc_ncmatimpl
    from . import _common as nc_common
    import pathlib


    #Trigger gemmi/spglib import error already here (and with nicer SystemExit):
    if not args.showcif:
        nc_cifutils._import_gemmi( sysexit = True )
        nc_ncmatimpl._import_spglib( sysexit = True )

    if args.via_ase:
        nc_ncmatimpl._import_ase( sysexit = True )

    with nc_common.WarningSpy( block = do_quiet):
        if args.valplot:
            _unit_test = nc_common.ncgetenv_bool('CIF2NCMAT_UNITTEST_NOPLOT')
            pdf_target = None
            if args.output and args.output.lower().endswith('pdf'):
                if args.output.lower() in ('pdf','.pdf'):
                    #auto-detect name:
                    if len(args.valplot)==1:
                        pdf_target = pathlib.Path(args.valplot[0]).name + '.pdf'
                elif args.output.endswith('.pdf'):
                    pdf_target = args.output
                if not pdf_target:
                    pdf_target = 'ncrystal_onlinedb_validate.pdf'
                assert len(pdf_target)>4 and pdf_target.endswith('.pdf')
            if _unit_test:
                print( 'UNITTEST: Would invoke nc_cifutils.produce_validation_plots with pdf_target=',pdf_target )
                return
            nc_cifutils.produce_validation_plots( args.valplot,
                                                  verbose_lbls = False,
                                                  pdf_target = pdf_target )
            return
        assert len(args.CIFFILE)==1
        cifsrc = nc_cifutils.CIFSource( args.CIFFILE[0], allow_fail = True )
        if args.via_ase:
            if not cifsrc.filepath:
                ase_input_data = cifsrc.load_data( quiet = do_quiet )
            else:
                ase_input_data = cifsrc.filepath

            if cifsrc.is_remote:
                #ase_input_data = cifsrc.load_data( quiet = do_quiet )
                ase_format = 'cif'
            else:
                #ase_input_data = args.CIFFILE[0]
                ase_format = None
            ase_output_cifdata = nc_ncmatimpl._cifdata_via_ase( ase_input_data,
                                                                ase_format = ase_format,
                                                                quiet = do_quiet )
            cifsrc = nc_cifutils.CIFSource( ase_output_cifdata )
        if cifsrc.invalid:
            raise SystemExit('Failed to load input')

        if args.showcif:
            import sys
            sys.stdout.write( cifsrc.load_data() )
            return

        lc = nc_cifutils.CIFLoader( cifsrc,
                                    merge_equiv = not args.nomerge,
                                    quiet = do_quiet,
                                    override_spacegroup = args.spacegroup )
        composer = lc.create_ncmat_composer( top_comments=['<<disableautogennotice>>'],
                                             uiso_temperature = args.uisotemp,
                                             quiet = do_quiet,
                                             remap = args.remap,
                                             skip_dyninfo = bool(args.dynamics) )

        if args.debyetemp and args.debyetemp > 0.0:
            composer.allow_fallback_dyninfo( args.debyetemp )
        if args.dynamics:
            composer.transfer_dyninfo_objects( args.dynamics )

        for elem, data in args.atomdata:
            composer.update_atomdb( elem, data )

        out, metadata = composer.create_ncmat( meta_data = True)
        autofn = nc_cifutils._suggest_filename( metadata, lc )

    fn = args.output or autofn

    if not args.noval:
        if not do_quiet:
            print("Verifying that resulting ncmat data can be loaded")
        from .core import directMultiCreate
        directMultiCreate(out,'vdoslux=0;dcutoff=0.2')

    if fn=='stdout':
        import sys
        sys.stdout.write(out)
    else:
        if not do_quiet:
            print(f"Writing {fn}")
        nc_common.write_text(fn,out)
