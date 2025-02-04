
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
                        print, warn )

def parseArgs( progname, arglist, return_parser=False ):

    progname_tool='nctool'
    descr=f"""

Script which can be used to create input files with reflections for McStas
crystalline sample components like PowderN and Single_crystal, based on NCrystal
cfg-strings (usually referring to NCMAT files with crystalline single-phase
materials).

Example invocations:

          $> {progname} "Y2SiO5_sg15_YSO.ncmat" --format=laz -o YSO.laz
          $> {progname} "Y2SiO5_sg15_YSO.ncmat" --format=lau -o YSO.lau

     It might be useful to set certain parameters in the cfg-string in addition
     to the input file name. For more information about available parameters,
     refer to https://github.com/mctools/ncrystal/wiki/CfgRefDoc, or run the
     cmd:

          $> {progname_tool} --doc

     For instance, to change the material temperature from room temperature, and
     to leave out reflections at dspacing<0.5Aa, run:

          $> {progname} "Y2SiO5_sg15_YSO.ncmat;dcutoff=0.5;temp=200K" --format=laz -o YSO_200K_largedsp.laz
          $> {progname} "Y2SiO5_sg15_YSO.ncmat;dcutoff=0.5;temp=200K" --format=lau -o YSO_200K_largedsp.lau

     To browse materials which are available in default NCrystal installations,
     refer to https://github.com/mctools/ncrystal/wiki/Data-library. To list
     actually available files in your installation, run:

          $> {progname_tool} --browse

"""
    from argparse import RawTextHelpFormatter
    import textwrap
    parser = create_ArgumentParser(prog = progname,
                                   description=descr,
                                   formatter_class=RawTextHelpFormatter)
    helpw = 60
    def wrap(t):
        return textwrap.fill(t,width=helpw)

    parser.add_argument('CFGSTR', type=str,
                        help=wrap("NCrystal cfg-string representing the material"
                              " to load. For the purposes of the present"
                              " script, it might be particularly relevant to"
                              " set the temp and dcutoff parameters"))
    parser.add_argument('--format','-f',default=None,choices=('laz','lau'),
                        help=wrap('What types of output files to produce. Choices'
                              ' are "lau" for .lau files for McStas'
                              ' PowderN or "laz" for .laz files for'
                              ' McStas Single_crystal.'))
    parser.add_argument('--output','-o',type=str,default=None,
                        help="Name of output file (default: stdout)")
    parser.add_argument('--quiet','-q',default=False,action='store_true',
                        help=wrap('Silence non-error output'
                              ' (automatic if --output=stdout).'))

    if return_parser:
        return parser

    _warned_outfile = False
    for i in range(len(arglist)):
        if arglist[i].startswith('--outfile'):
            arglist[i] = '--output' + arglist[i][len('--outfile'):]
            if not _warned_outfile:
                _warned_outfile = True
                warn('The --outfile option has been renamed to --output')


    args=parser.parse_args( arglist )

    if not args.format:
        parser.error('Please specify output type with --format,'
                     ' either --format=laz or --format=lau (run'
                     ' with --help for details)')

    #TODO: We could use a common function for the following logic (including the
    #stdout logic). Perhaps we should simply have a special "type" for the
    #argparse option? Also, in general --quiet and --force are connected to
    #--output.

    if args.output is not None and args.output!='stdout':
        import pathlib
        of = pathlib.Path(args.output)
        if of.is_dir():
            parser.error(f'Output destination is directory: {of}')
        if of.exists():
            parser.error(f'Output destination already exists: {of}')
        ofr = of.resolve().absolute()
        if not ofr.parent.exists() or not ofr.parent.is_dir():
            parser.error('Output destination is in non-existing'
                         f' directory: {ofr.parent}')
        args.output = ofr
    else:
        args.output = 'stdout'
    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)


@cli_entry_point
def main( progname, arglist ):

    #parse hidden option (needed for unit tests):
    override_prec = None
    override_opt = '--override-prec='
    hidden_args, standard_args = [], []
    for a in arglist:
        (hidden_args
         if a.startswith(override_opt)
         else standard_args).append( a )
    arglist = standard_args
    if hidden_args:
        override_prec = int(hidden_args[-1][len(override_opt):])

    args = parseArgs( progname, arglist )
    do_quiet = ( args.quiet or args.output == 'stdout' )

    def invoke():
        _main_impl(args,do_quiet,override_prec)

    if do_quiet:
        from ._common import modify_ncrystal_print_fct_ctxmgr
        with modify_ncrystal_print_fct_ctxmgr('block'):
            invoke()
    else:
        invoke()

def _main_impl( args, do_quiet, override_prec ):

    from . import mcstasutils
    full_info = True
    kwargs = dict(cfgstr = args.CFGSTR,
                  tgtformat = args.format,
                  verbose = full_info)
    if override_prec:
        kwargs['fp_format'] = f'%.{override_prec}g'

    content_iter = mcstasutils.cfgstr_2_hkl( **kwargs )
    #TODO: common helper function?
    if args.output is None or args.output == 'stdout':
        import sys
        for line in content_iter:
            sys.stdout.write(line + '\n')
    else:
        from ._common import write_text
        write_text( args.output,
                    ( '\n'.join( content_iter ) + '\n' ) )
        if not do_quiet:
            import pathlib
            print(f"Wrote {pathlib.Path(args.output).name}")
