
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

longopt_metadata = '--mdata'
metavar_elastic = 'MODE'
metavar_matname = 'NAME'
metavar_metadata = 'DATA'
longopt_elastic = '--elas'
longopt_matname = '--name'
longopt_datenow = '--now'
longopt_othertemps = '--othertemps'

#Examples here so they can be unit tested:
examples = [
    ['Al_sg225.ncmat;temp=350K'],
    ['Si_sg227.ncmat;temp=293.6K','-m','MATNUM:Si:99',longopt_datenow],
    ['ZnO_sg186_ZincOxide.ncmat;temp=293.15K','-n','ZnO',
     '-e','scaled','-m','MATNUM:Zn:101,O:102'],
    ['Bi_sg166.ncmat;comp=inelas;temp=77K','-m','AUTH:J. Doe']
]

def _parseArgs( progname, arglist, return_parser=False ):
    from .ncmat2endf import ( available_elastic_modes,
                              default_smin_value,
                              default_emax_value )
    from ._common import print
    from argparse import RawTextHelpFormatter
    import textwrap
    import json
    import shlex

    helpw = 60
    descrw = helpw + 22
    descr_sections = [
        """Script for creating a set of ENDF-6 thermal scattering files for the
        material described by a particular NCrystal cfg-string.
        """,
        """
        The script uses the endf-parserpy package from IAEA to format and check
        the syntax of the ENDF-6 file:
        """,
        """
        G. Schnabel, D. L. Aldama, R. Capote,
        https://doi.org/10.48550/arXiv.2312.08249
        """,
        f"""
        Note that while the handling of multiple temperatures in one ENDF-6
        file is supported via the {longopt_othertemps} keyword, it is not
        recommended. This is because NCrystal computes an optimal (alpha, beta)
        grid for each material and temperature, while the ENDF format imposes
        the same grid on all temperatures.
        """,
        ]

    # NOTICE ^^^^^^^^^^^
    #
    # When updating ncmat2endf there are 2 main doc texts that have to be
    # checked for updates:
    #
    #  1) The ncmat2endf function doc-string in ncmat2endf.py
    #  2) The ncmat2endf CLI --help text in _cli_ncmat2endf.py
    #

    descr = '\n\n'.join(textwrap.fill(' '.join(e.strip().split()),descrw)
                        for e in descr_sections)

    def exquote(e):
        #prefer " for quoting (for Windows compatibility)
        e = shlex.quote(e)
        if '"' not in e and "'" in e:
            e = e.replace("'",'"')
        if '"' not in e and ',' in e:
            #Add some quotes that shlex did not think necessary:
            e = '"%s"'%e
        return e

    descr += "\n\nExample invocations:\n\n"
    exw = descrw - len(progname) - 7
    expre=f'    $> {progname} '
    for example in examples:
        s=['']
        for e in example:
            e = exquote(e)
            if len(s[-1]+e) > exw:
                s[-1] += ' \\'
                s.append('')
            s[-1] += ' %s'%e
        descr += '%s%s\n'%(expre,s[0])
        for e in s[1:]:
            descr += '%s%s\n'%(' '*len(expre),e)
        descr += '\n'

    usagestr = (
        f'{progname} CFGSTR [{longopt_elastic} {metavar_elastic}]'
        + f' [{longopt_matname} {metavar_matname}]'
        + f' [{longopt_metadata} {metavar_metadata}]\n'
        + (' '*(len(progname)+8))
        + '[<<additional options described below>>]'
    )

    parser = create_ArgumentParser( prog = progname,
                                    description=descr.strip()+'\n',
                                    usage=usagestr,
                                    formatter_class=RawTextHelpFormatter )
    def wrap(t):
        return textwrap.fill(t,width=helpw)

    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('CFGSTR',
                               help=wrap('NCrystal cfg-string defining the'
                                         ' material.'))

    ba = parser.add_argument_group('Commonly used arguments')
    ba.add_argument( '-n', longopt_matname, metavar=metavar_matname,
                     help=wrap('Name of the material to be processed.'
                               'If set ENDF files will be named '
                               f'tsl_element_in_<{metavar_matname}>.endf.'))
    elasmode_default = 'scaled'
    assert elasmode_default in available_elastic_modes
    elasmode_other = list(e for e in available_elastic_modes
                          if e != elasmode_default )
    assert len(elasmode_other)==2
    ba.add_argument('-e', longopt_elastic,metavar=metavar_elastic,
                    help=wrap('Approximation used for the elastic component'
                              f' (default "{elasmode_default}, other options'
                              f' are "{elasmode_other[0]}" and'
                              f' "{elasmode_other[1]}").'
                              ' See DOI:10.1016/j.nima.2021.166227 for'
                              ' meaning of modes.'),
                    type=str, choices=available_elastic_modes,
                    default=elasmode_default)
    ba.add_argument(longopt_metadata,default={},
                    help=wrap('JSON dictionary containing ENDF-6'
                              f' metadata. Run with {longopt_metadata}=help '
                              'for more information.'))
    ba.add_argument('-m',metavar='KEY:VAL', dest='mdata_kvlist',
                    action='append', nargs='+',
                    help=wrap('Add metadata entries. Run with '
                              f'{longopt_metadata}=help for more info.'))
    ba.add_argument(longopt_datenow,action='store_true',
                    help=wrap('Set metadata fields EDATE, DDATE and RDATE'
                              ' to current date.'))

    parser.add_argument('-v','--verbose', action='count',default=0,
                        help=wrap('Increase verbosity. Specify twice'
                                  ' for additional verbosity.'))
    parser.add_argument('--quiet','-q',default=False,action='store_true',
                        help=wrap('Silence non-error output.'))
    ba.add_argument('-d', '--dir', default = '.', metavar='PATH', dest='outdir',
                    help=wrap('Directory for output files (default: current).'))
    ba.add_argument('-i','--index', default = '',
                    metavar='FILE', dest='jsonindex',
                    help=wrap('Story summary of output in FILE (JSON format).'))
    ba.add_argument('-f', '--force',action='store_true',
                    help=wrap('Overwrite output files if'
                              ' they already exist (danger!)'))

    expert_args = parser.add_argument_group('Advanced expert-only arguments')
    expert_args.add_argument(longopt_othertemps,metavar='TVALS',
                             nargs='+',
                             type=float,
                             help=wrap('Additional temperatures to process. As'
                                       ' noted above this is not normally'
                                       ' recommended, and it is preferred'
                                       ' to invoke the script for each'
                                       ' temperature independently using the'
                                       ' "temp" keyword in the cfg-string.') )
    expert_args.add_argument('--smin',metavar='VALUE',
                             type=float, default=default_smin_value,
                             help=wrap('Minimum value of S(alpha, beta) stored'
                                       f' (default: {default_smin_value})'))
    expert_args.add_argument('--emax',
                             type=float, default=default_emax_value,
                             help=wrap('Maximum neutron energy covered by'
                                       ' the kernels'
                                       f' (default: {default_emax_value:g}eV)')
                             )
    expert_args.add_argument('--asymsab',action='store_true',
                             help=wrap('Store S(a,b) in asymmetric form.'))
    expert_args.add_argument('--totsab',action='store_true',
                             help=wrap('Store S(a,b) branches for positive'
                                       ' and negative beta'))

    if return_parser:
        return parser

    #Avoid annoying CFGSTR-missing error when ppl use --mdata=help:
    is_mdata_help = False
    if f'{longopt_metadata}=help' in arglist:
        is_mdata_help = True
    elif longopt_metadata in arglist and 'help' in arglist:
        if arglist.index(longopt_metadata)+1==arglist.index('help'):
            is_mdata_help=True
    if is_mdata_help:
        arglist = [f'{longopt_metadata}=help','dummy']

    args=parser.parse_args(arglist)
    if args.mdata:
        if args.mdata == 'help':
            print(gen_metadata_doc())
            raise SystemExit(0)
        try:
            args.mdata = json.loads(args.mdata)
        except json.JSONDecodeError:
            parser.error(f'Argument to {longopt_metadata} must be a JSON'
                         ' dictionary of key, value pairs')
        else:
            if not isinstance(args.mdata,dict):
                parser.error(f'Argument to {longopt_metadata} must be a JSON'
                             ' dictionary of key, value pairs')
    assert isinstance(args.mdata,dict)
    for ee in args.mdata_kvlist or []:
        for e in ee:
            kv = list(_.strip() for _ in e.split(':',1))
            if not len(kv)==2 or not kv[0]:
                parser.error(f'Invalid parameter for -m: {repr(e)}')
            args.mdata[kv[0]] = kv[1]
    args.m = None

    #map verbosity to 0...3 needed for Python API:
    if args.quiet:
        if args.verbose:
            parser.error('Inconsistent usage of --quiet and --verbose flags')
    else:
        args.verbose = min( 3, args.verbose+1 )

    if args.jsonindex:
        import pathlib
        args.jsonindex = pathlib.Path(args.jsonindex)
        if args.jsonindex.is_file():
            if not args.force:
                parser.error('File already exists (run with --force to'
                             f' overwrite): {args.jsonindex}')
            args.jsonindex.unlink()
            assert not args.jsonindex.is_file()
        if not args.jsonindex.parent.is_dir():
            parser.error('Directory does not exist: {args.jsonindex.parent}')
        args.jsonindex = args.jsonindex.absolute()

    return args

def create_argparser_for_sphinx( progname ):
    return _parseArgs(progname,[],return_parser=True)

@cli_entry_point
def main( progname, arglist ):
    args = _parseArgs( progname, arglist )
    if args.quiet:
        from ._common import ( modify_ncrystal_print_fct_ctxmgr,
                               WarningSpy )
        with modify_ncrystal_print_fct_ctxmgr('block'):
            with WarningSpy( block = True ):
                _main_impl(args)
    else:
        _main_impl(args)

def _main_impl( args ):
    from .ncmat2endf import EndfMetaData, ncmat2endf
    metadata = EndfMetaData()
    if args.mdata:
        metadata.update_from_dict(args.mdata)
    if args.now:
        metadata.set_all_dates_as_now()
    lasym = 0
    if args.totsab:
        lasym = 1
    if args.asymsab:
        lasym += 2
    r = ncmat2endf( args.CFGSTR,
                    material_name = args.name,
                    endf_metadata = metadata,
                    othertemps = args.othertemps,
                    elastic_mode = args.elas,
                    force = args.force,
                    smin = args.smin,
                    emax = args.emax,
                    lasym = lasym,
                    verbosity = args.verbose,
                    outdir = args.outdir )
    if args.jsonindex:
        import json
        r_json = json.dumps( r, indent = 4 ).rstrip() + '\n'
        print(f'Writing index file: {args.jsonindex.name}')
        args.jsonindex.write_text( r_json )

def gen_metadata_doc():
    from ._ncmat2endf_impl import _impl_get_metadata_params_and_docs
    import textwrap

    d = _impl_get_metadata_params_and_docs()
    assert 'LIBNAME' in d
    assert 'ALAB' in d
    txt = ''
    w = 80
    def section( x ):
        return textwrap.fill(' '.join(x.strip().split()),w)

    txt += section(
        f"""Meta-data for ENDF can be provided by the {longopt_metadata}
        option, by specifying a JSON dictionary like:"""
    )
    txt+=('''\n\n  %s='{ "LIBNAME" : "MySuperLib"'''%longopt_metadata
          +''', "ALAB" : "MySuperLab" }'\n\n''')
    txt += section(
        """Or by adding individual items with the -m option like:"""
    )
    txt+=('''\n\n  -m LIBNAME:MySuperLib -m AUTH:"J. Chadwick"\n\n''')
    txt += section('The list of supported meta-data'
                   ' keys and their meaning is:')
    txt += '\n\n'
    kmax = max(len(k) for k in d)
    for k, v in d.items():
        s = f'   {k.rjust(kmax)} : '
        for i,e in enumerate(textwrap.fill( v, width=w-len(s) ).splitlines()):
            txt += (' '*len(s) if i else s) + e + '\n'
    return txt
