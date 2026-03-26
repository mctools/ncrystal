
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

def climod_metadata():
    return dict(
        displaygroup = 'main',
        displayorder = 30,
        descr=("Investigate multiple scattering and geometry effects via "
               "NCrystal's own Monte Carlo simulation framework.")
    )


def parseArgs( progname, arglist, return_parser=False ):
    from . import minimc as ncmmc
    from ._mmc_doc import ( doc_subjects,
                            _scenariocfg_examples )
    from argparse import RawTextHelpFormatter
    import textwrap
    import shlex
    import pathlib

    tallylists = ncmmc.tally_info()['tallylists']
    tallyhistavail_str = ', '.join(tallylists['ALLHISTS'])
    #--full-help keyword needs special workarounds:
    is_fullhelp = any( ( a.startswith('--fu')
                         and '--full-help'.startswith(a) )
                       for a in arglist )
    if is_fullhelp:
        arglist.append('-h')

    #Detect if -h/--help is specified and only basic options should be shown:
    def _bhflag(a):
        if 'h' not in a or len(a)<2 or a[0]!='-':
            return False
        if a in ('--h','--he','--hel','--help'):
            return True
        return a[1:].isalpha()
    is_basic_help = ( any( _bhflag(a) for a in arglist )
                      and not ( is_fullhelp or return_parser ) )

    def quote( s ):
        return shlex.quote(s).replace("'",'"')

    width = 80
    helpw = 58
    def hwrap(t):
        return textwrap.fill(' '.join(t.split()),width=helpw)
    descr_examples = ''
    for i,(descr, cfgstr, scenario,key) in enumerate(_scenariocfg_examples):
        assert descr.endswith('.')
        sstr = ' %s'%quote( scenario ) if scenario else ''
        cmd = '%s %s%s'%( progname, quote( cfgstr ), sstr )
        if is_fullhelp:
            descr = descr[:-1]
            descr_examples += textwrap.fill( f'{i+1}. {descr}:',
                                             width=width,
                                             initial_indent='',
                                             subsequent_indent='   ')
            descr_examples += f'\n  $> {cmd}\n\n'
        else:
            descr_examples += f'\n  {i+1}. $> {cmd}'

    #NOTE: Keep the description below synchronised with doc-string of the
    #hfg2ncmat python API function:
    def wrap( x ):
        return textwrap.fill(' '.join(x.strip().split()),width=width)

    intro="""Run a simulation with NCrystal's builtin MiniMC framework, in order
    to investigate neutron scattering in various geometries, including effects
    of multiple scattering."""


    if not is_fullhelp:
        descr=f"""{wrap(intro)}

Quick examples of how to use:

{descr_examples}

Run with --full-help for full list of options and more detailed instructions.
        """.strip()
    else:
        tallyflags_hists = ''
        tallyflags_other = ''
        tally_hists = set(tallylists['ALLHISTS'])
        tally_defaults = set(tallylists['DEFAULT'])
        for i,h in enumerate(sorted(tally_hists)):
            tallyflags_hists += '%s%s%s'%( ',' if i else '',
                                           h,
                                           '*' if h in tally_defaults else '')
        i=0
        for f in sorted(tallylists['ALL']):
            if f not in tally_hists:
                tallyflags_other += '%s%s'%( ',' if i else '', f )
                i += 1
        descr=f"""

Run a simulation with NCrystal's builtin MiniMC framework, in order to
investigate effects of geometry and multiple scattering with a particular
NCrystal material and neutron source.

When using this command-line interface, neutron quantities are always tallied in
histograms as the neutrons exit the geometry. By default these results are
simply shown as interactive plots, but options below allow for more advanced
usage such as printing or persistification into files.

For more info in general about the MiniMC or the various configuration strings,
visit: https://github.com/mctools/ncrystal/wiki/minimc

Alternatively, access documentation via the --doc option, or simply input the
string "help" into any option taking a string. For instance, both --doc=src and
--srccfg=help will result in the detailed documentation for source configuration
being shown.

Here are a few examples of usage based on an approach with SCENARIO strings:

{descr_examples}

The --input and --output options allow persistification of results to JSON
format.
"""

    parser = create_ArgumentParser(prog = progname,
                                   description=descr,
                                   formatter_class=RawTextHelpFormatter)
    #bg = parser.add_argument_group('basic options')
    bg = parser #use default group, since '-h/--help' will go there

    bg.add_argument('CFGSTR',type=str, nargs='?',
                        help=hwrap("NCrystal cfg-string defining"
                                  " sample material."))
    bg.add_argument('SCENARIO',type=str,nargs='?',
                    help=hwrap("String defining the simulation geometry and"
                               " neutron source in terms of a simple scenario."
                               " For  more detailed control, use the --geomcfg"
                               " and --srccfg flags instead."))

    bg.add_argument('--full-help',action='store_true',
                        help=hwrap("Provide more complete usage instructions."))

    bg.add_argument("--tally",'-t',metavar='TALLIES',type=str,default='',
                    help=hwrap(f"""Comma separated list of tally quantities
                    (select from: {tallyhistavail_str}). Using this option is
                    just a handy alternative to setting the "tally" parameter in
                    the enginecfg string. The exception is if running with
                    --input, in which case --tally can be used to select which
                    tallies from the file to show. Use --tally=help for more
                    information about tallies, and the --enginecfg option for
                    more detailed control of e.g. histogram binnings."""))

    bg.add_argument("--plot",'-p',action='store_true',
                    help=hwrap("""Launch interactive matplotlib plots of
                    results (default if no other action specified)."""))

    if not is_basic_help:
        ag = parser.add_argument_group('advanced options')
        ag.add_argument("--geomcfg",'-g',default=None,type=str,metavar='STR',
                        help=hwrap("""String defining simulation geometry (use
                        --geomcfg=help for more info)."""))
        ag.add_argument("--srccfg",'-s',default=None,type=str,metavar='STR',
                        help=hwrap("""String defining options for the
                         neutron source (use --srccfg=help for more info)."""))
        ag.add_argument("--enginecfg",'-e',default=None,type=str,metavar='STR',
                        help=hwrap("""String defining options for the simulation
                        engine (use --enginecfg=help for more info)."""))
        ag.add_argument("--dump",'-d',action='store_true',
                        help=hwrap("""Provide summary of results to stdout."""))
        ag.add_argument("--output",'-o',type=str, dest='outputfile',
                        help=hwrap("""Name of file in which to store results.
                        Must end with .json.gz or .json. As a special case it
                        can be "stdout" if no other output is requested."""))
        ag.add_argument('--force','-f',action='store_true',
                        help=hwrap("""Use with --output to allow overwriting an
                        existing output file."""))
        ag.add_argument("--input",'-i',type=str,dest='inputfile',
                        help=hwrap("""Input file name. Used to load previous
                        results (from --output) instead of running new
                        simulations. In this case, no SCENARIO, --geomcfg,
                        --srccfg, or --enginecfg should be specified."""))
        ag.add_argument('--quiet','-q',default=False,action='store_true',
                        help=wrap('Silence non-error output (automatic if'
                                  ' --output=stdout).'))
        assert ( set(doc_subjects) == set(['engine','src','geom',
                                           'scenario']) ),"update --doc text"
        ag.add_argument('--doc',choices=doc_subjects,
                        help=hwrap("""Show documentation for the configuration
                        strings available for geometry, source, engine, or
                        scenarios. Alternatively, the same documentation is
                        generated if any of the cfg-strings are equal to the
                        string "help"."""))
        ag.add_argument("--decode",action='store_true', help=hwrap("""Show
                        the final cfgstr, geomcfg, srccfg, and enginecfg."""))

    if return_parser:
        return parser
    args=parser.parse_args(arglist)

    #Accept magic "help" cfgstrings for quick help on any subject:
    if args.CFGSTR=='help':
        parser.error("Use -h/--help for instructions"
                     " (or --full-help for even more)")
    if args.SCENARIO=='help':
        args.doc = 'scenario'
        return args
    elif args.enginecfg=='help':
        args.doc = 'engine'
        return args
    elif args.srccfg=='help':
        args.doc = 'src'
        return args
    elif args.geomcfg=='help':
        args.doc = 'geom'
        return args
    elif args.tally=='help':
        return args
    is_mode_doc = args.doc is not None

    def check_infile( f ):
        if not f:
            return
        p = pathlib.Path(f)
        if not p.is_file():
            parser.error(f'Missing file: {f}')
        return p.absolute().resolve()

    def check_outfile( f ):
        if not f:
            return
        p = pathlib.Path(f)
        if p.is_file() and not args.force:
            parser.error(f'Output file aready exists: {f}')
        if p.is_dir():
            parser.error(f'Output file is a directory: {f}')
        if not p.parent or not p.parent.is_dir():
            parser.error(f'Output directory not found: {p.parent}')
        p = p.absolute()
        stem, sufs = p.stem, p.suffixes
        if ( not stem or stem.startswith('.') or stem.startswith('-')
             or ( sufs[-2:]!=['.json','.gz'] and sufs[-1:] != ['.json'] ) ):
             parser.error('Invalid filename (must end with'
                          f' .json or .json.gz): {p.name}')
        return p.resolve()

    args.inputfile = check_infile( args.inputfile )
    if args.outputfile=='stdout':
        for a in ['dump','decode','doc','plot']:
            if getattr(args,a,False):
                parser.error(f'--output=stdout is incompatible with --{a}')
    else:
        args.outputfile = check_outfile( args.outputfile )

    n_geomsrc = sum(int(e is not None) for e in (args.geomcfg,args.srccfg))
    #n_geomsrcengine = int(args.enginecfg is not None) + n_geomsrc

    is_mode_stdcfg = bool(n_geomsrc>0)
    is_mode_scenario = args.SCENARIO is not None
    is_mode_input = args.inputfile is not None
    is_mode_doc = args.doc is not None

    if args.decode:
        if args.dump:
            parser.error('Incompatible options: --decode and --dump')
        if is_mode_doc:
            parser.error('Incompatible options: --decode and --doc')
        if args.plot:
            parser.error('Incompatible options: --decode and --plot')
    if is_mode_input:
        if args.CFGSTR is not None:
            parser.error('Do not specify CFGSTR when using --input')
        assert not is_mode_scenario#not possible unless CFGSTR is set
        #parser.error('Do not specify SCENARIO string when using --input')
        if is_mode_doc:
            parser.error('Incompatible options: --doc and --input')
        if is_mode_stdcfg:
            parser.error('Do not specify --geomcfg or --srccfg'
                         ' when using --input')
        if args.enginecfg is not None:
            parser.error('Do not specify --enginecfg'
                         ' when using --input (use --tally instead if you'
                         ' are trying to filter tallies).')

    n_modes = sum([is_mode_stdcfg,is_mode_scenario,is_mode_input,is_mode_doc])
    if n_modes==0 and args.CFGSTR is not None:
        #a lone CFGSTR => scenario mode:
        args.SCENARIO=''
        n_modes, is_mode_scenario = 1, True
    if is_mode_scenario and is_mode_stdcfg:
        parser.error('Do not supply --srccfg or --geomcfg if'
                     ' also supplying a SCENARIO string')
    if n_modes!=1:
        parser.error(('%s arguments. Use -h, --help or --full-help for'
                      ' information about proper'
                      ' usage.')%( 'Inconsistent' if n_modes>1 else 'Missing'))

    if ( is_mode_stdcfg or is_mode_scenario ) and args.CFGSTR is None:
        parser.error('Missing CFGSTR argument')
    #already tested previously:
    assert not ( is_mode_input and args.CFGSTR is not None )
    if is_mode_doc and args.CFGSTR is not None:
        parser.error('Do not supply CFGSTR argument with --doc')

    if n_geomsrc == 1:
        parser.error('Options --srccfg and --geomcfg must always be supplied'
                     ' together (--enginecfg is optional and will default to'
                     ' an empty string).')
    if ( is_mode_stdcfg or is_mode_scenario ) and args.enginecfg is None:
        args.enginecfg = ''

    if is_mode_doc and args.outputfile:
        parser.error('Option --doc and --output can not be used together.')

    if not any( e for e in (is_mode_doc,args.plot,args.dump,
                            args.decode,args.outputfile) ):
        #Default is to plot if nothing else is requested
        args.plot = True

    if not (args.tally or '').strip():
        args.tally=None
    else:
        #At this point, merely test that all entries are valid TallyFlags:
        t = set([ e.strip() for e in args.tally.split(',')] )
        a = set(tallylists['ALL'])
        if t-a:
            parser.error('Unsupported tally flag: %s'%((t-a).pop()))
        args.tally=sorted(t)

    if args.tally and args.enginecfg and 'tally' in args.enginecfg:
        if ';tally=' in (';'+''.join(args.enginecfg.split())):
            parser.error('Can not use --tally when the --enginecfg'
                         ' also contains "tally=..."')

    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)


@cli_entry_point
def main( progname, arglist ):
    args = parseArgs( progname, arglist )
    do_quiet = ( args.quiet or args.outputfile == 'stdout' )
    if args.doc:
        do_quiet = False
    if do_quiet:
        from ._common import ( modify_ncrystal_print_fct_ctxmgr,
                               WarningSpy )
        with modify_ncrystal_print_fct_ctxmgr('block'):
            with WarningSpy( block = True ):
                _main_impl(progname,args)
    else:
        _main_impl(progname,args)

def print_tally_help():
    from . import minimc as ncmmc
    tallylists = ncmmc.tally_info()['tallylists']
    h = set(tallylists['ALLHISTS'])
    a = set(tallylists['ALL'])
    hi = ncmmc.tally_info()['tallyhistinfo']
    print("Available tally quantities:")
    maxe = max(len(e) for e in h)
    for e in sorted(h):
        line = f'  {e.rjust(maxe)} : {hi[e]["short_descr"]}'
        unit = hi[e]["unit"]
        if unit:
            line += f' ({unit})'
        print(line)
    print("Other tally flags (affects all histograms if set):")
    for e in sorted(a-h):
        print(f'  {e}')

def _main_impl( progname, args ):
    from . import minimc as ncmmc
    if args.tally == 'help':
        print_tally_help()
        return
    if args.doc is not None:
        ncmmc.gen_doc(args.doc,'print')
        return

    if args.SCENARIO is not None:
        _ = ncmmc.decode_scenario(args.CFGSTR,args.SCENARIO)
        args.geomcfg = _['geomcfg']
        args.srccfg = _['srccfg']

    if args.tally and not args.inputfile:
        args.enginecfg += '%stally=%s'%( ';' if args.enginecfg else '',
                                         ','.join(args.tally) )

    res = None
    if args.inputfile:
        res = ncmmc.MMCResults( args.inputfile.read_bytes() )

    if args.decode:
        import shlex
        print("Simulation setup:")
        if args.inputfile:
            assert res is not None
            _c = res.setup['material']['cfgstr']
            _g = res.setup['geom']['cfgstr']
            _s = res.setup['src']['cfgstr']
            _e = res.setup['engine']['cfgstr']
        else:
            _c,_g,_s,_e = ( args.CFGSTR,
                            args.geomcfg, args.srccfg, args.enginecfg )
        #(re)normalise cfg-strings to be sure to present them as such (and as an
        #added test of inputfile compatiblity):
        _g = ncmmc.decode_cfgstr(_g,'geom')['cfgstr']
        _s = ncmmc.decode_cfgstr(_s,'src')['cfgstr']
        _e = ncmmc.decode_cfgstr(_e,'engine')['cfgstr']
        from .cfgstr import normaliseCfg
        _c = normaliseCfg(_c)
        print('  material cfgstr: "%s"'%_c)
        print('          geomcfg: "%s"'%_g)
        print('           srccfg: "%s"'%_s)
        print('        enginecfg: "%s"'%_e)
        print('Command to reproduce:')
        _exargs = [progname,_c,'-g',_g,'-s',_s,'-e',_e]
        cmdquoted = shlex.join(_exargs)
        #prefer " over ' for windows:
        assert '"' not in ''.join(_exargs)
        assert "'" not in ''.join(_exargs)
        cmdquoted = cmdquoted.replace("'",'"')
        print('  %s'%cmdquoted)
        return

    if res is None and any( e is not None
                            for e in (args.plot, args.outputfile, args.dump) ):
        assert args.CFGSTR is not None
        res = ncmmc.run( cfgstr = args.CFGSTR,
                         geomcfg = args.geomcfg,
                         srccfg = args.srccfg,
                         enginecfg = args.enginecfg )

    #If running with an input file we are using --tally to filter what to show.
    tally_show_filter = None
    if args.tally and args.inputfile:
        _missing = set(args.tally)-set(res.tally_names)
        _select = set(args.tally)-_missing
        if _missing:
            from ._common import warn as ncwarn
            ncwarn('Indicated tallies missing in input:'
                   ' %s'%(' '.join(sorted(_missing))))
        def tally_show_filter( tname ):
            return tname in _select

    if args.plot:
        assert res is not None
        tallies=res.tallies
        if not res.tallies:
            from ._common import warn as ncwarn
            ncwarn('No tallies are available to plot.')
        for t in tallies:
            if ( tally_show_filter is None ) or tally_show_filter(t.name):
                print(f"Plotting tally: {t.name}")
                t.plot()
            else:
                print(f"Skipping tally: {t.name}")

    if args.dump:
        assert res is not None
        res.dump( tally_filter_fct = tally_show_filter )

    if args.outputfile:
        assert res is not None
        data = res.to_json()
        if args.outputfile == 'stdout':
            #Do not use print(..) which is affected by --quiet:
            import sys
            sys.stdout.write(data)
        else:
            data = data.encode('utf8')
            if args.outputfile.name.endswith('.gz'):
                import gzip
                data = gzip.compress( data, mtime=0 )
            else:
                if args.outputfile.is_file():
                    if args.force:
                        args.outputfile.unlink()
                    else:
                        raise SystemExit('ERROR: conflicting output file '
                                         'suddenly  appeared:'
                                         f' {args.outputfile}')
            assert args.outputfile.parent.is_dir()
            args.outputfile.write_bytes(data)
            print("Wrote: %s"%args.outputfile.name)
