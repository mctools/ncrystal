
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

#Fixme unit test examples in mmcscenario.py test!
#Examples, in form of (description, cfgstr, scenariocfg):

_scenariocfg_examples = [
    (
        """Pencil beam of neutrons impinging centrally on a diameter=1mfp (mean
        free path between scatterings) sphere of aluminium. Beam energy is
        chosen to be hopefully interesting for the material, based on Bragg
        threshold and temperature.""",
        "Al_sg225.ncmat",
        ""
    ),
    (
        """Pencil beam of 2Aa neutrons impinging centrally on a diameter=2mm
        sphere of 300K aluminium.""",
        "Al_sg225.ncmat;temp=300K",
        "2Aa pencil on 2mm sphere"
    ),
    (
        """A zero-divergence beam of 10meV neutrons uniformly illuminating a
        diameter=2mm sphere of 80K beryllium.""",
        "Be_sg194.ncmat;temp=80K",
        "10meV on 2mm sphere"
    ),
    (
        """1.8Aa neutrons impinging at right incidence on an infinite slab of
        thickness 10cm filled with humid air.""",
        "gasmix::air/0.9relhumidity",
        "1.8Aa on 10cm slab"
    ),
    (
        """Neutrons at a wavelength which is 99% of the Bragg threshold of PG,
        uniformly illuminating a PG filled sphere whose diameter is 2 times the
        mean free path length between scatterings.""",
        "C_sg194_pyrolytic_graphite.ncmat",
        "0.99BT on 2mfp"
    ),
    (
        """Neutrons at a wavelength which is 80% of the Bragg threshold of copper,
        uniformly illuminating a warm copper sphere whose diameter is 1 times the
        mean free path length between scatterings.""",
        "Cu_sg225.ncmat;temp=400K",
        ""
    ),
]

def _parseArgs( progname, arglist, return_parser=False ):
    if progname=='sb_nccmd_minimc':
        progname='ncrystal_minimc'#fixme if needed for unit test?

    from . import _mmc as ncmmc
    from argparse import RawTextHelpFormatter
    import textwrap
    import shlex
    import pathlib

    tallylists = ncmmc.available_tallies()

    #--full-help keyword needs special workarounds:
    is_fullhelp = any( a.startswith('--f') for a in arglist )
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
    for i,(descr, cfgstr, scenario) in enumerate(_scenariocfg_examples):
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
    to investigate neutron scattering in various geometries including effects of
    multiple scattering."""


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

Run a simulation with NCrystal's builtin MiniMC framework, in order to to
investigate bla...

    results are collected in ...

{descr_examples}

Available energy units: Aa, meV, eV, mfp (mfp is "mean-free-path")
Available length units: mm, cm, m
Available angular units: mm, cm, m

In case of a material cfg-string defining an anisotropic material (e.g. single
crystal), it is important to note that default beam travels along (0,0,1).

  2. Circular uniform beam of 2Aa neutrons illuminating entire 2mm sphere:
  $> {progname} "Al_sg225.ncmat;temp=300K "2Aa on R=2mm sphere"
  3. Pencil beam of 2Aa neutrons centrally on a 2mm sphere:
  $> {progname} "Al_sg225.ncmat;temp=300K "2Aa on D=2mm sphere"

  * 25meV neutrons

  $> {progname} "2Aa,2mm,0.2deg"


Alternatively save or load results of previous runs (fixme).

FIXME: We need to document the available tallies, available via the enginecfg
        tally keyword:
        {tallyflags_hists}
        {tallyflags_other}

"""

    #fixme: can run with just a single (cfgstr) argument, which then defaults to
    #a quick arg of something. I guess the default is a 0.1deg divergent beam
    #illuminating entire sphere uniformly?

    parser = create_ArgumentParser(prog = progname,
                                   description=descr,
                                   formatter_class=RawTextHelpFormatter)
    #bg = parser.add_argument_group('basic options')
    bg = parser #use default group, since '-h/--help' will go there

    bg.add_argument('CFGSTR',type=str, nargs='?',
                        help=hwrap("NCrystal cfg-string defining"
                                  " sample material."))
    bg.add_argument('SCENARIO',type=str,nargs='?',
                    help=hwrap("String defining the simulation setup. For "
                               "more detailed control, use the --geomcfg,"
                               " --srccfg, and --enginecfg flags instead."))

    bg.add_argument('--full-help',action='store_true',
                        help=hwrap("Provide more complete usage instructions."))

    bg.add_argument("--geomcfg",'-g',default=None,type=str,metavar='STR',
                        help=hwrap("String defining sample geometry."))
    bg.add_argument("--srccfg",'-s',default=None,type=str,metavar='STR',
                        help=hwrap("String defining options for the"
                                   " neutron source"))
    bg.add_argument("--enginecfg",'-e',default=None,type=str,metavar='STR',
                        help=hwrap("String defining options for the simulation"
                                   " engine."))

    if not is_basic_help:
        ag = parser.add_argument_group('advanced options')
        ag.add_argument("--decode",action='store_true', help=hwrap("""Show
                        resulting cfgstr, geomcfg, srccfg, and enginecfg."""))

        ag.add_argument("--tally",'-t',metavar='TALLIES',type=str,default='',
                        help=hwrap("""Comma separated list of quantities to
                        tally in the simulations (use --tally==help for a
                        list). Using this is an alternative to selecting the
                        tallies via the "tally" parameter in the enginecfg. If
                        in --input mode and not performing a new simulation, the
                        list is instead used as a filter for which tallies from
                        the file to investigate."""))

        ag.add_argument("--plot",'-p',action='store_true',
                        help=hwrap("""Launch interactive matplotlib plots of
                        results (default if no other action specified)."""))

        ag.add_argument("--dump",'-d',action='store_true',
                        help=hwrap("""Provide summary of results to stdout."""))
        ag.add_argument("--output",'-o',type=str, dest='outputfile',
                        help=hwrap("""Output file name. Must end with .json.gz
                        or .json. As a special case it can be "stdout" if no
                        other output is requested."""))

        #FIXME: Todo: mode which can merge --merge (N args, last one is the new
        #file) results. We should also save the seed in each case (actually it
        #should simply be an enginecfg option!). And finally we should have a
        #--doc mode.
        ag.add_argument("--input",'-i',type=str,dest='inputfile',
                        help=hwrap("""Input file name. Used to load previous
                        results (from --output) instead of running new
                        simulations. In this case, no SCENARIO, --geomcfg,
                        --srccfg, or --enginecfg should be specified."""))

        ag.add_argument("--ref",'-r',type=str,dest='reffile',
                        help=hwrap("""Reference file name. Used with --plot to
                        overlay previous reference results."""))

        ag.add_argument('--force','-f',action='store_true',
                        help=hwrap("""Use with --output to allow overwriting an
                        existing output file."""))
    if return_parser:
        return parser
    args=parser.parse_args(arglist)

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
        p = p.absolute()
        if not p.parent or not p.parent.is_dir():
            parser.error(f'Output directory not found: {p.parent}')
        stem, sufs = p.stem, p.suffixes
        if ( not stem or stem.startswith('.') or stem.startswith('-')
             or ( sufs[-2:]!=['.json','.gz'] and sufs[-1:] != ['.json'] ) ):
             parser.error('Invalid filename (must end with'
                          f' .json or .json.gz): {p.name}')
        return p.resolve()

    args.inputfile = check_infile( args.inputfile )
    args.reffile = check_infile( args.reffile )
    if args.outputfile=='stdout':
        if args.dump or args.decode:
            parser.error('--output=stdout is incompatible with --decode/--dump')
    else:
        args.outputfile = check_outfile( args.outputfile )

    n_geomsrc = sum(int(e is not None) for e in (args.geomcfg,args.srccfg))
    n_geomsrcengine = int(args.enginecfg is not None) + n_geomsrc

    is_mode_stdcfg = bool(n_geomsrcengine>0)
    is_mode_scenario = args.SCENARIO is not None
    is_mode_input = args.inputfile is not None

    if is_mode_input:
        if is_mode_scenario:
            parser.error('Do not specify SCENARIO string when using --input')
        if is_mode_stdcfg:
            parser.error('Do not specify --geomcfg, --srccfg, or --enginecfg'
                         ' when using --input')
        if args.CFGSTR:
            parser.error('Do not specify CFGSTR when using --input')

    is_mode_doc = False#fixme!!! args.doc is not None
    n_modes = sum([is_mode_stdcfg,is_mode_scenario,is_mode_input,is_mode_doc])
    if n_modes==0 and args.CFGSTR is not None:
        #a lone CFGSTR => scenarion mode:
        args.SCENARIO=''
        n_modes, is_mode_scenario = 1, True
    if n_modes!=1:
        parser.error(('%s arguments. Use -h, --help or --full-help for'
                      ' information about proper'
                      ' usage.')%( 'Inconsistent' if n_modes>1 else 'Missing'))
    if ( is_mode_stdcfg or is_mode_scenario ) and args.CFGSTR is None:
        parser.error('Missing CFGSTR argument')
    if is_mode_input and args.CFGSTR is not None:
        parser.error('Do not supply CFGSTR argument with --input')
    if is_mode_doc and args.CFGSTR is not None:
        parser.error('Do not supply CFGSTR argument with --doc')
    if is_mode_scenario and is_mode_stdcfg:
        parser.error('Do not supply any of --srccfg/--geomcfg/--enginecfg if'
                     ' also supplying a SCENARIO string')

    if n_geomsrcengine>0 and n_geomsrc<2:
        parser.error('Options --srccfg and --geomcfg must always be supplied'
                     ' together (--enginecfg is optional and will default to'
                     ' an empty string).')
    if n_geomsrcengine and args.enginecfg is None:
        args.enginecfg = ''

    if is_mode_doc:
        for v in ('decode','plot','dump'):
            if getattr(args,v):
                parser.error('Option --doc and '
                             f'--{v} can not be used together.')

    if args.reffile and not args.plot:
        parser.error('Option --ref requires --plot.')

    if not any( e for e in (is_mode_doc,args.plot,args.dump,
                            args.decode,args.outputfile) ):
        #Default is to plot if nothing else is requested
        args.plot = True

    if not (args.tally or '').strip():
        args.tally=None
    elif args.tally.strip()=='help':
        h = set(tallylists['ALLHISTS'])
        a = set(tallylists['ALL'])
        #fixme: flags should have a description
        print("Available tally quantities:")
        for e in sorted(h):
            print(f'  {e}')
        print("Other tally flags:")
        for e in sorted(a-h):
            print(f'  {e}')
        raise SystemExit
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
    return _parseArgs(progname,[],return_parser=True)

@cli_entry_point
def main( progname, arglist ):
    args = _parseArgs( progname, arglist )

    from . import _mmc as ncmmc
    if args.SCENARIO is not None:
        _ = ncmmc.minimc_decode_scenario(args.CFGSTR,args.SCENARIO)
        args.geomcfg = _['geomcfg']
        args.srccfg = _['srccfg']
        args.enginecfg = _['enginecfg']
        args.SCENARIO = None

    if args.tally and not args.inputfile:
        args.enginecfg += '%stally=%s'%( ';' if args.enginecfg else '',
                                         ','.join(args.tally) )

    res = None
    if args.inputfile:
        from ._mmc_utils import MMCResults
        res = MMCResults( args.inputfile.read_bytes() )

    #Transfer tallies to enginecfg!


    if args.decode:
        import shlex
        print("Simulation setup:")
        if args.inputfile:
            _c,_g,_s,_e = ( args.CFGSTR, args.geomcfg, args.srccfg, args.enginecfg )
        else:
            _c,_g,_s,_e = ( args.CFGSTR, args.geomcfg, args.srccfg, args.enginecfg )
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
        del _c, _g, _s, _e,_exargs,cmdquoted

    if res is None and any( e is not None
                            for e in (args.plot, args.outputfile, args.dump) ):
        assert args.CFGSTR is not None
        res = ncmmc.run_minimc( cfgstr = args.CFGSTR,
                                geomcfg = args.geomcfg,
                                srccfg = args.srccfg,
                                enginecfg = args.enginecfg )



    if args.plot:
        assert res is not None
        for t in res.tallies:
            #fixme filter in case of input!!!
            t.plot()
            #FIXME: plots are missing a title! We should also get some output,
            #like the provided/src stats.

    if args.dump:
        assert res is not None
        res.dump()

    if args.outputfile:
        assert res is not None
        data = res.to_json()
        print("HEJSA",repr(data))
        if args.outputfile == 'stdout':
            #fixme: test no other output!! (or enable quiet!)
            import sys
            sys.stdout.write(data)
        else:
            data = data.encode('utf8')
            if args.outputfile.name.endswith('.gz'):
                import gzip
                data = gzip.compress( data, mtime=0 )#fixme: test if we need
                                                     #lower compresslevel.
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

