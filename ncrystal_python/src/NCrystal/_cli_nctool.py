
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

import os
import math
from ._cliimpl import ( create_ArgumentParser,
                        cli_entry_point,
                        print, warn )
from . import _common as nccommon

@cli_entry_point
def main( progname, arglist ):
    if '--bench' in arglist:
        benchmark_mode(progname,arglist)
    elif '--pytrace' in arglist:
        arglist.remove('--pytrace')
        import trace
        trace.Trace().runfunc(std_main,progname,arglist)
    else:
        std_main( progname, arglist )

#Sphinx doc function. Signature always the following:
def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)

#in unit tests we dont display interactive images and we reduce cpu consumption
#by watching for special env var:
def _is_unittest():
    return nccommon.ncgetenv_bool('TOOL_UNITTESTS')

#Function for importing required python modules which may be missing, to provide
#a somewhat more helpful error to the user:
def import_optpymod(name):
    import importlib
    errmsg = None
    try:
        themod = importlib.import_module(name)
    except ImportError as exc:
        errmsg = 'ERROR: Could not import a required python module: %s'%name
        if maybeThisIsConda() and name in ('matplotlib','numpy'):
            errmsg += (' (looks like you are using conda so you might solve it'
                       ' by running "conda install [-c conda-forge] %s")'%name)
        raise SystemExit(errmsg) from exc
    return themod

def parseArgs( progname, arglist, *, return_parser = False ):
    import argparse
    from . import core as nccore

    if '--mc' in arglist and ( '--help' in arglist or '-h' in arglist ):

        from ._common import print
        print(f"""
usage: {progname} --mc SRCCFG GEOMCFG MATCFG

Invoke the embedded stepping Monte Carlo to display a 4pi diffraction pattern of
neutrons going through a simple sample. SRCCFG describes the source of neutrons,
GEOMCFG the shape of the sample volume, and MATCFG is a normal NCrystal
cfg-string describing the sample material.

In the simplest case, simply supply a neutron energy value (e.g. "20meV",
"1.8Aa", ...)  as SRCCFG and a geometry thickness as GEOMCFG (e.g. "1mm",
"2.5cm", "1mfp", ...). This will result in a pencil beam of neutrons of that
energy impinging centrally on a sphere with that thickness (diameter). The
thickness unit "mfp" means mean-free-path between scattering interactions. The
number of neutrons simulated will be automatically determined, so the simulation
can finish in less than a second.

For advanced users, it is also possibly to provide more options in both SRCCFG
and GEOMCFG, but that is currently considered experimental and not documented
further here.

As an example, to get the diffraction patterm of a a 2Aa neutron through a
sphere of T=250K aluminium with diameter 2cm:

$> {progname} --mc "2Aa" "2cm" "Al_sg225.ncmat;temp=250K"

You can use --logy/--liny to override the y-axis logarithmic setting of the
plot, and if you supply --pdf, a pdf file will produced instead of an
interactive plot being launched.
""")
        return None


    descr="""

The most common usage of this tool is to load input data (usually .ncmat files)
with NCrystal (v%s) and plot resulting isotropic cross sections for thermal
neutrons. This is done by specifying one or more configurations ("cfg-strings"),
which indicates data names (e.g. file names) and optionally cfg parameters
(e.g. temperatures). Specifying more than one configuration, results in a single
comparison plot of the total scattering cross section based on the different
materials. Specifying just a single file, results in a more detailed cross
section plot as well as a 2D plot of generated scatter angles. Other behaviours
can be obtained by specifying flags as indicated below.

"""%nccore.get_version()

    descr=descr.strip()

    epilog="""
examples:
  $ %(prog)s Al_sg225.ncmat
  plot aluminium cross sections and scatter-angles versus neutron wavelength.
  $ %(prog)s Al_sg225.ncmat Ge_sg227.ncmat --common temp=200
  cross sections for aluminium and germanium at T=200K
  $ %(prog)s "Al_sg225.ncmat;dcutoff=0.1" "Al_sg225.ncmat;dcutoff=0.4" "Al_sg225.ncmat;dcutoff=0.8"
  effect of d-spacing cut-off on aluminium cross sections
  $ %(prog)s "Al_sg225.ncmat;temp=20" "Al_sg225.ncmat;temp=293.15" "Al_sg225.ncmat;temp=600"
  effect of temperature on aluminium cross sections
  $ %(prog)s "phases<0.65*Al_sg225.ncmat&0.35*MgO_sg225_Periclase.ncmat>;temp=100K"
  investigate multiphase material at 100K"""

    parser = create_ArgumentParser(prog=progname,
                                   description=descr,
                                   epilog=epilog,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_cfgs', metavar='CFGSTR', type=str, nargs='*',
                        help="""Input data (cfg-strings) to investigate. This
                        can just be simple file-names or full-blown cfg-strings
                        in the usual NCrystal syntax (see also examples
                        below).""")
    parser.add_argument('--version', action='version', version=str(nccore.get_version()))
    parser.add_argument('-d','--dump', action='count',default=0,
                        help=('Dump derived information rather than displaying'
                              ' plots. Specify multiple times to increase'
                              ' verbosity.'))
    parser.add_argument('--mc', nargs=2,metavar=('SRCCFG','GEOMCFG'),
                        help=('Run embedded Monte Carlo app to produce a'
                              ' diffraction pattern of material. Run'
                              ' --mc --help for detailed instructions.'))
    parser.add_argument('--common','-c', metavar='CFG', type=str, default=[],
                        help='Common configuration items that will be applied to all input cfg strings',action='append')
    parser.add_argument('--coh_elas','--bragg', action='store_true',
                        help="""Only generate coherent-elastic (Bragg diffraction) component""")
    parser.add_argument('--incoh_elas', action='store_true',
                        help="""Only generate incoherent-elastic component""")
    parser.add_argument('--sans', action='store_true',
                        help="""Only generate SANS component""")
    parser.add_argument('--elastic', action='store_true',
                        help="""Only generate elastic components (including SANS)""")
    parser.add_argument('--inelastic', action='store_true',
                        help="""Only generate inelastic components""")
    parser.add_argument('-a','--absorption', action='store_true',
                        help="""Include absorption in cross section plots""")
    parser.add_argument('--phases', action='store_true',
                        help="""Show cross section breakdown of a single multiphase material by phase rather than physics process""")
    parser.add_argument('-x','--xrange', type=str,nargs='?',
                        help='Override plot range, e.g. "1e-5:1e2" or "0:10"')
    parser.add_argument('--logy', action='store_true',
                        help='Force y-axis to use logarithmic scale.')
    parser.add_argument('--liny', action='store_true',
                        help='Force y-axis to use linear scale.')
    parser.add_argument('-e','--energy', action='store_true',
                        help="""Show plots versus neutron energy rather than wavelength""")
    parser.add_argument('-p','--pdf', action='store_true',
                        help="""Generate PDF file rather than launching an interactive plot viewer.""")
    parser.add_argument('--test', action='store_true',
                        help="""Perform quick validation of NCrystal installation.""")
    parser.add_argument('--doc', action='count',default=0,
                        help="""Print documentation about the available cfg-str variables. Specify twice for more detailed help.""")
    dpi_default=200
    parser.add_argument('--dpi', default=-1,type=int,
                        help="""Change plot resolution. Set to 0 to leave matplotlib defaults alone.
                        (default value is %i, or whatever the NCRYSTAL_DPI env var is set to)."""%dpi_default)
    parser.add_argument('--cfg',action='store_true',
                        help='Print normalised cfg-string and dump meta-data about loaded physics processes.')
    parser.add_argument('--plugins', action='store_true',
                        help='List currently enabled loaded plugins.')
    parser.add_argument('-b','--browse', action='store_true',
                        help='List data available in standard locations (e.g. the files in the current directory or search path)')
    parser.add_argument('--extract', type=str, default=None, metavar="DATANAME",
                        help='''Extract contents of DATANAME (e.g. a file name) using the same lookup mechanism as used for data
                        specified in NCrystal cfg strings. This can therefore also be used to inspect
                        in-memory (or on-demand created) data.''')

    if return_parser:
        return parser

    args = parser.parse_args(arglist)
    del arglist

    if args.logy and args.liny:
        parser.error('Do not specify both --liny and --logy')
    if not args.logy and not args.liny:
        args.logy = 'auto'
    else:
        args.logy = bool(args.logy)
    assert args.logy in (True,False,'auto')
    delattr(args,'liny')

    if args.xrange:
        try:
            _=args.xrange.split(':')
            if not len(_)==2:
                raise ValueError
            _ = ( float(_[0]), float(_[1]) )
            if not ( _[0]>=0.0 and _[1]>_[0]):
                raise ValueError
        except ValueError:
            parser.error(f'Invalid --xrange argument: "{args.xrange}"')
        args.xrange = _

    has_single_cfgstr = args.input_cfgs and len(args.input_cfgs)==1
    if args.cfg and not has_single_cfgstr:
        parser.error('Option --cfg requires exactly one cfg-string to be specified.')

    if args.phases and not has_single_cfgstr:
        parser.error('Option --phase requires exactly one cfg-string to be specified.')

    if args.dump and not has_single_cfgstr:
        parser.error('Option --dump requires exactly one cfg-string to be specified.')

    if args.mc and not has_single_cfgstr:
        parser.error('Option --mc requires exactly one cfg-string to be specified.')

    if args.extract or args.plugins or args.doc or args.browse:
        return args

    if args.dpi==-1:
        args.dpi = nccommon.ncgetenv_int_nonneg('DPI',dpi_default)

    if args.dpi>3000:
        parser.error('Too high DPI value requested.')

    if args.test:
        if any((args.cfg,
                args.input_cfgs,
                args.dump,
                args.mc,
                args.coh_elas,
                args.incoh_elas,
                args.sans,
                args.elastic,
                args.inelastic,
                args.absorption,
                args.pdf,
                args.phases)):
            parser.error('Do not specify other arguments with --test.')

    ncomp_select = sum((1 if _ else 0) for _ in (args.coh_elas,args.incoh_elas,args.sans,args.elastic,args.inelastic))
    if ncomp_select > 1:
        parser.error('Do not specify more than one of: --coh_elas/--bragg, --incoh_elas, --sans, --elastic or --inelastic.')

    if args.coh_elas:
        args.comp = 'coh_elas'
    elif args.incoh_elas:
        args.comp = 'incoh_elas'
    elif args.sans:
        args.comp = 'sans'
    elif args.elastic:
        args.comp = 'elastic'
    elif args.inelastic:
        args.comp = 'inelastic'
    else:
        args.comp = 'all'

    if args.absorption and ncomp_select>0:
        parser.error('Do not specify --absorption with either of: --coh_elas/--bragg, --incoh_elas, --sans, --elastic or --inelastic.')

    if args.dump and (ncomp_select>0 or args.absorption or args.mc):
        parser.error('Do not specify --dump with either of: --mc, --coh_elas/--bragg, --incoh_elas, --sans, --elastic, --inelastic or --absorption.')

    if args.dump and len(args.input_cfgs)>1:
        parser.error('Do not specify more than one input cfg string with --dump [-d].')

    if args.mc and (ncomp_select>0 or args.absorption):
        parser.error('Do not specify --mc with either of: --coh_elas/--bragg, --incoh_elas, --sans, --elastic, --inelastic or --absorption.')

    if args.mc and len(args.input_cfgs)>1:
        parser.error('Do not specify more than one input cfg string with --mc [-d].')

    args.common=';'.join(args.common)
    return args


def std_main( progname, arglist ):
    from . import core as nccore
    args = parseArgs( progname, arglist )
    if args is None:
        return#ended (e.g. after --help)

    if args.doc:
        full = (args.doc >= 2)
        from .cfgstr import generateCfgStrDoc
        print(generateCfgStrDoc( 'txt_full' if full else 'txt_short' ),end='')
        if not full:
            print("NOTE: Condensed output generated. Specify --doc twice for more details.")
        return

    if args.extract:
        s = nccore.createTextData(args.extract).rawData
        if s is None:
            raise SystemExit('Error: unknown file "%s"'%args.extract)
        print(s,end='')
        return

    if args.plugins:
        from .plugins import browsePlugins
        browsePlugins(dump=True)
        return

    if args.browse:
        from .datasrc import browseFiles
        browseFiles(dump=True)
        return

    if args.test:
        from ._testimpl import test
        test()
        return

    if args.cfg or ( len(args.input_cfgs)==1 and not (args.dump or args.mc) ):
        #Dump cfg debug info if requested or running with just 1 file.
        assert len(args.input_cfgs)==1
        origcfg=args.input_cfgs[0]
        if args.common:
            origcfg += f';{args.common}'
        print(f'==> Debugging cfg-string: "{origcfg}"')
        _ = comp2cfgpars(args.comp)
        if _:
            assert args.comp != 'all'
            print(f'==> Adding due to --{args.comp} flag specified: "{_}"')
            origcfg += ';' + _
        from .cfgstr import normaliseCfg
        normcfg = normaliseCfg(origcfg)
        print(f'==> Normalised cfg-string : "{normcfg}"')
        abs_obj=nccore.createAbsorption(normcfg)
        sc_obj = nccore.createScatter(normcfg)
        print( '==> Absorption process (code level objects):')
        proc_a = abs_obj
        proc_a.dump(' '*27)
        print( '==> Scattering process (code level objects):')
        proc_s = sc_obj
        proc_s.dump(' '*27)
        if args.cfg:
            return

    _mpldpi[0] = args.dpi

    cfgs=[Cfg(e,args.common) for e in args.input_cfgs]
    cfgs_normalisedstrings = [c.cfgstr for c in cfgs]
    for cstr in set(cfgs_normalisedstrings):
        if cfgs_normalisedstrings.count(cstr)!=1:
            warn("Configuration specified more than once: \"%s\""%cstr)

    if not cfgs:
        raise SystemExit('Error: nothing selected.'
                         ' Please run with --help for usage instructions.')

    if args.dump:
        assert len(cfgs)==1
        cfgs[0].get_info().dump(verbose = int(args.dump)-1 )
        return

    if args.mc:
        plot_mmc(cfgs[0].cfgstr,
                 srccfg = args.mc[0],
                 geomcfg = args.mc[1],
                 logy = args.logy,
                 do_pdf = args.pdf )
    else:
        plot_xsect( cfgs, comp  = args.comp, absorption = args.absorption, pdf=args.pdf,
                    versus_energy=args.energy, xrange = args.xrange, logy = args.logy,
                    breakdown_by_phases = args.phases )

        if len(cfgs)==1 and not nccommon.ncgetenv_bool('TOOL_NO2DSCATTER'):
            plot_2d_scatangle( cfgs[0], comp = args.comp, pdf=args.pdf, versus_energy=args.energy, xrange = args.xrange )

        if len(cfgs)==1 and not nccommon.ncgetenv_bool('TOOL_NOAUTOMCPLOT'):
            from . import _mmc as nc_mmc
            auto_params = nc_mmc.quick_diffraction_pattern_autoparams( cfgs[0].cfgstr )
            plot_mmc( cfgs[0].cfgstr,
                      auto_params['neutron_energy_str'],
                      auto_params['material_thickness_str'],
                      logy=args.logy,
                      do_pdf=args.pdf )
    if args.pdf:
        _,_,pdf = import_npplt(True)
        import datetime
        try:
            d = pdf.infodict()
        except AttributeError:
            d={}
        d['Title'] = 'Plots made with NCrystal nctool from file%s %s'%('' if len(args.input_cfgs)==1 else 's',
                                                                       ','.join(os.path.basename(f) for f in args.input_cfgs))
        d['Author'] = 'NCrystal %s (via nctool)'%nccore.get_version()
        d['Subject'] = 'NCrystal plots'
        d['Keywords'] = 'NCrystal'
        d['CreationDate'] = datetime.datetime.today()
        d['ModDate'] = datetime.datetime.today()
        pdf.close()
        print("created %s"%_pdffilename)

def create_ekins(npoints,range_override):
    from ._numpy import _np_geomspace
    if range_override:
        if range_override[0]<=0.0:
            range_override = ( range_override[1]*1e-10, range_override[1] )
        return _np_geomspace(*range_override,npoints)
    else:
        return _np_geomspace(1e-5,1e2,npoints)

def create_wavelengths(np,cfgs,npoints,range_override):
    if range_override:
        wls_min,wls_max = range_override
    else:
        bragg_thresholds = [c.braggthreshold() or 0.0 for c in cfgs]
        fallback = 10.0#materials with no bragg threshold
        max_bragg_threshold = ( max(bragg_thresholds) if bragg_thresholds else None ) or fallback
        wls_max = 1.2 * ( float(int(max_bragg_threshold*1.01+1.0)) if not math.isinf(max_bragg_threshold) else fallback )
        wls_min = 1e-4
    from ._numpy import _np_linspace
    return _np_linspace( wls_min, wls_max, npoints )

_mpldpi=[None]
_pdffilename='ncrystal.pdf'
_npplt = None
def import_npplt(pdf=False):
    #TODO: For interactive usage, we must have consistent global matplotlib
    #manager for all of NCrystal, and we must use a context manager to switch to
    #the agg backend when we need pdf plots. We should also only every change
    #rcParams with the matplotlib.pyplot.rc_context context manager.  See also
    #the discussion at: https://github.com/matplotlib/matplotlib/issues/26362

    #Maybe we should also ensure that all of our matplotlib plots simply return
    #standalone Figure() objects?
    #https://matplotlib.org/devdocs/gallery/user_interfaces/web_application_server_sgskip.html
    global _npplt
    if _npplt:
        #pdf par must be same as last call:
        if bool(pdf) != bool(_npplt[2] is not None):
            raise SystemExit('Currently we do not support switching back and '
                             'forth between pdf plotting in the same process '
                             '(see also '
                             'https://github.com/mctools/ncrystal/issues/202)')
        return _npplt
    np = import_optpymod('numpy')
    mpl = import_optpymod('matplotlib')

    if _mpldpi[0]:
        mpl.rcParams['figure.dpi']=_mpldpi[0]

    #ability to quit plot windows with Q:
    if 'keymap.quit' in mpl.rcParams and 'q' not in mpl.rcParams['keymap.quit']:
        mpl.rcParams['keymap.quit'] = tuple(list(mpl.rcParams['keymap.quit'])+['q','Q'])

    if _is_unittest() or pdf:
        mpl.use('agg')
    if pdf:
        try:
            from matplotlib.backends.backend_pdf import PdfPages
        except ImportError:
            PdfPages = None
        if not PdfPages:
            raise SystemExit("ERROR: Your installation of matplotlib does"
                             " not have the required support for PDF output.")
    plt = import_optpymod('matplotlib.pyplot')

    _npplt = (np,plt,PdfPages(_pdffilename) if pdf else None)

    return _npplt

#functions for creating labels and title:
def _remove_common_keyvals(dicts):
    """remove any key from the passed dicts which exists with identical value in all
    the dicts. Returns a single dictionary with entries thus removed."""
    sets=[set((k,v) for k,v in list(d.items())) for d in dicts]
    common = dict(set.intersection(*sets)) if sets else {}
    for k in list(common.keys()):
        for d in dicts:
            d.pop(k,None)
    return common

def _serialise_cfg(part):
    #transform to list of tuples [(key,value),..] where entries can be
    #(parname,value) or compname/filename.

    mpstart = 'phases<'
    if part._cfg.cfgstr.startswith(mpstart):
        assert part._cfg.cfgstr.count('>')==1
        _=part._cfg.cfgstr[7:].split('>')
        main,trailing_common_cfg = _ if len(_)>1 else (_[0],'')
    else:
        _=part._cfg.cfgstr.split(';',1)
        main,trailing_common_cfg = _ if len(_)>1 else (_[0],'')
    ll = [ ('[FILENAME]', main.strip() ), #using '[]' in special keys avoids clashes
           ('[COMPNAME]', part._compname ) ] #(cfg strs can't contain such chars)
    for e in trailing_common_cfg.split(';'):
        e=e.strip()
        if e == 'ignorefilecfg':
            raise SystemExit('ERROR: The ignorefilecfg keyword is no longer supported')
        elif e:
            k,v=e.split('=')
            ll.append( (k.strip(),v.strip()) )
    return ll

def _classify_differences(parts):
    ll=[]
    for p in parts:
        ll.append( dict(_serialise_cfg(p)) )
    common = _remove_common_keyvals(ll)
    return common,ll

def _cfgdict_to_str(cfgdict):
    fn = cfgdict.pop('[FILENAME]','')
    if '*' in fn:
        #multiphase, fix up a bit
        _=''
        for phase in fn.split('&'):
            fraction,phcfg = phase.split('*')
            if _:
                _ += ' + '
            multsymb = '\u00D7'
            _ += '%s%s(%s)'%(fraction,multsymb,phcfg)
        fn = '{%s}'%_
    o = [fn] if fn else []
    cn = cfgdict.pop('[COMPNAME]','')
    if cfgdict:
        o += [', '.join('%s=%s'%(k,v) for k,v in sorted(cfgdict.items()))]
    if cn:
        o += [ { 'coh_elas':'Coherent elastic',
                 'incoh_elas':'Incoherent elastic',
                 'elastic':'Elastic',
                 'inelastic':'Inelastic',
                 'sans':'SANS',
                 'absorption':'Absorption',
                 'all':'Total scattering',
                 'scattering+absorption':'Total scattering+Absorption'}[cn] ]
    return ' '.join(b%a for a,b in zip(o,['%s','[%s]','(%s)']))

def create_title_and_labels(parts):
    partscfg_common,partscfg_unique = _classify_differences(parts)
    return _cfgdict_to_str(partscfg_common),[(_cfgdict_to_str(uc) or 'default') for uc in partscfg_unique]

def _end_plot(plt,pdf):
    if _is_unittest():
        with open(os.devnull,'wb') as fh:
            plt.savefig(fh,format='raw')
        plt.close()
    elif pdf:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

def comp2cfgpars(comp):
    assert comp in ('coh_elas','incoh_elas','elastic','inelastic','sans','all')
    return { 'coh_elas' : 'incoh_elas=0;inelas=0;sans=0',
             'incoh_elas' : 'coh_elas=0;inelas=0;sans=0',
             'elastic' : 'inelas=0',
             'inelastic' : 'elas=0',
             'sans' : 'coh_elas=0;incoh_elas=0;inelas=0',
             'all' : '' }[comp]


def _frexp10(x):
    exp = int(math.floor(math.log10(abs(x))))
    return x / 10**exp, exp

def _latex_format(x):
    if len('%f'%x)<6:
        return '%f'%x
    b,e = _frexp10(x)
    e=f'10^{e}'
    if b==1:
        return e
    return f'{b:g}'+r'\cdot'+e

def plot_mmc(cfgstr,srccfg,geomcfg,logy,do_pdf):
    assert logy in (True,False,'auto')
    if logy=='auto':
        logy=True
    np,plt,pdf = import_npplt(do_pdf)
    from . import _mmc as nc_mmc
    quick_mode = False
    if ';' not in srccfg and ';' not in geomcfg:
        #NB: this logic is in principle not quite robust in all cases
        #(e.g. srccfg='sphere' would register as quick_mode, and then fail):
        quick_mode = True

    if quick_mode:
        res = nc_mmc.quick_diffraction_pattern(cfgstr,
                                               neutron_energy = srccfg,
                                               material_thickness = geomcfg,
                                               nstat = 'auto' )
        nstat = res.setup_info['nstat']
        title = (f'${_latex_format(nstat)}$ {srccfg} neutrons'
                 f' through {geomcfg} diameter sphere')
    else:
        res = nc_mmc.runsim_diffraction_pattern( cfgstr,
                                                 geomcfg = geomcfg,
                                                 srccfg = srccfg,
                                                 tally_detail_lvl = 2 )
        title = (f'src="{srccfg}", geom="{geomcfg}"')

    res.plot_breakdown( rebin_factor=10,#todo: hardcoded for now
                        logy=logy,
                        title = title,
                        plt = plt,
                        do_show = False
                       )
    _end_plot(plt,pdf)

def plot_xsect(cfgs,comp,absorption,pdf,versus_energy,xrange,logy,breakdown_by_phases):
    assert logy in (True,False,'auto')
    if logy=='auto':
        logy=versus_energy
    assert comp in ('coh_elas','incoh_elas','elastic','inelastic','sans','all')
    assert not absorption or comp=='all'
    scalefactors = None

    if breakdown_by_phases:
        assert len(cfgs)==1
        if cfgs[0].nPhases()==0:
            warn("--phases ignored for a single phase material")
            breakdown_by_phases = False
        else:
            mothercfg = cfgs[0]
            scalefactors = list(mothercfg.getChildPhaseNumberFraction(i) for i in range(mothercfg.nPhases()))
            assert abs(sum(scalefactors)-1.0)<1e-10
            cfgs = list(mothercfg.getChildPhaseCfg(i) for i in range(mothercfg.nPhases()))

    np,plt,pdf = import_npplt(pdf)
    if versus_energy:
        plt.xlabel('Neutron energy [eV]')
    else:
        plt.xlabel('Neutron wavelength [angstrom]')
    plt.ylabel('Cross section [barn/atom]')
    if len(cfgs)==1 and comp in ('all','elastic'):
        if comp=='all':
            parts=[cfgs[0].get_scatter('coh_elas'),cfgs[0].get_scatter('incoh_elas'),
                   cfgs[0].get_scatter('inelastic'),cfgs[0].get_scatter('sans')]
        else:
            assert comp=='elastic'
            parts=[cfgs[0].get_scatter('coh_elas'),cfgs[0].get_scatter('incoh_elas'),cfgs[0].get_scatter('sans')]
        if absorption:
            assert comp=='all'
            parts += [cfgs[0].get_absorption()]
    else:
        if absorption:
            assert comp=='all'
            parts=[c.get_totalxsect() for c in cfgs]
        else:
            parts=[c.get_scatter(comp) for c in cfgs]

    if breakdown_by_phases:
        sc_sans = mothercfg.get_scatter('sans')
        if not sc_sans._nullprocess:
            scalefactors += [1.0]
            parts += [sc_sans]

    if not breakdown_by_phases:
        #trim unused process types (but always show all in case of breakdown_by_phases):
        parts = [p for p in parts if not p._nullprocess]

    if not breakdown_by_phases:
        title,labels = create_title_and_labels(parts)
        if len(set(labels))!=len(labels):
            warn("Comparing identical setups?")
    else:
        title,labels = create_title_and_labels(parts)

    npts = 3000
    if versus_energy:
        ekins = create_ekins(npts,xrange)
    else:
        wavelengths = create_wavelengths(np,cfgs,npts,xrange)
        from .constants import wl2ekin as nc_wl2ekin
        ekins = nc_wl2ekin(wavelengths)
    need_tot = (len(cfgs)==1 and len(parts)>1) or breakdown_by_phases
    xsects_tot = None
    max_len_label = 0

    #colors inspired by http://www.mulinblog.com/a-color-palette-optimized-for-data-visualization/
    col_red = "#F15854"
    partcols = [
        "#5DA5DA", # (blue)
        "#FAA43A", # (orange)
        "#60BD68", # (green)
        #"#B2912F", # (brown)
        "#B276B2", # (purple)
        #"#DECF3F", # (yellow)
        #"#F17CB0", # (pink)
        "#4D4D4D", # (gray)
        ]
    if not need_tot:
        partcols = [col_red]+partcols

    linewidth = 2.0

    xvar = ekins if versus_energy else wavelengths

    ydatarange = {}
    def update_ydatarange(xsects):
        ynz = xsects[np.nonzero(xsects)]
        y0nonzero,y0,y1 = ( ynz.min() if len(ynz) else None), xsects.min(), xsects.max()
        if y0 < ydatarange.get('ymin',float('inf')):
            ydatarange['ymin'] = y0
        if y0nonzero is not None and y0nonzero < ydatarange.get('ymin_nonzero',float('inf')):
            ydatarange['ymin_nonzero'] = y0nonzero
        if y1 > ydatarange.get('ymax',float('-inf')):
            ydatarange['ymax'] = y1

    for ipart,part in enumerate(parts):
        #cfg = part._cfg
        #compname = part._compname
        if part.isOriented():
            raise SystemExit("ERROR: This utility can not produce quick"
                             " cross-section plots for oriented processes"
                             " (but you can still inspect the material with"
                             " --dump)")
        xsects = part.crossSectionIsotropic(ekins)
        if scalefactors is not None:
            xsects *= scalefactors[ipart]
        if need_tot:
            if xsects_tot is None:
                xsects_tot = np.zeros(len(xsects))
            xsects_tot += xsects
        label=labels[ipart]
        max_len_label = max(max_len_label,len(label))
        ls={0:'-',1:'--',2:':'}.get(ipart//len(partcols),'-.')
        plt.plot(xvar,xsects,label=label,color=partcols[ipart%len(partcols)],lw=linewidth,ls=ls)
        update_ydatarange(xsects)
    if need_tot:
        if comp=='all' and len(cfgs)==1 and not breakdown_by_phases:
            #sanity check
            if absorption:
                xsects_tot_direct = cfgs[0].get_totalxsect().crossSectionIsotropic(ekins)
            else:
                xsects_tot_direct = cfgs[0].get_scatter('all').crossSectionIsotropic(ekins)
            xsects_discrepancy = xsects_tot_direct-xsects_tot
            discr_lvl = max(abs(xsects_discrepancy))
            if discr_lvl > 1e-10:
                warn("WARNING: Discrepancy in breakdown into components detected (at the"
                     +f" {discr_lvl} level)!! Some plugins might be incorrectly programmed.")
                plt.plot(xvar,xsects_discrepancy,
                         label='WARNING: Discrepancy!',
                         color="cyan",lw=linewidth)

        update_ydatarange(xsects_tot)
        plt.plot(xvar,xsects_tot,
                 label={'all':'Total','elastic':'Total elastic'}[comp],
                 color="#F15854",lw=linewidth)#red-ish colour (see above)
    leg_fsize = 'large'
    if max_len_label > 40:
        leg_fsize = 'medium'
    if max_len_label > 60:
        leg_fsize = 'small'
    if max_len_label > 80:
        leg_fsize = 'smaller'
    try:
        if len(parts)>1:
            leg=plt.legend(loc='best',fontsize=leg_fsize)
            if hasattr(leg,'set_draggable'):
                leg.set_draggable(True)
            else:
                leg.draggable(True)
    except TypeError:
        plt.legend(loc='best')
    plt.grid()
    single_yval = bool(ydatarange.get('ymin','n/a')==ydatarange.get('ymax','n/a'))
    if single_yval and ydatarange.get('ymin','n/a')==0.0:
        if logy:
            warn('Could not set log scale since curves are 0.0 everywhere')
        plt.gca().set_ylim(0.0,1.0)
    elif logy:
        if ydatarange.get('ymin',1.0) <= 0.0:
            _=ydatarange.get('ymin_nonzero',ydatarange.get('ymax',1.0)*1e-10)
            if _:
                plt.gca().set_ylim(_,None)
        plt.gca().set_yscale('log')
    else:
        plt.gca().set_ylim(0,None)

    if versus_energy:
        plt.gca().set_xlim(ekins[0],ekins[-1])
        plt.gca().set_xscale('log')
    else:
        if wavelengths[0]*100<wavelengths[-1]:
            plt.gca().set_xlim(0.0,wavelengths[-1])
        else:
            plt.gca().set_xlim(wavelengths[0],wavelengths[-1])
    if title:
        if len(title)>30:
            plt.title(title,fontsize='smaller')
        else:
            plt.title(title)
    _end_plot(plt,pdf)

def plot_2d_scatangle(cfg,comp,pdf,versus_energy,xrange):
    assert comp in ('coh_elas','incoh_elas','elastic','inelastic','sans','all')
    part=cfg.get_scatter(comp)

    np,plt,pdf = import_npplt(pdf)

    #reproducible plots:
    import random
    random.seed(123456)

    #higher granularity wavelengths than for 1D plot to avoid artifacts:
    npts = 100 if _is_unittest() else 30000
    if versus_energy:
        ekins = create_ekins(npts,xrange)
    else:
        wavelengths = create_wavelengths(np,[cfg],npts,xrange)
        from .constants import wl2ekin as nc_wl2ekin
        ekins = nc_wl2ekin(wavelengths)

    #get title (label should be uninteresting for a single part):
    title,labels = create_title_and_labels([part])

    #First figure out how many points to put at each wavelength (or energy)
    if not part._nullprocess:
        xsects = part.crossSectionIsotropic(ekins)
        n2d = 100 if _is_unittest() else 25000
        sumxs = np.sum(xsects)
        if sumxs:
            n_at_xvar = np.random.poisson(xsects*n2d/np.sum(xsects))
        else:
            n_at_xvar = np.zeros(len(xsects))
    else:
        n_at_xvar = np.zeros(len(wavelengths))

    xvar = ekins if versus_energy else wavelengths

    n2d=int(np.sum(n_at_xvar))#correction for random fluctuations
    if n2d>0:
        from .constants import wl2ekin as nc_wl2ekin
        plot_angles = np.zeros(n2d)
        plot_xvar = np.zeros(n2d)
        j = 0
        for i,n in enumerate(n_at_xvar):
            i,n = int(i),int(n)
            evalue = xvar[i] if versus_energy else nc_wl2ekin(xvar[i])
            ekinfinal,mu = part.sampleScatterIsotropic(evalue,repeat=int(n))
            plot_angles[j:j+n] = np.arccos(mu)
            plot_xvar[j:j+n].fill(xvar[i])
            j+=n
        plot_angles *= (180./np.pi)
    else:
        plot_angles = None
        plot_xvar = None

    if plot_xvar is not None:
        try:
            plt.scatter(plot_xvar,plot_angles,alpha=0.2,marker='.',edgecolor=None,color='black',s=2,zorder=1)
        except ValueError:
            plt.scatter(plot_xvar,plot_angles,alpha=0.2,edgecolor=None,color='black',s=2,zorder=1)

    if versus_energy:
        plt.gca().set_xlim(ekins[0],ekins[-1])
        plt.gca().set_xscale('log')
    else:
        if wavelengths[0]*100<wavelengths[-1]:
            plt.gca().set_xlim(0.0,wavelengths[-1])
        else:
            plt.gca().set_xlim(wavelengths[0],wavelengths[-1])
    plt.gca().set_ylim(0,180)
    if versus_energy:
        plt.xlabel('Neutron energy [eV]')
    else:
        plt.xlabel('Neutron wavelength [angstrom]')
    plt.ylabel('Scattering angle [degrees]')
    if title:
        plt.title(title)
    plt.grid()
    _end_plot(plt,pdf)

class XSSum:
    #Combine scatter+absorption processes (hence no sampleScatterIsotropic method).
    def __init__(self,*processes):
        self._p = processes[:]
    def crossSectionIsotropic(self,ekin):
        return sum(p.crossSectionIsotropic(ekin) for p in self._p)
    def isOriented(self):
        return any(p.isOriented() for p in self._p)

class Cfg:
    def __init__(self,cfgstr, commoncfgstr):
        from .cfgstr import normaliseCfg
        self._cfgstr = normaliseCfg('%s;%s'%(cfgstr,commoncfgstr))
        self._sc = {}
        self._abs = None
        self._totxs = None
        self._info = None
        self._bt = 'not_init'
        self._iphase = None

    def nPhases(self):
        return len(self.get_info().phases)

    def getChildPhaseCfg(self,iphase):
        assert( iphase < self.nPhases() )
        childcfg = Cfg( self._cfgstr,'phasechoice=%i'%iphase )
        childcfg._iphase = iphase
        return childcfg

    def getChildPhaseNumberFraction(self,iphase):
        #fraction of atoms in phase
        i = self.get_info()
        volfrac,iph = i.phases[iphase]
        return volfrac*iph.numberdensity / i.numberdensity

    def braggthreshold(self):
        """in Aa (or None). Largest BT value for any crystalline phase."""
        if self._bt != 'not_init':
            return self._bt
        def largestbtrecursive(info):
            if info.isSinglePhase():
                return info.braggthreshold
            bts = [ largestbtrecursive(ph[1]) for ph in info.phases ]
            bts = [ e for e in bts if e ]
            return max(bts) if bts else None
        self._bt = largestbtrecursive( self.get_info() )
        return self._bt

    def get_scatter(self,comp = 'all', allowfail = False):
        if comp not in self._sc:
            from . import core as nccore
            extra_cfg = comp2cfgpars(comp)
            cstr = ';'.join([self._cfgstr,extra_cfg])
            try:
                sc = nccore.createScatter(cstr)
            except nccore.NCException:
                if allowfail:
                    return None
                else:
                    raise
            sc._nullprocess = sc.isNull()
            sc._cfg = self
            sc._compname = comp
            self._sc[comp] = sc
        return self._sc[comp]

    def get_absorption(self):
        if not self._abs:
            from . import core as nccore
            a = nccore.createAbsorption(self._cfgstr)
            a._nullprocess = a.isNull()
            a._cfg = self
            a._compname = 'absorption'
            self._abs = a
        return self._abs

    def get_totalxsect(self):
        if not self._totxs:
            a,s = self.get_absorption(),self.get_scatter('all')
            t = XSSum(a,s)
            t._nullprocess = a._nullprocess and s._nullprocess
            t._cfg = self
            t._compname = 'scattering+absorption'
            self._totxs = t
        return self._totxs

    def get_info(self):
        if not self._info:
            from . import core as nccore
            self._info = nccore.createInfo(self._cfgstr)
        return self._info

    @property
    def cfgstr(self):
        return self._cfgstr

def benchmark_mode( progname, arglist):

    parser = create_ArgumentParser(prog=progname,
                                   description='nctool benchmark mode'
                                   ' (triggered by presence of --bench). '
                                   'Measured timing is printed (in '
                                   'milliseconds).')
    parser.add_argument('--bench', action='store_true',
                        help='Flag triggering the benchmark mode.')
    parser.add_argument('cfgstr', metavar='CFGSTR', type=str,nargs=1,
                        help="NCrystal material cfg-string to investigate.")
    parser.add_argument('--onlyinfo', action='store_true',
                        help='Unless this flag is specified, both createInfo()'
                        ' and createScatter() are benchmarked.')
    parser.add_argument('-n','--nrepeat', default=1,type=int,
                        help="""Number of times to measure.""")
    parser.add_argument('--threads', default=1,type=int,
                        help="Number of threads in NCrystal's"
                        " factory thread pool.")
    args=parser.parse_args( arglist )
    nccommon.ncsetenv('FACTORY_THREADS',args.threads)
    from . import misc as nc_misc
    assert len(args.cfgstr)==1
    dt = nc_misc._benchloadcfg( args.cfgstr[0],
                                do_scatter = not args.onlyinfo,
                                nrepeat=args.nrepeat )
    if _is_unittest():
        dt = 0.0
    print(f'{dt*1000.0:.2f}ms')

def maybeThisIsConda():
    import os
    import sys
    return ( os.environ.get('CONDA_PREFIX',None) or
             os.path.exists(os.path.join(sys.base_prefix, 'conda-meta')) )
