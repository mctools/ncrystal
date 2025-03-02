
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

"""Implementation of the built-in tests."""

__all__ = ['test','test_cmdline','test_all']
from ._common import print as _nc_print

def test( verbose = False ):
    """Quick test that NCrystal works as expected in the current installation."""
    _actualtest( verbose = verbose )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def test_cmdline( verbose = False ):
    """Quick test that NCrystal command-line scripts are present and working in
    the current installation."""
    _actual_test_cmdline( verbose = verbose )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def test_all( verbose = False ):
    """Run both test() and test_cmdline()."""
    _actualtest( verbose = verbose )
    _actual_test_cmdline( verbose = verbose )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def _get_prfct( verbose ):
    if verbose and verbose != 'quiet':
        return lambda *a, **kw : _nc_print('::NCrystalTest::',*a,**kw)
    else:
        return lambda *a, **kw : None

def _actualtest( verbose ):
    prfct = _get_prfct( verbose )
    prfct('starting standard Python-API testing')
    try:
        import numpy as _np
    except ImportError:
        _np = None
    from . import core as NC
    from .constants import wl2ekin
    from .plugins import hasFactory

    def require(b):
        if not b:
            raise RuntimeError('check failed')
    def flteq(a,b,rtol=1.0e-6,atol=1.0e-6):
        return abs(a-b) <= 0.5 * rtol * (abs(a) + abs(b)) + atol
    def require_flteq(a,b):
        if not flteq(a,b):
            raise RuntimeError('check failed (%.16g != %.16g, diff %g)'%(a,b,a-b))
        return True
    require(hasFactory('stdncmat'))
    from . import _common as nc_common
    require( nc_common.prettyFmtValue(0.25) == '1/4' )
    require( nc_common.prettyFmtValue(0.25 + 1e-12) != '1/4' )
    require( nc_common.prettyFmtValue(0.25 + 1e-15) == '1/4' )

    _cfgstr='stdlib::Al_sg225.ncmat;dcutoff=1.4'
    prfct(f'Trying to createInfo("{_cfgstr}")')
    al = NC.createInfo(_cfgstr)
    prfct('Verifying loaded Info object')
    require(al.hasTemperature() and require_flteq(al.getTemperature(),293.15))
    require_flteq(al.getXSectFree(),1.39667)
    require_flteq(al.getXSectAbsorption(),0.231)
    require_flteq(al.getDensity(),2.69864547673)
    require_flteq(al.getNumberDensity(),0.06023238256131625)
    require(al.hasDebyeTemperature())

    require(al.hasStructureInfo())
    si=al.getStructureInfo()
    require( si['spacegroup'] == 225 )
    require_flteq(si['a'],4.04958)
    require_flteq(si['b'],4.04958)
    require_flteq(si['c'],4.04958)
    require( si['alpha'] == 90.0 )
    require( si['beta'] == 90.0 )
    require( si['gamma'] == 90.0 )
    require( si['n_atoms'] == 4 )
    require_flteq(si['volume'],66.4094599932)
    require( al.hasHKLInfo() )
    require( al.nHKL() == 3 )
    require_flteq(al.hklDLower(),1.4)
    require( al.hklDUpper() > 1e36 )
    expected_hkl = { 0  : (1, 1, 1, 8, 2.3380261031049243, 1.7731590373262052),
                     1  : (2, 0, 0, 6, 2.02479, 1.7317882793764163),
                     2  : (2, 2, 0, 12, 1.4317427394787092, 1.5757351418107877) }
    for idx,hkl in enumerate(al.hklList()):
        h,k,l_,mult,dsp,fsq = hkl
        require(idx<len(expected_hkl))
        e = expected_hkl[idx]
        require( list(e)[0:4] == [h,k,l_,mult] )
        require_flteq(dsp, e[4])
        require_flteq(fsq, e[5])

    _cfgstr2='stdlib::Al_sg225.ncmat;dcutoff=1.4;incoh_elas=0;inelas=0'
    prfct(f'Trying to createScatter("{_cfgstr2}")')
    #We do all createScatter... here with independent RNG, for reproducibility
    #and to avoid consuming random numbers from other streams.
    alpc = NC.createScatterIndependentRNG(_cfgstr2)
    prfct('Verifying loaded Scatter object')
    require( alpc.name == 'PowderBragg' )
    require( isinstance(alpc.name,str) )
    require( alpc.refCount() in (1,2) )
    require( type(alpc.refCount()) == int ) # noqa E721
    require( alpc.isNonOriented() )
    #_nc_print(alpc.xsect(wl=4.0))
    require_flteq(1.632435821586171,alpc.crossSectionIsotropic(wl2ekin(4.0)) )
    require_flteq(1.632435821586171,alpc.crossSection(wl2ekin(4.0),(1,0,0)))
    require( alpc.crossSectionIsotropic(wl2ekin(5.0)) == 0.0 )

    require( alpc.rngSupportsStateManipulation() )
    require(alpc.getRNGState()=='a79fd777407ba03b3d9d242b2b2a2e58b067bd44')

    alpc.setRNGState('deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    require(alpc.getRNGState()=='deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    alpc_clone = alpc.clone()
    require(alpc.getRNGState()=='deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    require(alpc_clone.getRNGState()=='e0fd16d42a2aced7706cffa08536d869b067bd44')
    alpc_clone2 = alpc_clone.clone(for_current_thread=True)
    require(alpc_clone2.getRNGState()=='cc762bb1160a0be514300da860f6d160b067bd44')
    alpc_clone3 = alpc_clone.clone(rng_stream_index = 12345 )
    require(alpc_clone3.getRNGState()=='3a20660a10fd581bd7cddef8fc3f32a2b067bd44')

    #Pick Nickel at 1.2 angstrom, to also both vdos + incoherent-elastic + coherent-elastic:
    _cfgstr3='stdlib::Ni_sg225.ncmat;dcutoff=0.6;vdoslux=2'
    _seed=2543577
    prfct(f'Trying to createScatterIndependentRNG("{_cfgstr3}",{_seed})')
    nipc = NC.createScatterIndependentRNG(_cfgstr3,_seed)
    prfct('Verifying loaded Scatter object')
    nipc_testwl = 1.2
    #_nc_print(nipc.xsect(wl=nipc_testwl),nipc.xsect(wl=5.0))
    require_flteq(16.76474410391571,nipc.xsect(wl=nipc_testwl))
    require_flteq(16.76474410391571,nipc.xsect(wl=nipc_testwl,direction=(1,0,0)))
    require_flteq(5.958467463288343,nipc.xsect(wl=5.0))

    require( nipc.name == 'ProcComposition' )

    expected = [ ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, 0.8391056916029316 ),
                 ( 0.03200142524676351, -0.37261037212010517 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, -0.15963879335683306 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.0779984939788784, -0.5293576625127443 ),
                 ( 0.05348552589207497, -0.09540771817962344 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.041255667101120046, -0.21139471030502716 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, -0.3915665520281999 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.05750889239879721, -0.5221309343148964 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.08122761653652728, -0.9893394211150188 ),
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.05809677932650489, -0.9514020394895405 ),
                 ( 0.056808478892590906, 0.3042167239859003 ),
                 ( 0.056808478892590906, 0.7378808571510718 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.08045215149882884, -0.8062011016624717 ),
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.06930589080120417, 0.079019907465779 ),
                 ( 0.04019429812207957, -0.9619814414415885 ),
                 ( 0.08983559328581395, -0.5822087429342399 ) ]

    if _np is None:
        ekin,mu=[],[]
        for i in range(30):
            _ekin,_mu=nipc.sampleScatterIsotropic(wl2ekin(nipc_testwl))
            mu += [_mu]
            ekin += [_ekin]
    else:
        ekin,mu = nipc.sampleScatterIsotropic(wl2ekin(nipc_testwl),repeat=30)

    for i in range(len(ekin)):
        #print ( f'    ( {ekin[i]}, {mu[i]} ),');continue
        require_flteq(ekin[i],expected[i][0])
        require_flteq(mu[i],expected[i][1])

    expected = [ ( 0.056808478892590906, (0.07228896531453344, -0.5190173207165885, 0.8517014302500192) ),
                 ( 0.056808478892590906, (-0.9249112255344181, -0.32220112076758217, -0.20180600252850442) ),
                 ( 0.056808478892590906, (-0.15963879335683306, -0.8486615569734178, 0.5042707778277745) ),
                 ( 0.04922198429225973, (-0.9779858857774598, 0.14099376149839138, 0.1538322672218415) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.7905105193171594, -0.6081672667471253) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.8869759070713323, -0.4504882066969593) ),
                 ( 0.056808478892590906, (0.07228896531453344, -0.39741541395284924, -0.914787021249449) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.9768880366798581, -0.1880309758785167) ),
                 ( 0.02561081364848724, (-0.8847741369311427, -0.465745536939693, 0.015994418980024606) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.539797243436807, 0.16202909009269678) ),
                 ( 0.07443151255169597, (-0.6036845347910699, -0.06282202590029042, -0.7947442201839992) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.10854661864786977, 0.5530389874487663) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.7795115518549294, 0.3238994199452849) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.746175597107444, 0.6618128767069312) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.4247181868490453, 0.8996001033001911) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.5555760611065321, -0.6355189486093415) ),
                 ( 0.05736877062456004, (-0.17262993734116835, -0.6849866797325108, 0.7078079918470932) ),
                 ( 0.056808478892590906, (0.3042167239859003, -0.8706122815482211, -0.3866347631352975) ),
                 ( 0.056808478892590906, (-0.7384733804796917, 0.6322144258925643, -0.23443972789660028) ),
                 ( 0.056808478892590906, (-0.15963879335683306, 0.21525619037302965, -0.9634211063505222) ),
                 ( 0.056808478892590906, (0.41359447569500096, 0.4927058865194684, 0.7656242675514158) ),
                 ( 0.056808478892590906, (0.25796367721315083, 0.48520231047621615, 0.8354839670198411) ),
                 ( 0.056808478892590906, (0.5785005938702705, 0.8104481067271115, -0.09225469740985966) ),
                 ( 0.04320250494907263, (-0.03036895176557113, -0.49547892839001373, 0.8680889115120317) ),
                 ( 0.054287027970592844, (-0.360243154961136, -0.9063878964988544, 0.22064870356299168) ),
                 ( 0.056808478892590906, (0.36219866365374176, -0.8822186430862216, 0.3008361577978114) ),
                 ( 0.056808478892590906, (0.7680722413286334, 0.5975216576265994, -0.23028873347945303) ),
                 ( 0.056808478892590906, (0.32922859149927786, -0.9426419619170849, 0.0550878042084668) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.2489220191768986, -0.9631737706493833) ),
                 ( 0.0670453578395921, (-0.8979763975977056, 0.34669277926553477, 0.2710027788835134) ) ]

    for i in range(30):
        out_ekin,outdir = nipc.sampleScatter(wl2ekin(nipc_testwl),(1.0,0.0,0.0))
        #print ( f'    ( {out_ekin}, {outdir} ),');continue
        require_flteq(out_ekin,expected[i][0])
        require_flteq(outdir[0],expected[i][1][0])
        require_flteq(outdir[1],expected[i][1][1])
        require_flteq(outdir[2],expected[i][1][2])

    _cfgstr4="""stdlib::Ge_sg227.ncmat;dcutoff=0.5;mos=40.0arcsec
                            ;dir1=@crys_hkl:5,1,1@lab:0,0,1
                            ;dir2=@crys_hkl:0,-1,1@lab:0,1,0""".replace(' ','').replace('\n','')
    _seed4 = 3453455
    prfct(f'Trying to createScatterIndependentRNG("{_cfgstr4}",{_seed4})')
    gesc = NC.createScatterIndependentRNG(_cfgstr4,_seed4)
    prfct('Verifying loaded Scatter object')
    require_flteq(591.0263476502018,gesc.crossSection(wl2ekin(1.540),( 0., 1., 1. )))
    require_flteq(1.667600586136298,gesc.crossSection(wl2ekin(1.540),( 1., 1., 0. )))
    prfct('standard Python-API testing done')


class CallInspector:

    @property
    def name( self ):
        return self.__name

    @property
    def the_real_inspected_object( self ):
        return self.__realobj

    def __init__( self, name, *, realobj = None, subfcts = None, fmtcall = None ):
        self.__name = name
        self.__realobj = realobj
        self.__subfcts = {}
        self.__fmtcall = fmtcall or _fmtcall
        for sf in (subfcts or []):
            if isinstance(sf,str):
                self.__subfcts[sf] = None
            else:
                assert len(sf)==2 and isinstance(sf,tuple)
                self.__subfcts[sf[0]] = sf[1]

    def __call__( self, *a, **kwargs ):
        return self.__getattrimpl('__call__')(*a,**kwargs)

    def __getattr__( self, attrname ):
        return self.__getattrimpl(attrname)

    def __getattrimpl( self, attrname ):
        if attrname not in self.__subfcts:
            raise RuntimeError(f'Not allowed to access attribute {attrname} of {self.__name} objects.')
        sf_kwargs = self.__subfcts[attrname] or {}
        thefmtcall = sf_kwargs.get('fmtcall') or self.__fmtcall
        def wrapper( *args, **kwargs ):
            _n = f'{self.__name}.{attrname}' if attrname!='__call__' else self.__name
            _nc_print("CALLING %s"%thefmtcall(_n,args,kwargs))
            res = getattr(self.__realobj,attrname)(*args,**kwargs) if self.__realobj else None
            return CallInspector( name = 'ResultOf[%s(..)]'%_n,
                                  realobj = res, **sf_kwargs ) if sf_kwargs else None
        return wrapper

_fmtvalue_default_ndigits = 3
def _fmtvalue( x, *, ndigits = _fmtvalue_default_ndigits ):
    import numbers
    if isinstance(x,numbers.Real) and not isinstance(x,numbers.Integral):
        fmtstring = f'%.{ndigits}g'
        s=fmtstring%x#only 3 to avoid too many spurious test issues due to FPE irrep.
        return s+'.0' if ( s.isdigit() or s[0]=='-' and s[1:].isdigit() ) else s
    return repr(x)

def _fmtcall(fctname,args=tuple(),kwargs={}):
    import numbers
    def _fmt(a):
        if isinstance(a,numbers.Real):
            return _fmtvalue(a)
        def pruneaddr(s):
            while ' object at 0x' in s:
                _ = s.split(' object at 0x',1)
                while _[1] and _[1][0].isalnum():
                    _[1] = _[1][1:]
                s = ' object at SNIPADDR'.join(_)
            return s
        return 'Object[%s]'%a.name if isinstance(a,CallInspector) else pruneaddr(repr(a))
    ll = [ _fmt(a) for a in args ]
    ll += [ '%s=%s'%(k,_fmt(v)) for k,v in sorted(kwargs.items()) ]
    a=','.join(ll)
    return f'{fctname}({a})'

def _create_pyplot_inspector( pass_calls_to_real_plt ):
    import copy
    if pass_calls_to_real_plt:
        import matplotlib.pyplot as realplt
    else:
        realplt = None
    def shorten( x ):
        if hasattr(x,'shape') and len(x.shape)==2:
            return 'Array(shape=%s,content=%s)'%(x.shape,shorten(x.flatten()))
        if isinstance(x,str) or not hasattr(x,'__len__'):
            return x
        #Fewer digits to guard against annoying test errors:
        maxxval = max(x)
        def _fmtthislistval(val):
            if val==0.0:
                return '0.0'
            if abs(val)<abs(maxxval)*1e-13:
                return 'TINYVAL'
                #return _fmtvalue( val, ndigits = max(1,(_fmtvalue_default_ndigits-2)) )
            elif abs(val)<abs(maxxval)*1e-8:
                return 'SMALLVAL'
                #return _fmtvalue( val, ndigits = max(1,(_fmtvalue_default_ndigits-1)) )
            else:
                return _fmtvalue( val )
        if len(x) <= 10:
            return list(_fmtthislistval(e) for e in x)
        else:
            return list(_fmtthislistval(e) for e in x[0:3])+['...']+list(_fmtthislistval(e) for e in x[-3:])
    def _create_shortening_fmtcall( nargs_to_shorten, kwargs_to_shorten = None ):
        def fmtcall_pltplot( name, args, kwargs ):
            plot_args, plot_kwargs = args, kwargs
            if nargs_to_shorten>0:
                plot_args = [ shorten(e) for e in args[0:nargs_to_shorten] ] + list(args[nargs_to_shorten:])
            if any( e in kwargs for e in (kwargs_to_shorten or [])):
                plot_kwargs = copy.deepcopy(kwargs)
                for xy in ('x','y'):
                    if xy in plot_kwargs:
                        plot_kwargs[xy] = shorten(plot_kwargs[xy])
            return _fmtcall(name, plot_args,plot_kwargs)
        return fmtcall_pltplot

    return CallInspector( name = 'plt',
                          realobj = realplt,
                          subfcts = [ 'title','xlabel','ylabel','semilogx','semilogy',
                                      'legend','grid','show','figure',
                                      'savefig','ylim','xlim','colorbar',
                                      'tight_layout','close',
                                      ( 'plot', dict( fmtcall=_create_shortening_fmtcall(2,('x','y')) ) ),
                                      ( 'pcolormesh',dict( subfcts=['set_clim'],
                                                           fmtcall=_create_shortening_fmtcall(3) ) ),
                                     ] )

def _create_pdfpages_inspector( real_pdfpages ):
    return CallInspector( name = 'PdfPages',
                          realobj = real_pdfpages,
                          subfcts = [ ('__call__',dict(subfcts=['savefig','close'])) ] )

def _run_cmd( cmd ):
    import sys
    import subprocess
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        p = subprocess.Popen( cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE )
        output = p.communicate()[0]
        ok = p.returncode==0
    except OSError:
        ok = False
    sys.stdout.flush()
    sys.stderr.flush()
    return ok, output

def _test_cmdline_script_availablity( prfct ):
    import subprocess
    import shutil
    from .cli import cli_tool_list, cli_tool_lookup
    for t in cli_tool_list():
        c = cli_tool_lookup(t)['shellcmd']
        prfct(f"Testing availability of command: {c}")
        if not shutil.which(c):
            raise RuntimeError(f'Command {c} not found!')
        ev = subprocess.run( [c,'--help'], capture_output = True )
        if ev.returncode != 0:
            raise RuntimeError(f'Command "{c} --help" did not run succesfully!')

def _actual_test_cmdline( verbose ):
    from .cli import cli_tool_lookup
    import shlex
    prfct = _get_prfct( verbose )
    prfct('starting testing of cmd-line utilities')
    _test_cmdline_script_availablity( prfct )
    cmds = ['ncrystal-config --help',
            'ncrystal-config -s',
            'nctool --version',
            'nctool --help',
            'nctool --test',
            'ncrystal_endf2ncmat --help',
            'ncrystal_hfg2ncmat --help',
            'ncrystal_ncmat2cpp --help',
            'ncrystal_ncmat2hkl --help',
            'ncrystal_cif2ncmat --help',
            'ncrystal_vdos2ncmat --help',
            'ncrystal_verifyatompos --help',
            'ncrystal-config --show cmakedir']
    for cmd in cmds:
        cmd = shlex.split(cmd)
        cmd[0] = cli_tool_lookup( cmd[0] )['shellcmd']
        cmd = shlex.join(cmd)
        prfct('Trying to run:',cmd)
        ok, output = _run_cmd(cmd)
        if not ok:
            raise RuntimeError('Command failed: %s'%cmd)
    prfct('testing of cmd-line utilities done')
