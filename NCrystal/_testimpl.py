"""

Implementation of the built-in unit tests.

"""

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

__all__ = ['test','test_cmdline','test_cmake','test_extra','test_all']
from ._common import print as _nc_print
import contextlib as _contextlib

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

def test_cmake( verbose = False ):
    """Quick test that the NCrystal installation supports down-stream CMake/C++
    projects (possibly after using ncrystal-config to get CMAKE_PREFIX_PATH and
    (DY)LD_LIBRARY_PATHS setup).
    """
    _actual_test_cmake( verbose = verbose, ignore_if_absent = False )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def test_extra( verbose = False ):
    """Run both test() and test_cmdline(). If CMake is present, also run test_cmake()."""
    _actualtest( verbose = verbose )
    _actual_test_cmdline( verbose = verbose )
    _actual_test_cmake( verbose = verbose, ignore_if_absent = True )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def test_all( verbose = False ):
    """Run both test(), test_cmdline(), and test_cmake()."""
    _actualtest( verbose = verbose )
    _actual_test_cmdline( verbose = verbose )
    _actual_test_cmake( verbose = verbose, ignore_if_absent = False )
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
        h,k,l,mult,dsp,fsq = hkl
        require(idx<len(expected_hkl))
        e = expected_hkl[idx]
        require( list(e)[0:4] == [h,k,l,mult] )
        require_flteq(dsp, e[4])
        require_flteq(fsq, e[5])

    _cfgstr2='stdlib::Al_sg225.ncmat;dcutoff=1.4;incoh_elas=0;inelas=0'
    prfct(f'Trying to createScatter("{_cfgstr2}")')
    #We do all createScatter... here with independent RNG, for reproducibility
    #and to avoid consuming random numbers from other streams.
    alpc = NC.createScatterIndependentRNG(_cfgstr2)
    prfct('Verifying loaded Scatter object')
    require( alpc.name == 'PCBragg' )
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
    require_flteq(16.76322537767633,nipc.xsect(wl=nipc_testwl))
    require_flteq(16.76322537767633,nipc.xsect(wl=nipc_testwl,direction=(1,0,0)))
    require_flteq(5.958094744249944,nipc.xsect(wl=5.0))

    require( nipc.name == 'ProcComposition' )

    expected = [ ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, 0.8391056916029316 ),
                 ( 0.032001313096074194, -0.3726211530494784 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, -0.15963879335683306 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.07799891342244511, -0.5293689067509328 ),
                 ( 0.05348597239804991, -0.09542895891111686 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.04125596596989546, -0.2114086959439411 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, -0.3915665520281999 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.05750930296126236, -0.5221485562238027 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.0812282226184265, -0.9893430983107131 ),
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.058097184485747494, -0.951408724433637 ),
                 ( 0.056808478892590906, 0.3042167239859003 ),
                 ( 0.056808478892590906, 0.7378808571510718 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.08045270890565197, -0.8062090812584963 ),
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.06930621305459778, 0.0790013056899895 ),
                 ( 0.04019454300337846, -0.9619857378392679 ),
                 ( 0.08983663449895618, -0.5822245509723732 ) ]

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
                 ( 0.0492224669270449, (-0.9779916301402904, 0.140975569549056, 0.15381241876342955) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.7905105193171594, -0.6081672667471253) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.8869759070713323, -0.4504882066969593) ),
                 ( 0.056808478892590906, (0.07228896531453344, -0.39741541395284924, -0.914787021249449) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.9768880366798581, -0.1880309758785167) ),
                 ( 0.025610418737037184, (-0.884783168996826, -0.4657283985847863, 0.015993830422524287) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.539797243436807, 0.16202909009269678) ),
                 ( 0.07443181306716587, (-0.6036941700256581, -0.06282145095069493, -0.7947369466543514) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.10854661864786977, 0.5530389874487663) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.7795115518549294, 0.3238994199452849) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.746175597107444, 0.6618128767069312) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.4247181868490453, 0.8996001033001911) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.5555760611065321, -0.6355189486093415) ),
                 ( 0.05736918247226906, (-0.17265143322057086, -0.684984059616355, 0.7078052844380157) ),
                 ( 0.056808478892590906, (0.3042167239859003, -0.8706122815482211, -0.3866347631352975) ),
                 ( 0.056808478892590906, (-0.7384733804796917, 0.6322144258925643, -0.23443972789660028) ),
                 ( 0.056808478892590906, (-0.15963879335683306, 0.21525619037302965, -0.9634211063505222) ),
                 ( 0.056808478892590906, (0.41359447569500096, 0.4927058865194684, 0.7656242675514158) ),
                 ( 0.056808478892590906, (0.25796367721315083, 0.48520231047621615, 0.8354839670198411) ),
                 ( 0.056808478892590906, (0.5785005938702705, 0.8104481067271115, -0.09225469740985966) ),
                 ( 0.04320288783696474, (-0.030385860971878235, -0.49547867364837517, 0.8680884651996275) ),
                 ( 0.05428746831317616, (-0.3602629693021255, -0.9063804616575692, 0.22064689364463652) ),
                 ( 0.056808478892590906, (0.36219866365374176, -0.8822186430862216, 0.3008361577978114) ),
                 ( 0.056808478892590906, (0.7680722413286334, 0.5975216576265994, -0.23028873347945303) ),
                 ( 0.056808478892590906, (0.32922859149927786, -0.9426419619170849, 0.0550878042084668) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.2489220191768986, -0.9631737706493833) ),
                 ( 0.0670456981700729, (-0.8979842133391931, 0.34668021323231085, 0.2709929562678522) ) ]

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
    require_flteq(591.025731949681,gesc.crossSection(wl2ekin(1.540),( 0., 1., 1. )))
    require_flteq(1.666984885615526,gesc.crossSection(wl2ekin(1.540),( 1., 1., 0. )))
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
    l = [ _fmt(a) for a in args ]
    l += [ '%s=%s'%(k,_fmt(v)) for k,v in sorted(kwargs.items()) ]
    a=','.join(l)
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

def _run_cmd( cmd, env = None ):
    import sys
    import subprocess
    sys.stdout.flush()
    sys.stderr.flush()
    if env and sys.platform == 'darwin':
        #osx's system protection #$%#^ is so damn annoying:
        for name in ['DYLD_LIBRARY_PATH','LD_LIBRARY_PATH']:
            v = env.get(name)
            if v is not None:
                from shlex import quote
                #NB: "export A=B" is posix, not bash, so hopefully OK on all OSX
                #shells:
                cmd = 'export %s=%s && %s'%(name,quote(v),cmd)
    try:
        p = subprocess.Popen( cmd, shell=True, env=env,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE )
        output = p.communicate()[0]
        ok = p.returncode==0
    except OSError:
        ok = False
    sys.stdout.flush()
    sys.stderr.flush()
    return ok, output

def _actual_test_cmdline( verbose ):
    prfct = _get_prfct( verbose )
    prfct('starting testing of cmd-line utilities')
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
        prfct('Trying to run:',cmd)
        ok, output = _run_cmd(cmd)
        if not ok:
            raise RuntimeError('Command failed: %s'%cmd)
    prfct('testing of cmd-line utilities done')

@_contextlib.contextmanager
def _work_in_tmpdir():
    """Context manager for working in a temporary directory (automatically created+cleaned) and then switching back"""
    import os
    import tempfile
    the_cwd = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            yield
    finally:
        os.chdir(the_cwd)

def _actual_test_cmake( verbose = False, ignore_if_absent = False ):
    prfct = _get_prfct( verbose )
    prfct('starting testing of compiled downstream cmake-based projects')
    import shutil
    cmakecmd = shutil.which('cmake')
    if not cmakecmd:
        if ignore_if_absent:
            prfct('Skipping CMake test since "cmake" command not found.')
            return
        raise RuntimeError('cmake command not found')

    prfct('Using cmake command: %s'%cmakecmd)


    _singleline2str = lambda x : ((x.decode() if hasattr(x,'decode') else x) or '').strip()

    _cmd = 'ncrystal-config --show cmakedir'
    ok, ncrystal_cmakedir = _run_cmd(_cmd)
    if not ok:
        raise RuntimeError('Command failed: %s'%_cmd)
    ncrystal_cmakedir = _singleline2str(ncrystal_cmakedir)
    prfct('ncrystal-config gives cmakedir: %s'%ncrystal_cmakedir)
    _cmd = 'ncrystal-config --show libdir'
    ok, ncrystal_libdir = _run_cmd(_cmd)
    if not ok:
        raise RuntimeError('Command failed: %s'%_cmd)
    ncrystal_libdir = _singleline2str(ncrystal_libdir)
    prfct('ncrystal-config gives libdir: %s'%ncrystal_libdir)

    def prepend_to_path( envdict, pathvar, entry ):
        from shlex import quote
        entry = quote(((entry.decode() if hasattr(entry,'decode') else entry) or '').strip())
        if not entry:
            return
        l = [ entry ]
        orig = envdict.get(pathvar,'')
        if orig:
            l.append( orig )
        envdict[pathvar] = ':'.join(l)

    import os
    cmake_env = os.environ.copy()
    runtime_env = os.environ.copy()
    prfct('Testing with cmakedir added to CMAKE_PREFIX_PATH and libdir added to (DY)LD_LIBRARY_PATH')
    prepend_to_path( cmake_env, 'CMAKE_PREFIX_PATH', ncrystal_cmakedir )
    prepend_to_path( runtime_env, 'LD_LIBRARY_PATH', ncrystal_libdir )
    prepend_to_path( runtime_env, 'DYLD_LIBRARY_PATH', ncrystal_libdir )
    def runcmake(args,apply_CMAKE_ARGS=True):
        cmake_args = cmake_env.get('CMAKE_ARGS','') if apply_CMAKE_ARGS else ''
        cmd=f'{cmakecmd} {cmake_args} {args}'
        ok, output = _run_cmd(cmd, env=cmake_env)
        prfct(f'Launching: {cmd}')
        if not ok:
            raise RuntimeError(f'Test command failed: {cmd}')

    downstream_cmake = """
cmake_minimum_required(VERSION 3.10...3.24)
project(MyExampleProject LANGUAGES CXX)
find_package(NCrystal REQUIRED)
file(WRITE "${PROJECT_BINARY_DIR}/cppsrc.cpp" "
#include \\"NCrystal/NCrystal.hh\\"
#include <iostream>
int main() {
  auto pc = NCrystal::createScatter( \\"Al_sg225.ncmat;dcutoff=0.5;temp=25C\\" );
  auto wl = NCrystal::NeutronWavelength{1.8};
  auto xsect = pc.crossSectionIsotropic( wl );
  std::cout << \\"Powder Al x-sect at \\" << wl << \\" is \\" << xsect << std::endl;
  return 0;
}
")
add_executable(exampleapp "${PROJECT_BINARY_DIR}/cppsrc.cpp")
target_link_libraries( exampleapp NCrystal::NCrystal )
install( TARGETS exampleapp DESTINATION bin )
"""

    import pathlib
    with _work_in_tmpdir():
        td = pathlib.Path('.').absolute()
        ( td / 'src' ).mkdir()
        ( td / 'src' / 'CMakeLists.txt' ).write_text(downstream_cmake)
        ( td  / 'bld' ).mkdir()
        os.chdir( td / 'bld' )
        runcmake('../src -DCMAKE_INSTALL_PREFIX=../install')
        runcmake('--build . --target install --config Release',apply_CMAKE_ARGS=False)
        os.chdir( td )
        app = ( td / 'install' / 'bin' / 'exampleapp' )
        if not app.exists():
            raise RuntimeError('Could not produce example app in downstream cmake project')
        prfct(f'Launching: {app}')
        ok, output = _run_cmd(str(app),env=runtime_env)
        if not ok:
            raise RuntimeError('Could not launch example app compiled in downstream cmake project')
        #envrun = os.environ.copy()
        expected_output='Powder Al x-sect at 1.8Aa is 1.44816barn'
        if not output.decode().strip()==expected_output:
            raise RuntimeError('Example app compiled in downstream cmake project produces unexpected output')
        prfct('App produced expected output')
    prfct('testing of compiled downstream cmake-based projects done')
