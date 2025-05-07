
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

"""

Internal implementation details for utilities in ncmat.py

"""

__all__ = []

def matsrc_detect_fmt( class_MaterialSource, data ):
    if isinstance(data,class_MaterialSource):
        return 'NCrystal.MaterialSource'
    from . import core as _nc_core
    for c in ('Info','Scatter','Absorption','LoadedMaterial','TextData'):
        if isinstance(data,getattr(_nc_core,c)):
            return f'NCrystal.{c}'
    if hasattr(data,'__fspath__'):
        return 'filepath'
    if isinstance(data,bytes) or isinstance(data,str):
        compat_str = (lambda s : s.encode()) if isinstance(data,bytes) else (lambda s : s)
        if compat_str('\n') in data and data.startswith(compat_str('NCMAT')):
            return 'ncmat'
        else:
            return 'cfgstr'
    if hasattr(data,'load') and hasattr(data,'create_ncmat'):
        from .ncmat import NCMATComposer as _composer
        if isinstance(data,_composer):
            return 'NCrystal.NCMATComposer'

def matsrc_allfmts( as_str = False ):
    _= ('cfgstr','filepath','ncmat','NCrystal.TextData','NCrystal.NCMATComposer',
        'NCrystal.Info','NCrystal.Scatter','NCrystal.Absorption',
        'NCrystal.LoadedMaterial','NCrystal.MaterialSource')
    return _ if not as_str else '"%s"'%('", "'.join(_))

def matsrc_initfmt( data, fmt, cfg_params ):
    from .exceptions import NCBadInput
    cfg_params = (cfg_params or '').strip() or None
    from . import core as _nc_core
    def to_str(s):
        return s.decode() if hasattr(s,'decode') else s
    def cfgcombine( cfgstr, extra_params ):
        ll = list( e.strip()
                   for e in (to_str(cfgstr),to_str(extra_params))
                   if (e or '').strip() )
        return ';'.join(ll) if ll else None
    if fmt == 'NCrystal.MaterialSource':
        if not cfg_params:
            return data
        wrapped_loadfct = data.get('loadfct')
        if not wrapped_loadfct:
            raise NCBadInput('Can not add cfg_params to MaterialSource'
                             ' which is of "preloaded" type.')
        def loadfct( extra_cfg, doInfo, doScatter, doAbsorption ):
            return wrapped_loadfct( cfgcombine(cfg_params,extra_cfg), doInfo, doScatter, doAbsorption )
        descr = data.get('description')
        plotlabel = data.get('plotlabel')
        if plotlabel:
            plotlabel+=f'+ ";{cfg_params}"'
        return dict( description = f'{descr} + ";{cfg_params}"',
                     plotlabel = plotlabel,
                     loadfct = loadfct )
    if fmt == 'filepath':
        import pathlib
        p = pathlib.Path(data)
        if not p.exists():
           raise NCBadInput(f'Missing file: {p}')
        try:
            _textdata = p.read_text()
        except OSError:
            raise NCBadInput(f'Could not read data from file: {p}')
        except UnicodeDecodeError:
            raise NCBadInput(f'Problems interpreting data in file as text: {p}')
        if _textdata.startswith('NCMAT'):
            _dtype = 'ncmat'
        else:
            _ = p.suffix
            if _ and _.startswith('.'):
                _ = _[1:]
            _dtype = _ or None
        descr = p.name
        if cfg_params:
            descr += ';' + cfg_params
        def loadfct( extra_cfg, doInfo, doScatter, doAbsorption ):
            return _nc_core.directMultiCreate( _textdata,
                                               cfg_params = cfgcombine(cfg_params,extra_cfg),
                                               dtype = _dtype,
                                               doInfo = doInfo,
                                               doScatter = doScatter,
                                               doAbsorption = doAbsorption )
        return dict( description = descr, loadfct = loadfct )
    if fmt=='cfgstr':
        cfgstr = to_str(data)
        if cfg_params:
            cfgstr += ';'+to_str(cfg_params)
        from .cfgstr import normaliseCfg
        cfgstr = normaliseCfg( cfgstr )
        def loadfct( extra_cfg, doInfo, doScatter, doAbsorption ):
            cs = cfgcombine(cfgstr,extra_cfg)
            i = _nc_core.createInfo(cs) if doInfo else None
            s = _nc_core.createScatter(cs) if doScatter else None
            a = _nc_core.createAbsorption(cs) if doAbsorption else None
            return _nc_core.LoadedMaterial.fromExistingObjects(info=i,scatter=s,absorption=a)
        return dict( description = cfgstr, loadfct = loadfct )
    if fmt=='NCrystal.TextData':
        descr = data.dataSourceName or 'Anonymous TextData'
        if cfg_params:
            descr += ';' + cfg_params
        def loadfct( extra_cfg, doInfo, doScatter, doAbsorption ):
            return _nc_core.directMultiCreate( data,
                                               cfg_params = cfgcombine(cfg_params,extra_cfg),
                                               doInfo = doInfo,
                                               doScatter = doScatter,
                                               doAbsorption = doAbsorption )
        return dict( description = descr, loadfct = loadfct )
    if fmt=='ncmat':
        data = to_str(data)
        def loadfct( extra_cfg, doInfo, doScatter, doAbsorption ):
            return _nc_core.directMultiCreate( data,
                                               dtype = 'ncmat',
                                               cfg_params = cfgcombine(cfg_params,extra_cfg),
                                               doInfo = doInfo,
                                               doScatter = doScatter,
                                               doAbsorption = doAbsorption )
        return dict( description = 'Anonymous NCMAT data', loadfct = loadfct )
    if fmt=='NCrystal.NCMATComposer':
        txtdata = data.create_ncmat()# NCMATComposer is not immutable so we simply
                                     # capture the state by invoking create_ncmat
                                     # already.
        plotlabel = data.plotlabel
        def loadfct( extra_cfg, doInfo, doScatter, doAbsorption ):
            return _nc_core.directMultiCreate( txtdata,
                                               dtype = 'ncmat',
                                               cfg_params = cfgcombine(cfg_params,extra_cfg),
                                               doInfo = doInfo,
                                               doScatter = doScatter,
                                               doAbsorption = doAbsorption )
        return dict( description = 'NCMATComposer object', loadfct = loadfct, plotlabel = plotlabel )
    if fmt=='NCrystal.LoadedMaterial':
        return dict( description = 'NCrystal.LoadedMaterial object',
                     preloaded = data )
    if fmt=='NCrystal.Info':
        return dict( description = 'NCrystal.Info object',
                     preloaded = _nc_core.LoadedMaterial.fromExistingObjects( info = data ) )
    if fmt=='NCrystal.Scatter':
        return dict( description = 'NCrystal.Scatter object',
                     preloaded = _nc_core.LoadedMaterial.fromExistingObjects( scatter = data ) )
    if fmt=='NCrystal.Absorption':
        return dict( description = 'NCrystal.Absorption object',
                     preloaded = _nc_core.LoadedMaterial.fromExistingObjects( absorption = data ) )
    raise _nc_core.NCLogicError(f'unhandled fmt option: {fmt}')

def matsrc_init( class_MaterialSource, matsrc_extract_d_fct, data, fmt, cfg_params ):
    from .exceptions import NCBadInput
    assert fmt is None or isinstance(fmt,str)
    assert data is not None
    assert cfg_params is None or isinstance(cfg_params,str)
    if not fmt:
        fmt = matsrc_detect_fmt( class_MaterialSource, data )
        if not fmt:
            raise NCBadInput('Could not detect format of MaterialSource.'
                             ' If this is merely a detection problem, set'
                             ' the fmt parameter to an accepted value (one of:'
                             ' %s)'%matsrc_allfmts(as_str=True))
    if fmt not in matsrc_allfmts():
        raise NCBadInput( f'Unknown format "{fmt}" (must be one of:'
                          ' %s)'%matsrc_allfmts(as_str=True) )
    d = matsrc_initfmt( matsrc_extract_d_fct(data) if fmt=='NCrystal.MaterialSource' else data, fmt, cfg_params )
    assert len(d) in (2,3)
    assert 'description' in d
    assert int('preloaded' in d) + int('loadfct' in d) == 1
    if 'preloaded' in d and (cfg_params or '').strip():
        raise NCBadInput('MaterialSource objects that are wrapping preloaded'
                         ' physics objects can not be initialised with'
                         ' additional cfg parameters.')
    return d

def anytextdata_init( data, *, is_path, name ):
    if isinstance(name,bytes):
        name = name.decode()
    assert name is None or isinstance(name,str)
    res_keepalive = None
    if hasattr(data,'rawData') and hasattr(data,'dataSourceName'):
        #TextData:
        res_keepalive = data
        res_name = name or data.dataSourceName
        res_content = data.rawData
        return res_name, res_content, res_keepalive
    if hasattr( data, '__fspath__' ):
        import pathlib
        p = pathlib.Path(data)
        res_name = name or p.name
        res_content = p.read_text()
        return res_name, res_content, res_keepalive
    if isinstance(data,bytes):
        data = data.decode()
    if not isinstance(data,str):
        from .exceptions import NCBadInput
        raise NCBadInput('Invalid text data / text file data (got type %s)'%type(data))
    is_path = is_path if ( is_path is not None ) else ( '\n' not in data )
    if is_path:
        from ._common import _lookup_existing_file
        p = _lookup_existing_file( data )
        if not p:
            raise NCBadInput(f'Path not found: {data}')
        res_name = name or p.name
        res_content = p.read_text()
    else:
        res_content = data
        res_name = name#Ok if None
    return res_name or None, res_content, res_keepalive

def _anyvdos_detect_fmt( class_AnyVDOS, data ):
    if isinstance(data,class_AnyVDOS):
        return 'NCrystal.AnyVDOS'
    from .core import Info
    if isinstance( data, Info.DI_VDOS ):
        return 'NCrystal.Info.DI_VDOS'
    if isinstance( data, Info.DI_VDOSDebye ):
        return 'NCrystal.Info.DI_VDOSDebye'
    import numbers
    if isinstance(data, numbers.Real):
        return 'debyetemperature'
    if hasattr(data,'__len__') and len(data)==2:
        eg,dens = data
        if hasattr(eg,'__len__') and len(eg) >= 2 and hasattr(dens,'__len__') and len(dens)>=2 and ( len(eg)==len(dens) or len(eg)==2 ):
            return 'arrays'

def _anyvdos_preinit( data, fmt ):
    allfmts = ('NCrystal.Info.DI_VDOS','NCrystal.Info.DI_VDOSDebye','debyetemperature','arrays','NCrystal.AnyVDOS')
    from ._numpy import _ensure_numpy, _np
    _ensure_numpy()

    if fmt =='NCrystal.Info.DI_VDOS':
        return dict( dos_orig = ( data.vdosOrigEgrid(), data.vdosOrigDensity() ),
                     dos = data.vdosData() )
    if fmt =='NCrystal.Info.DI_VDOSDebye':
        return dict( dos = data.vdosData(),
                     debyetemp = data.debyeTemperature() )
    if fmt =='debyetemperature':
        debyetemp = float(data)
        from .vdos import createVDOSDebye
        return dict( dos = createVDOSDebye( debyetemp ),
                     debyetemp = debyetemp )
    if fmt =='arrays':
        eg, dens = data
        dens = _np.asarray( dens, dtype = float )
        if len(eg)==2 and len(dens)>2:
            eg = (float(eg[0]),float(eg[1]))
        else:
            eg = _np.asarray( eg, dtype = float )
        return dict( dos = ( eg, dens ) )
    assert not fmt == 'NCrystal.AnyVDOS'#should have been handled in calling code
    if fmt not in allfmts:
        from .exceptions import NCBadInput
        s='"%s"'%('", "'.join(allfmts))
        if fmt is None:
            raise NCBadInput('Could not detect fmt of VDOS specification, consider specifying it explictly (must be one of %s)'%(s))
        raise NCBadInput('Invalid VDOS fmt "%s". Possible options are %s'%(fmt,s))
    assert False, "unhandled fmt case: %s"%fmt

def _anyvdos_initfmt( data, fmt ):
    p = _anyvdos_preinit( data, fmt )
    #find derived quantities:
    d = {}
    from ._numpy import _ensure_numpy, _np_linspace, _np, _np_trapezoid
    _ensure_numpy()
    def _needsexpand( egrid, dos ):
        return len(egrid)==2 and len(dos)>2
    def _expandegrid( egrid, dos ):
        if _needsexpand(egrid,dos):
            return _np_linspace( egrid[0], egrid[-1], len(dos) )
        else:
            return _np.asarray(egrid,dtype=float)
    def integrate( expanded_egrid, dos ):
        x,y = expanded_egrid, dos
        if not x[0] > 0.0:
            #simply ignore points <= 0.0
            i = _np.argmax(x>0.0)
            x,y = x[i:],y[i:]
        parabola_contrib = ( x[0] * y[0] ) / 3.0#integral of parabola through (0,0) and (x[0],y[0]) over [0,x[0]]
        return parabola_contrib + _np_trapezoid(x=x,y=y)
    for k in ('dos','dos_orig'):
        if k not in p:
            continue
        eg,dos = p[k]
        egexpand = _expandegrid(eg,dos) if _needsexpand(eg,dos) else None
        d[k+'_expandedegrid'] = egexpand
        d[k+'_integral'] = integrate( ( eg if egexpand is None else egexpand), dos )

    #combine and return:
    import copy
    res = {'derived':d}
    for k,v in p.items():
        res[k] = copy.deepcopy(v)
    return res

def _anyvdos_init( class_AnyVDOS, anyvdos_extract_d, data, fmt, label ):
    if label is None and hasattr(data,'_plotlabel'):
        from .core import Info
        if isinstance(data,Info.DynamicInfo):
            label = data._plotlabel()
    fmt = _anyvdos_detect_fmt(class_AnyVDOS,data) if fmt is None else fmt
    if fmt == 'NCrystal.AnyVDOS':
        d = anyvdos_extract_d( data )
        if label is None or d.get('label')==label:
            return d
        else:
            import copy
            d = dict( (k,copy.deepcopy(v)) for k,v in d.items() if k != 'label' )
    else:
        d = _anyvdos_initfmt( data, fmt )
    d['label'] = label or 'anonymous VDOS curve from "%s"'%fmt
    if d.get('debyetemp'):
        _='TDebye=%gK'%d['debyetemp']
        if not d['label'].endswith(')'):
            d['label'] += f' ({_})'
        else:
            d['label'] = d['label'][:-1]+f', {_})'
    from types import MappingProxyType#read-only dict
    d['derived'] = MappingProxyType(d['derived'])
    return MappingProxyType( d )

def detect_scatcomps( standard_comp_types, matsrc ):
    from .exceptions import NCBadInput, NCCalcError
    if matsrc.is_preloaded:
        from .exceptions import NCBadInput
        raise NCBadInput('detect_scattering_components can not be used'
                         ' with preloaded material sources')
    res=[]
    for ct in standard_comp_types:
        def load(extra = None):
            p = f'comp={ct}'
            if extra is not None:
                p += f';{extra}'
            return matsrc.load( extra_cfg_params = p,
                                doInfo = False,
                                doAbsorption = False )
        m = None
        if ct == 'inelas':
            #For inelas component first we try to see if we can get away with
            #using a potentially smaller vdoslux than is actually used. This is
            #not a great solution, but in the absence of deeper changes to
            #NCrystal processes, this seems a reasonable compromise.
            try:
                m = load('vdoslux=0')
            except NCCalcError as e:
                if not str(e).startswith('VDOS expansion too slow'):
                    raise e
                m = None
        if not m:
            m = load()
        if not m.scatter.isNull():
            res.append(ct)
    return res
