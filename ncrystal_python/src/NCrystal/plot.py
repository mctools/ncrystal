
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

Utility functions related to plotting.

"""

def plot_xsect( material, *, mode='wl', do_show = True, do_newfig = True,
                npts=5000, xmin = None, xmax = None, ymin=None, ymax=None,
                logx = None, logy=None,
                scatter_breakdown = True, xsmode='peratom',
                show_absorption = True, show_scattering = True,
                do_legend = True, do_grid = True, color = None,
                labelfct = None, only_total = False, extra_cfg=None,
                title = None ):
    """Quick plot (with matplotlib) showing the cross sections produced by
    current material. The mode parameter can be either 'wl' or 'ekin', and
    logx/logy can be set to None or a boolean to control whether the
    corresponding axis should be logarithmic (None indicates False when
    mode=='wl' and True when mode=='ekin'). Set scatter_breakdown to False
    to show the total scattering contribution rather than a breakdown into
    components (only affects materials that are not preloaded).

    The extra_cfg parameter can be used to append cfg parameters
    (e.g. "comp=bragg,incoh_elas;temp=50K").

    Use the plot_xsects (note the trailing "s") function instead to compare
    cross sections from multiple materials.

    """
    from . import misc as nc_misc
    matsrc = nc_misc.MaterialSource( material, cfg_params = extra_cfg )
    assert mode in ('wl','ekin')
    assert logx is None or logx in (True,False)
    assert xsmode in ('macroscopic','peratom')
    loadedmat = matsrc.load()
    if not loadedmat.scatter and not loadedmat.absorption:
        from .exceptions import NCBadInput
        raise NCBadInput('Can not produce plots for material source which contains no physics processes')
    do_absn = show_absorption and bool(loadedmat.absorption)
    do_scat = show_scattering and bool(loadedmat.scatter)
    labelfct = labelfct or ( lambda x : x )
    plt = _import_matplotlib_plt()
    if do_newfig:
        plt.figure()
    from ._numpy import _np_linspace, _np_geomspace, _ensure_numpy
    _ensure_numpy()
    from .constants import wl2ekin
    from ._common import _palette_Few as _palette
    wlmax = _estimate_longest_interesting_wavelength( loadedmat.info )
    if mode=='wl':
        xmin = 0.0 if xmin is None else xmin
        xmax = wlmax if xmax is None else xmax
        logxydefault = False
    else:
        xmin = ( wl2ekin(wlmax)*1e-3 ) if ( xmin is None ) else ( xmin or 1e-7 )
        xmax = xmax or xmin*1e5
        logxydefault = True
    assert xmax > xmin
    logx = logxydefault if logx is None else bool(logx)
    logy = logxydefault if logy is None else bool(logy)
    x = (_np_geomspace if logx else _np_linspace)( xmin, xmax, npts )
    xsectargs = dict(wl=x) if mode=='wl' else dict(ekin=x)
    if mode=='macroscopic' and not loadedmat.info:
        from .exceptions import NCBadInput
        raise NCBadInput('Can not produce macroscopic cross section'
                         ' plots for material source which contains'
                         ' no NCrystal.Info object (and therefore no material density).')
    xsectfactor = ( loadedmat.info.factor_macroscopic_xs ) if ( xsmode=='macroscopic' and loadedmat.info ) else 1.0
    do_scatter_breakdown = do_scat and scatter_breakdown and not matsrc.is_preloaded and not loadedmat.scatter.isNull()
    xs_s = xsectfactor * loadedmat.scatter.xsect(**xsectargs) if do_scat else None
    xs_a = xsectfactor * loadedmat.absorption.xsect(**xsectargs) if do_absn else None
    xs_tot = xs_s + xs_a if ( do_scat and do_absn ) else ( xs_s if do_scat else xs_a )
    nullxs = not ( xs_tot.max() > 0.0 )
    if nullxs and logy:
        from ._common import warn
        warn('Could not set log scale since curves are 0.0 everywhere')
        logy=False

    plotcalls = []
    def add_plot(x,y,*a,**kw):
        if y is not None:
            #if logy:
            #    nonzeroy = ( y > 0.0 )
            #    x,y = x[nonzeroy],y[nonzeroy]
            plotcalls.append( (x,y,a,kw) )
    if do_scatter_breakdown:
        lblmap = { 'coh_elas':'Coherent elastic',
                   'incoh_elas':'Incoherent elastic',
                   'inelas':'Inelastic',
                   'sans':'SANS'}
        colmap = { 'coh_elas':'blue',
                   'incoh_elas':'yellow',
                   'inelas':'green',
                   'sans':'pink'}#nb unused: brown, gray (but brown looks similar to yellow)

        assert set(lblmap.keys())==set(nc_misc.standard_comp_types),"update lblmap"
        assert set(colmap.keys())==set(nc_misc.standard_comp_types),"update colmap"
        for comp in nc_misc.detect_scattering_components( matsrc ):
            _xs_comp = xsectfactor * matsrc.load(f'comp={comp}',doInfo=False,doAbsorption=False).scatter.xsect(**xsectargs)
            add_plot(x,_xs_comp,label=labelfct(lblmap[comp]),
                     color = color or _palette[colmap[comp]] )
    else:
        add_plot(x,xs_s,label=labelfct('Scattering'),
                 color = color or _palette['blue'])
    add_plot(x,xs_a,label=labelfct('Absorption'),color = color or _palette['purple'])

    if only_total:
        plotcalls=[]
    if only_total or len(plotcalls) > 1:
        totlabel = 'Total' if ( do_absn and do_scat ) else ( 'Total scattering' if do_scat else 'Total absorption' )
        add_plot(x,xs_tot,label=labelfct(totlabel),
                 color=color or _palette['red'])
    for x,y,a,kw in plotcalls:
        plt.plot(x,y,*a,**kw)

    plt.xlabel('Neutron wavelength (angstrom)' if mode=='wl' else 'Neutron energy (eV)')
    plt.ylabel('Macroscopic cross section (1/cm)' if xsmode=='macroscopic' else 'Cross section per atom (barn)')

    auto_ymin = ( None if logy else 0.0)
    auto_ymax = ( 1.0 if nullxs else None )
    ymin = ( auto_ymin if ymin is None else ymin )
    ymax = ( auto_ymax if ymax is None else ymax )
    if ymin is not None or ymax is not None:
        if ymin is not None and ymax is not None and not ymax>ymin:
            from ._common import warn
            def _fmt(x):
                return ( 'auto' if x is None else x )
            warn('ymin/ymax parameters would lead to a plot range of'
                 f' [{ymin},{ymax}]. Reverting to [{_fmt(auto_ymin)},{_fmt(auto_ymax)}].')
            ymin,ymax = auto_ymin,auto_ymax
        plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    if title is not False:
        plt.title(title or matsrc.plotlabel)
    _plt_final(do_grid,do_legend,do_show,logy=logy,logx=logx,plt=plt)

def plot_xsects( *materials, **plotkwargs ):
    """Compares cross sections of multiple materials. Accepts similar keywords
    as the plot_xsect function."""
    from .misc import MaterialSource
    mats = [ MaterialSource(m) for m in materials ]
    do={}
    for e in ['do_newfig','do_legend','do_grid','do_show','logy','logx']:
        do[e] = plotkwargs.get(e,None if e in ('logy','logx') else True)
        plotkwargs[e] = False
    title = plotkwargs.get('title')
    plotkwargs['title'] = False
    plotkwargs['scatter_breakdown']=False
    plotkwargs['only_total'] = True
    plt = _import_matplotlib_plt()
    if do['do_newfig']:
        plt.figure()
    col_ordered = _get_col_ordered()
    for i,m in enumerate(mats):
        plotkwargs['color'] = col_ordered[i%len(col_ordered)]
        plotkwargs['labelfct'] = lambda _ : m.plotlabel
        if i+1 == len(mats):
            #the last:
            plotkwargs['logy'] = do['logy']
            plotkwargs['logx'] = do['logx']
        plot_xsect(m,**plotkwargs)
    if title:
        plt.title(title)
    _plt_final(do['do_grid'],do['do_legend'],do['do_show'],plt=plt)

def plot_vdos( *vdos, unit='meV',
               show_orig_data = False,
               show_normalised = True,
               do_show = True, do_newfig = True,
               npts_parabola=5000, logy=False,
               do_legend = True, do_grid = True, color = None,
               labelfct = None ):
    """
    Quick plot (with matplotlib) showing the requested VDOS curve(s).
    """
    from . import misc as nc_misc
    from . import vdos as nc_vdos
    #from ._common import _palette_Few as _palette
    from ._numpy import _np_linspace,_ensure_numpy, _np
    unit_name, unit_value = nc_vdos._parsevdosunit( unit )

    #from ._numpy import _np_linspace, _np_geomspace, _np, _ensure_numpy
    #_ensure_numpy()
    #from .constants import wl2ekin

    vdoslist = [ ( nc_misc.AnyVDOS(v),False) for v in vdos ]
    if show_orig_data:
        ll=[]
        for v,_ in vdoslist:
            ll.append( (v, False ) )
            if v.has_orig:
                ll.append( (v, True ) )
        vdoslist = ll

    plt = _import_matplotlib_plt()
    if do_newfig:
        plt.figure()

    col_ordered = _get_col_ordered()

    for ivdos,(vdos,use_orig) in enumerate(vdoslist):
        lbl = vdos.label or 'VDOS#%i'%(ivdos+1)
        if use_orig:
            lbl += ' (original)'
        if labelfct:
            lbl = labelfct( lbl )
        x = vdos.egrid( orig = use_orig ) / unit_value
        y = vdos.dos( orig = use_orig, norm = show_normalised ) * unit_value
        color = col_ordered[ivdos%len(col_ordered)]
        if x[0] > 0.0 and npts_parabola:
            _ensure_numpy()
            xp = _np_linspace( 0.0, x[0], npts_parabola )
            plt.plot(xp, y[0] *(xp/x[0])**2, color = color, ls = ':' )
        if y[-1] > 0.0:
            _ensure_numpy()
            x = _np.append( x, [ x[-1] ] )
            y = _np.append( y, [ 0.0 ] )
        plt.plot(x,y,color = color, label = lbl, marker='.' )

    plt.xlabel(f'Frequency ({unit_name})')
    plt.ylabel('VDOS (arbitrary units)')
    _plt_final(do_grid,do_legend,do_show,logy=logy,plt=plt)

def plot_knl( kernel, do_newfig = True, do_show = True, do_grid = True, logz=False, phasespace_curves = None,
              clim=None, xlim = None, ylim = None):
    """Quick plot (with matplotlib) showing the requested scattering kernel in
    in S(alpha,beta) format). The kernel is assume to be a dictionary with keys
    'alpha', 'beta', and 'sab' associated with values being numpy arrays
    containing the relevant info. Optionally, keys 'egrid' and 'suggestedEmax'
    can be present, and will be used to estimate Emax - the highest intended
    neutron energy for which the table is to be used. Emax is only used to plot
    the kinematic bound of Emax.
    """
    #estimate emax:
    _egrid = kernel.get('egrid')
    _suggestedEmax = kernel.get('suggestedEmax')
    emax = None
    if _suggestedEmax and _suggestedEmax>0.0:
        emax = _suggestedEmax
    if emax is None and _egrid is not None:
        import numbers
        if isinstance(_egrid,numbers.Real):
            emax = float(_egrid)
        elif len(_egrid)==1:
            emax = _egrid[0]
        elif len(_egrid)==3:
            emax = _egrid[1]
        elif len(_egrid)>3:
            emax = _egrid[-1]

    alpha = kernel['alpha']
    beta = kernel['beta']
    sab = kernel['sab']
    assert len(alpha)>=5
    assert len(beta)>=5
    assert len(alpha)*len(beta) == len(sab)

    from ._numpy import _ensure_numpy, _np, _np_linspace
    _ensure_numpy()
    plt = _import_matplotlib_plt()
    if do_newfig:
        plt.figure()


    def plotedges(x):
        """given array [a,b,c,...,x,y,z] returns [a,(a+b)2,(b+c)/2,...,(y+z)/2,z]"""
        return _np.concatenate([[x[0]],0.5*(x[:-1]+x[1::1]),[x[-1]]])
    na,nb = len(alpha),len(beta)
    pa = plotedges(alpha)
    pb = plotedges(beta)
    x,y=_np.meshgrid(pa,pb)
    if logz:
        sab = _np.log(sab)
    sab=sab.reshape((nb,na))

    cmap = 'jet'#fallback for some special options
    #sab_max = sab.max()
    quadmesh = plt.pcolormesh( x,y,sab, clim=clim,cmap=cmap )

    def _real_mpl_object( o ):
        #in case we are running under NCRYSTAL_FAKEPYPLOT=log, we are dealing
        #with a wrapped object - which we can not pass directly to matplotlib
        #classes.
        return o.the_real_inspected_object if hasattr(o,'the_real_inspected_object') else o

    plt.colorbar(_real_mpl_object(quadmesh))
    if clim:
        quadmesh.set_clim(clim)

    plt.xlim(*(xlim or (0.0,alpha[-1])))
    plt.ylim(*(ylim or (beta[0],beta[-1])))

    def alpha_limits(E_div_kT,beta):
        c = E_div_kT
        cb = c + beta
        assert cb>=0.0,f'c+beta={cb}, beta={beta}, c={c}'
        a = cb + c
        b = 2*_np.sqrt( c * cb )
        assert not _np.isinf(a).any()
        assert not _np.isinf(b).any()
        return max(0.0,a - b), a + b
    alpha_minus = _np.vectorize(lambda *a : alpha_limits(*a)[0])
    alpha_plus = _np.vectorize(lambda *a : alpha_limits(*a)[1])
    def kin_curve( E_div_kT ):
        b0, b1 = b0,b1=max(beta.min(),-E_div_kT),beta.max()
        #Smoother curve by putting more points around beta=0, less points
        #further out (actually, we could improve this alg when
        #ekin_div_kT>>1):
        if b0<1.0 and b1>1.0:
            bvals=_np.append(_np_linspace(b0,1.0,2500), _np_linspace(1.0+1e-14,b1,2500))
        else:
            bvals=_np_linspace(b0,b1,5000)
        color='yellow'
        plt.plot(alpha_minus(E_div_kT,bvals),bvals,label=f'E/kT={E_div_kT}',color=color,lw=2)
        plt.plot(alpha_plus(E_div_kT,bvals).clip(0,alpha.max()),bvals,color=color,lw=2)

    from .constants import constant_boltzmann
    kT = kernel['temperature']*constant_boltzmann

    for ekin in (phasespace_curves or []):
        kin_curve( ekin / kT )

    plt.xlabel('alpha')
    plt.ylabel('beta')
    _plt_final(do_grid=do_grid,do_legend=False,do_show=do_show,plt=plt)

def plot_vdos_Gn( Gn, unit = 'meV', logy = False,
                  do_newfig = True, do_legend = True,
                  do_grid = True, do_show = True ):
    """Plot Sjolander Gn functions. Each Gn must be specified as a sequence like
    (egrid,Gnvals,n[,label]), where either len(egrid)==len(Gnvals) or
    egrid=(emin,emax). The Gn argument can either be a single Gn function, or a
    sequence of them.
    """
    def decode_Gn( x ):
        if len(x) not in (3,4):
            return None
        if len(x[1]) <= 4:
            return None
        if len(x) == 3:
            egrid, gnvals, n = x
            label = None
        elif len(x) == 4:
            egrid, gnvals, n, label = x
        else:
            return None
        from ._numpy import _ensure_numpy,_np,_np_linspace
        _ensure_numpy()
        if len(egrid)==2 and len(gnvals)>2:
            egrid = _np_linspace(egrid[0],egrid[1],len(gnvals))
        return dict( egrid = _np.asarray(egrid,dtype=float),
                     gnvals = _np.asarray(gnvals,dtype=float),
                     n = n,
                     label = f'{label} (G{n})' if label else f'G{n}' )
    _ = decode_Gn( Gn )
    if _ is not None:
        gnlist = [ _ ]
    else:
        gnlist = [ decode_Gn(e) for e in Gn ]
    if not gnlist or any( e is None for e in gnlist ):
        from .exceptions import NCBadInput
        raise NCBadInput('Invalid specification of Gn functions.')

    from . import vdos as nc_vdos
    #from ._common import _palette_Few as _palette
    unit_name, unit_value = nc_vdos._parsevdosunit( unit )

    plt = _import_matplotlib_plt()
    if do_newfig:
        plt.figure()

    for gn in gnlist:
        plt.plot( gn['egrid']/unit_value, gn['gnvals'], 'o', label = gn['label'] )

    plt.xlabel(f'Energy ({unit_name})')
    _plt_final(do_grid,do_legend,do_show,logy=logy,plt=plt)

def _get_col_ordered():
    from ._common import _palette_Few
    return [_palette_Few.get(k,k) for k in
            ('blue','orange','green','red','brown',
             'purple','yellow','pink','gray','black')]

def _plt_final(do_grid,
               do_legend,
               do_show, *,
               logx=False,
               logy=False,
               plt=None,
               extra_legend_kwargs = None ):
    plt = plt or _import_matplotlib_plt()
    if logx:
        plt.semilogx()
    if logy:
        plt.semilogy()
    if do_legend:
        leg=plt.legend(**(extra_legend_kwargs or {}))
        if do_legend=='draggable':
            leg.set_draggable(True)
    if do_grid:
        plt.grid()
    if do_show:
        plt.show()

_fakepyplot_mode_cache = [None]
def _fakepyplot_mode():
    if _fakepyplot_mode_cache[0] is None:
        from ._common import ncgetenv
        _ = ncgetenv('FAKEPYPLOT','')
        if not _ or _=='0':
            res = False
        else:
            res = 'log_and_plot' if ( _ == 'log' ) else 'log_and_block'
        _fakepyplot_mode_cache[0] = res
    return _fakepyplot_mode_cache[0]

_theplt=[None]
def _import_matplotlib_plt():
    if not _theplt[0]:
        mode = _fakepyplot_mode()
        if mode:
            from ._testimpl import _create_pyplot_inspector
            _theplt[0] = _create_pyplot_inspector( pass_calls_to_real_plt = ( mode == 'log_and_plot' ) )
        else:
            #first matplotlib alone, for a perhaps more pedagogical ImportError.
            import matplotlib # noqa F401
            import matplotlib.pyplot as plt
            _theplt[0]=plt
    return _theplt[0]

_thepdfpages=[None]
def _import_matplotlib_pdfpages():
    if not _thepdfpages[0]:
        def _raw_import_pdf_pages():
            try:
                from matplotlib.backends.backend_pdf import PdfPages
            except ImportError:
                raise ImportError("ERROR: Your installation of matplotlib does not have the required support for PDF output.")
            return PdfPages
        mode = _fakepyplot_mode()
        if mode!='log_and_block':
            import matplotlib
            matplotlib.use('agg')
        if mode:
            from ._testimpl import _create_pdfpages_inspector
            _thepdfpages[0] = _create_pdfpages_inspector( _raw_import_pdf_pages() if mode != 'log_and_block' else None )
        else:
            _thepdfpages[0] = _raw_import_pdf_pages()
    return _thepdfpages[0]

def _find_highest_bragg_edge( info ):
    if info.isSinglePhase():
        return info.braggthreshold
    ll = [_find_highest_bragg_edge(p) for frac,p in info.phases]
    ll = [e for e in ll if e]
    return max(ll) if ll else None

def _estimate_longest_interesting_wavelength( info ):
    #longest wavelength of interest in material, for the purpose of plotting. In
    #case of bragg edges, this will be the longest bragg edge found in the
    #material.
    _ = _find_highest_bragg_edge( info ) if info else None
    return _*1.2 if _ else 15.0
