
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

"""

Utility functions related to plotting.

"""

__all__=['plot_xsect','plot_xsects','plot_vdos','plot_knl','PlotContext']

def plot_xsect( material, *, mode='wl',
                npts=5000, xmin = None, xmax = None, ymin=None, ymax=None,
                logx = None, logy=None,
                scatter_breakdown = True, xsmode='peratom',
                show_absorption = True, show_scattering = True,
                do_legend = True, do_grid = True, color = None,
                labelfct = None, only_total = False, extra_cfg=None,
                title = None, **kw_plot ):
    """Quick plot (with matplotlib) showing the cross sections produced by
    the provided material. The mode parameter can be either 'wl' or 'ekin', and
    logx/logy can be set to None or a boolean to control whether the
    corresponding axis should be logarithmic (None indicates False when
    mode=='wl' and True when mode=='ekin'). Set scatter_breakdown to False
    to show the total scattering contribution rather than a breakdown into
    components (only affects materials that are not preloaded).

    The extra_cfg parameter can be used to append cfg parameters
    (e.g. "comp=bragg,incoh_elas;temp=50K").

    Use the plot_xsects (note the trailing "s") function instead to compare
    cross sections from multiple materials.

    Other arguments are hopefully self-explanatory, and any excess arguments are
    passed along to the PlotContext.

    """
    from . import misc as nc_misc
    from ._numpy import _np_linspace, _np_geomspace, _ensure_numpy
    from .constants import wl2ekin
    from ._common import _palette_Few as _palette

    _ensure_numpy()
    matsrc = nc_misc.MaterialSource( material, cfg_params = extra_cfg )
    assert mode in ('wl','ekin')
    assert logx is None or logx in (True,False)
    assert xsmode in ('macroscopic','peratom')
    loadedmat = matsrc.load()
    if not loadedmat.scatter and not loadedmat.absorption:
        from .exceptions import NCBadInput
        raise NCBadInput('Can not produce plots for material source which'
                         ' contains no physics processes')
    do_absn = show_absorption and bool(loadedmat.absorption)
    do_scat = show_scattering and bool(loadedmat.scatter)
    labelfct = labelfct or ( lambda x : x )
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
    if xsmode=='macroscopic' and not loadedmat.info:
        from .exceptions import NCBadInput
        raise NCBadInput('Can not produce macroscopic cross section'
                         ' plots for material source which contains'
                         ' no NCrystal.Info object (and therefore no'
                         ' material density).')
    xsectfactor = ( ( loadedmat.info.factor_macroscopic_xs )
                    if ( xsmode=='macroscopic' and loadedmat.info ) else 1.0 )
    do_scatter_breakdown = ( do_scat and scatter_breakdown
                             and not matsrc.is_preloaded
                             and not loadedmat.scatter.isNull() )
    xs_s = ( xsectfactor * loadedmat.scatter.xsect(**xsectargs)
             if do_scat else None )
    xs_a = ( xsectfactor * loadedmat.absorption.xsect(**xsectargs)
             if do_absn else None )
    xs_tot = ( (xs_s + xs_a ) if ( do_scat and do_absn )
               else ( xs_s if do_scat else xs_a ) )
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
    if not only_total:
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
                _xs_comp = xsectfactor * matsrc.load(f'comp={comp}',
                                                     doInfo=False,
                                                     doAbsorption=False).scatter.xsect(**xsectargs)
                add_plot(x,_xs_comp,label=labelfct(lblmap[comp]),
                         color = color or _palette[colmap[comp]] )
        else:
            add_plot(x,xs_s,label=labelfct('Scattering'),
                     color = color or _palette['blue'])
        add_plot(x,xs_a,label=labelfct('Absorption'),color = color or _palette['purple'])
    if only_total or len(plotcalls) > 1:
        totlabel = 'Total' if ( do_absn and do_scat ) else ( 'Total scattering' if do_scat else 'Total absorption' )
        add_plot(x,xs_tot,label=labelfct(totlabel),
                 color=color or _palette['red'])
    pctx = PlotContext(**kw_plot).check_unused()
    for x,y,a,kw in plotcalls:
        pctx.axis.plot(x,y,*a,**kw)

    pctx.axis.set_xlabel('Neutron wavelength (angstrom)' if mode=='wl' else 'Neutron energy (eV)')
    pctx.axis.set_ylabel('Macroscopic cross section (1/cm)' if xsmode=='macroscopic' else 'Cross section per atom (barn)')

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
        pctx.axis.set_ylim(ymin,ymax)
    pctx.axis.set_xlim(xmin,xmax)
    if title is not False:
        pctx.axis.set_title(title or matsrc.plotlabel)

    return pctx.finalise( do_grid = do_grid,
                          do_legend = do_legend,
                          logx=logx,
                          logy=logy )

def plot_xsects( *materials, **plotkwargs ):
    """Compares cross sections of multiple materials. Accepts similar keywords
    as the plot_xsect function."""
    #TODO: x/y limits are not working very well, since some of the logic resides
    #in the plot_xsect function.
    from .misc import MaterialSource
    pctx = PlotContext( **plotkwargs )
    plotkwargs = pctx.kwargs_unused

    mats = [ MaterialSource(m) for m in materials ]
    do={}
    for e in ['do_legend','do_grid','logy','logx']:
        do[e] = plotkwargs.get(e,None if e in ('logy','logx') else True)
        plotkwargs[e] = False
    title = plotkwargs.get('title')
    plotkwargs['title'] = False
    plotkwargs['scatter_breakdown']=False
    plotkwargs['only_total'] = True
    col_ordered = _get_col_ordered()
    for i,m in enumerate(mats):
        plotkwargs['color'] = col_ordered[i%len(col_ordered)]
        plotkwargs['labelfct'] = lambda _ : m.plotlabel
        if i+1 == len(mats):
            #the last:
            plotkwargs['logy'] = do['logy']
            plotkwargs['logx'] = do['logx']
        plot_xsect(m,**plotkwargs,**pctx.kwargs_subcontext())
    if title:
        pctx.axis.set_title(title)
    return pctx.finalise( do_grid = do['do_grid'],
                          do_legend = do['do_legend'] )

def plot_vdos( *vdos, unit='meV',
               show_orig_data = False,
               show_normalised = True,
               npts_parabola=5000, logy=False,
               do_legend = True, do_grid = True, color = None,
               labelfct = None, **kw_plot ):
    """Quick plot (with matplotlib) showing the requested VDOS curve(s).

    Most arguments are hopefully self-explanatory, and any excess arguments are
    passed along to the PlotContext.
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

    pctx = PlotContext(**kw_plot).check_unused()

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
            pctx.axis.plot(xp, y[0] *(xp/x[0])**2, color = color, ls = ':' )
        if y[-1] > 0.0:
            _ensure_numpy()
            x = _np.append( x, [ x[-1] ] )
            y = _np.append( y, [ 0.0 ] )
        pctx.axis.plot(x,y,color = color, label = lbl,
                       marker=('o' if use_orig else '.') )

    pctx.axis.set_xlabel(f'Frequency ({unit_name})')
    pctx.axis.set_ylabel('VDOS (arbitrary units)')
    return pctx.finalise(do_grid=do_grid,
                         do_legend=do_legend,
                         logy=logy)

def plot_knl( kernel, do_grid = True, logz=False, phasespace_curves = None,
              clim=None, xlim = None, ylim = None, cmap = 'jet', **kw_plot ):
    """Quick plot (with matplotlib) showing the requested scattering kernel in
    in S(alpha,beta) format). The kernel is assume to be a dictionary with keys
    'alpha', 'beta', and 'sab' associated with values being numpy arrays
    containing the relevant info. Optionally, keys 'egrid' and 'suggestedEmax'
    can be present, and will be used to estimate Emax - the highest intended
    neutron energy for which the table is to be used. Emax is only used to plot
    the kinematic bound of Emax.

    Most arguments are hopefully self-explanatory, and any excess arguments are
    passed along to the PlotContext.
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
    pctx = PlotContext( **kw_plot ).check_unused()

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

    #sab_max = sab.max()
    quadmesh = pctx.axis.pcolormesh( x,y,sab, clim=clim,cmap=cmap )

    def _real_mpl_object( o ):
        #in case we are running under NCRYSTAL_FAKEPYPLOT=log, we are dealing
        #with a wrapped object - which we can not pass directly to matplotlib
        #classes.
        return ( o.the_real_inspected_object
                 if hasattr(o,'the_real_inspected_object')
                 else o )

    fig = pctx.axis.get_figure() if hasattr(pctx.axis,'get_figure') else None
    if fig:
        fig.colorbar(_real_mpl_object(quadmesh))
    else:
        from ._common import warn
        warn('Could not get figure from axis. Will not add colorbar')
    if clim:
        quadmesh.set_clim(clim)

    pctx.axis.set_xlim(*(xlim or (0.0,alpha[-1])))
    pctx.axis.set_ylim(*(ylim or (beta[0],beta[-1])))

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
        pctx.axis.plot(alpha_minus(E_div_kT,bvals),bvals,label=f'E/kT={E_div_kT}',color=color,lw=2)
        pctx.axis.plot(alpha_plus(E_div_kT,bvals).clip(0,alpha.max()),bvals,color=color,lw=2)

    from .constants import constant_boltzmann
    kT = kernel['temperature']*constant_boltzmann

    for ekin in (phasespace_curves or []):
        kin_curve( ekin / kT )

    pctx.axis.set_xlabel('alpha')
    pctx.axis.set_ylabel('beta')
    return pctx.finalise(do_grid=do_grid,do_legend=False)

def plot_vdos_Gn( Gn, unit = 'meV', logy = False,
                  do_legend = True, do_grid = True, **kw_plot ):
    """Plot Sjolander Gn functions. Each Gn must be specified as a sequence like
    (egrid,Gnvals,n[,label]), where either len(egrid)==len(Gnvals) or
    egrid=(emin,emax). The Gn argument can either be a single Gn function, or a
    sequence of them.

    Most arguments are hopefully self-explanatory, and any excess arguments are
    passed along to the PlotContext.
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

    pctx = PlotContext( **kw_plot ).check_unused()

    for gn in gnlist:
        pctx.axis.plot( gn['egrid']/unit_value, gn['gnvals'], 'o',
                        label = gn['label'] )

    pctx.axis.set_xlabel(f'Energy ({unit_name})')
    return pctx.finalise(do_grid=do_grid,do_legend=do_legend,logy=logy)

def _get_col_ordered():
    from ._common import _palette_Few
    return [_palette_Few.get(k,k) for k in
            ('blue','orange','green','red','brown',
             'purple','yellow','pink','gray','black')]

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
            _theplt[0] = (
                _create_pyplot_inspector(
                    pass_calls_to_real_plt = ( mode == 'log_and_plot' )
                )
            )
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
                raise ImportError("ERROR: Your installation of matplotlib does"
                                  " not have the required support for PDF"
                                  " output.")
            return PdfPages
        mode = _fakepyplot_mode()
        if mode!='log_and_block':
            import matplotlib
            matplotlib.use('agg')
        if mode:
            from ._testimpl import _create_pdfpages_inspector
            _thepdfpages[0] = (
                _create_pdfpages_inspector(
                    _raw_import_pdf_pages() if mode != 'log_and_block' else None
                )
            )
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

class PlotContext:
    """Class which can be used to control the behaviour of various plotting
    functions in NCrystal. By default, these plots will simply use
    plt=matplotlib.pyplot to create a new figure with plt.subplots() and plot
    onto the resulting axis object. Finally, plt.show() will be called.

    This is often a suitable behaviour when using Python interactively, but it
    might not always be the case. For instance, one might with to change this
    when:

    1. Wanting to collect multiple items into the same plot (so a new
       figure/axis should not be created, and plt.show() should not be called).
    2. Wanting to embed plots into e.g. a Qt canvas.
    3. Wanting a nice experience when plotting in a matplotlib notebook.

    We provide this PlotContext utility to help with that. For now this is an
    under-developed and under-documented feature.
    """
    #TODO: add support for pdf output as well, and make sure we migrate ALL
    #NCrystal code (including in cifutils.py, _cli_nctool.py) to using this
    #PlotContext system.
    def __init__( self,
                  plot_context = None,
                  plt = None,
                  axis = None,
                  do_show = None,
                  **excess_args ):
        """Provide at most one of: plot_context, plt, and axis. Do not provide
           do_show with plot_context. If no arguments are given, the
           matplotlib.pyplot module will be used as the plt object, a new figure
           and axis will be created with plt.subplots(), the plots will be added
           into that axis, and at the end (when .finalise() is invoked),
           plt.show() will be called. The various arguments allow to override
           this behaviour.
        """
        if plot_context=='__no_init__':
            return

        def err(msg):
            from .exceptions import NCBadInput
            raise NCBadInput(msg)

        if plot_context is not None:
            if do_show is not None:
                return err('explicit do_show not supported with plot_context')
            if axis is not None:
                return err('do not provide both plot_context and axis')
            if plt is not None:
                return err('do not provide both plot_context and plt')
            #Done, just adopt contents ("move semantics"):
            self.__axisgen = plot_context.__axisgen
            self.__pltgen = plot_context.__pltgen
            self.__do_show = plot_context.__do_show
            self.__excess_args = plot_context.__excess_args
            plot_context._clear_data()
            return

        if (plt is not None) and (axis is not None):
            return err('do not provide both axis and plt')

        #Setup on-demand access to plt object:
        if plt:
            def plt_gen():
                return plt
        else:
            _cache_pltgen = [None]
            def plt_gen():
                if _cache_pltgen[0] is None:
                    _cache_pltgen[0] = _import_matplotlib_plt()
                return _cache_pltgen[0]
        #Setup on-demand access to axis object (might trigger use of plt to
        #create new figure with plt.subplots()):
        if axis:
            def axis_gen():
                return axis
        else:
            _cache_axisgen = [None]
            def axis_gen():
                if _cache_axisgen[0] is None:
                    plt = plt_gen()
                    subplotsfct = getattr(plt,'subplots',None)
                    if not subplotsfct:
                        from .exceptions import NCLogicError
                        raise NCLogicError('plt.subplots unavailable')
                    if hasattr(plt,'the_real_inspected_object'):
                        #unit tests
                        plt.subplots()
                        _axis = plt.gca()
                    else:
                        _fig, _axis = plt.subplots()
                    _cache_axisgen[0] = _axis
                return _cache_axisgen[0]

        if do_show is None:
            #if we got plt or nothing, default is to show. If we got axis,
            #default is to not show.
            do_show = bool(axis is None)

        #All done, setup data:
        self.__axisgen = axis_gen
        self.__pltgen = plt_gen
        self.__do_show = do_show
        self.__excess_args = excess_args

    def _clear_data(self):
        self.__axisgen = None
        self.__pltgen = None
        self.__do_show = None
        self.__excess_args = None

    def subcontext( self ):
        """Create a subcontext. This will have access to the same plt+axis
        objects, but will never possess any .kwargs_unused and will never
        trigger plt.show(). If a plotting function needs to call other plotting
        functions, it should pass along a subcontext of its own PlotContext
        object.
        """
        self._check_alive()
        c = PlotContext(plot_context='__no_init__')
        c.__axisgen = self.__axisgen
        c.__pltgen = self.__pltgen
        c.__do_show = False
        c.__excess_args = {}
        return c

    def kwargs_subcontext( self ):
        """Simply returning {'plot_context':obj.subcontext()} for convenience.
        """
        return dict( plot_context=self.subcontext() )

    def check_unused( self ):
        """If there should be no .kwargs_unused to pass along to other
        functions, this function can be used to check that all arguments are
        consumed. Intended usage:

        pctx = PlotContext(**kw_plot).check_unused()
        """
        for k,v in sorted(self.__excess_args.items()):
            from .exceptions import NCBadInput
            raise NCBadInput(f'Unsupported argument: {k}={repr(v)}')
        return self

    @property
    def kwargs_unused( self ):
        self._check_alive()
        return self.__excess_args

    @property
    def axis( self ):
        self._check_alive()
        return self.__axisgen()

    def _check_alive(self):
        if self.__axisgen is None:
            from .exceptions import NCLogicError
            raise NCLogicError(
                'Do not use PlotContext object after it has already been'
                ' consumed. It is consumed when its .finalise() method is'
                ' invoked or when it is passed directly as the value of a'
                ' plot_context argument (you most likely want to pass it as'
                ' plot_context=obj.subcontext() instead in this case).')

    def finalise( self, *,
                  do_grid = False,
                  do_legend = False,
                  logx=False,
                  logy=False,
                  extra_legend_kwargs = None ):
        """Finish up plot and possibly invoke plt.show depending on context."""
        #Consume (to avoid future usage and to release refs):
        self._check_alive()
        do_show, pltgen, axisgen = self.__do_show, self.__pltgen, self.__axisgen
        self._clear_data()
        if logx:
            axisgen().semilogx()
        if logy:
            axisgen().semilogy()
        if do_legend:
            leg=axisgen().legend(**(extra_legend_kwargs or {}))
            if do_legend=='draggable':
                leg.set_draggable(True)
        if do_grid:
            axisgen().grid()
        plt=None
        if do_show:
            plt = pltgen()
            plt.show()
