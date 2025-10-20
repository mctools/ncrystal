
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

# For backwards compatibility we keep __all__ here in this internal file. In
# NCrystal 4.0.0 we actually renamed the public API file from _mmc.py to mmc.py.

__all__ = ['runsim_diffraction_pattern',
           'quick_diffraction_pattern',
           'quick_diffraction_pattern_autoparams']

def runsim_diffraction_pattern( cfgstr, *,
                                geomcfg,
                                srccfg,
                                nthreads = 'auto',
                                tally_detail_lvl = 2 ):
    """Run embedded "mini-MonteCarlo" to produce diffraction pattern, including
    both effects of multiple scattering and absorption physics. This for now
    results in exit angle (angle of emitted neutrons wrt. the Z-axis)
    histograms, with a perfect 4pi detector.

    This is highly experimental, and not yet fully documented.

    Example geomcfg: 'sphere;r=0.1'. This describes a 0.1m=10cm radius sphere
    centered at (0,0,0).

    Example srccfg: 'constant;ekin=0.025;n=1e6;z=-0.1'. This starts 1e6 0.025eV
    neutrons at (0,0,-10cm) with a direction (0,0,1).

    tally_detail_lvl can be reduced to 1 or 0, if only the exit_angle histogram
    is needed. tally_detail_lvl=2 provides more details, including histograms
    for specific components (e.g. separating the contributions from single and
    multiple scattering, and inelastic/elastic scatterings.

    nthreads can be 'auto' or a specific integral value of threads to use.

    """
    assert tally_detail_lvl in (0,1,2)
    from ._numpy import _np, _ensure_numpy, _np_linspace
    _ensure_numpy()

    nthreads = 9999 if nthreads=='auto' else min(9999,max(1,int(nthreads)))

    setup_info = dict( nthreads = nthreads,
                       tally_detail_lvl = tally_detail_lvl,
                       material_cfgstr = cfgstr,
                       geomcfg = geomcfg,
                       srccfg = srccfg )

    from ._chooks import _get_raw_cfcts
    _rawfct = _get_raw_cfcts()
    cfct = _rawfct['runmmcsim_stdengine']

    ct,errsq,res_json = cfct( nthreads = nthreads,
                              tally_detail_lvl = int(tally_detail_lvl),
                              mat_cfgstr = str(cfgstr),
                              mmc_geomcfg = str(geomcfg),
                              mmc_srccfg = str(srccfg) )

    class Results:

        def __init__( self,
                      main_hist_content,
                      main_hist_errsq,
                      json_details,
                      cfgstr,
                      setup_info ):

            self.__setup_info = setup_info

            mainhist_dict = dict( title = 'MAIN',
                                  bindata = dict( content = main_hist_content,
                                                  errorsq = main_hist_errsq,
                                                  nbins = len(main_hist_content),
                                                  xmin = 0.0,
                                                  xmax = 180.0 ),
                                  stats = None )

            datadict = {}
            if json_details is not None:
                import json
                datadict = json.loads(json_details)
                mainhist_dict['stats'] = datadict.get('main_stats')

            from ._hist import Hist1D
            hists = [ Hist1D( mainhist_dict ) ]
            hists += [ Hist1D( h )
                       for h in (datadict.get('breakdown_hists',[])) ]
            self.__hists = hists
            self.__cfgstr = cfgstr

        @property
        def histograms( self ):
            return self.__hists

        @property
        def histogram_main( self ):
            return self.__hists[0]

        @property
        def histogram_breakdown( self ):
            return dict( (h.title,h) for h in self.__hists[1:] )

        def histogram_sum( self, *, select=None, exclude=None ):
            if isinstance(exclude,str):
                exclude=[exclude]
            if isinstance(select,str):
                select=[select]
            hl = self.__hists[1:]
            if not exclude and not select:
                return self.histogram_main
            if select:
                hl = [ h for h in hl if h.title in select ]
            if exclude:
                hl = [ h for h in hl if h.title not in exclude ]
            if len(hl) <= 1:
                return hl[0] if hl else None
            h = hl[0].clone()
            for o in hl[1:]:
                h.add_contents( o )
            return h

        @property
        def cfgstr( self ):
            return self.__cfgstr

        @property
        def setup_info( self ):
            return self.__setup_info

        def plot_main( self, **kwargs ):
            self.plot(hist=0,**kwargs)

        def plot_breakdown( self, **kwargs ):
            self.plot(hist='breakdown',**kwargs)

        def plot( self, *,
                  do_show = True,
                  do_newfig = True,
                  do_grid = False,
                  logy = False,
                  rebin_factor = 1,
                  hist = None,
                  title = None,
                  plt = None):

            breakdown_mode = hist=='breakdown'

            if breakdown_mode:
                hist = None
            else:
                from ._hist import Hist1D
                if hist is None:
                    hist = self.__hists[0]
                    assert isinstance(hist,Hist1D)
                elif not isinstance(hist,Hist1D):
                    #must be index:
                    hist = self.__hists[hist]
                assert isinstance(hist,Hist1D)

            from .plot import _import_matplotlib_plt, _plt_final
            plt = plt or _import_matplotlib_plt()
            if do_newfig:
                plt.figure()

            do_legend = False
            tot_integral_with_absorption = None#unknown

            if breakdown_mode:

                colors = {
                    'SINGLESCAT_ELAS' : 'blue',
                    'SINGLESCAT_INELAS' : 'orangered',
                    'MULTISCAT_PUREELAS' : 'cornflowerblue',
                    'MULTISCAT_OTHER' : 'orange',
                    'NOSCAT' : 'green',
                }

                do_legend = True
                hists = [h for h in self.histograms if h.title!='MAIN']
                hist_main = [h for h in self.histograms if h.title=='MAIN'][0]
                if not hists:
                    from .exceptions import NCCalcError
                    raise NCCalcError('Can not produce breakdown plots'
                                      ' without more detailed tally.')
                _hists = [h for h in hists if h.stats.get('integral',1.0) > 0.0]
                if _hists:
                    hists=_hists
                #Sort histograms for (possibly) more meaningful plot order:
                hists.sort( key = lambda h: ['MULTISCAT_OTHER',
                                             'SINGLESCAT_INELAS',
                                             'MULTISCAT_PUREELAS',
                                             'SINGLESCAT_ELAS',
                                             'NOSCAT'].index(h.title))

                integrals = [ h.stats.get('integral',-1.0) for h in hists ]
                tot_integral = ( sum(integrals)
                                 if all( i>=0.0 for i in integrals )
                                 else None )
                if 'nstat' in self.setup_info:#TODO: we should always have this!
                    #A bit hackish way to get the initial flux... works only
                    #with the current unweighted generators and quick mode:
                    tot_integral_with_absorption = self.setup_info['nstat']

                if rebin_factor != 1:
                    _ = []
                    for h in hists:
                        h = h.clone( rebin_factor = rebin_factor )
                        _.append( h )
                    hists = _
                    hist_main = hist_main.clone( rebin_factor = rebin_factor)

                h=hists[0].clone()
                curve = h._hist_curve()
                axis = plt.gca()
                ymax_non_NOSCAT = [0.0]
                def _fractionval_fmt(x):
                    return f'({x*100.0:.3g}%)'
                def _addplot( title, curve, refcurve = None, fraction = None):
                    maxval = 0.0
                    if refcurve is None:
                        y2 = _np.zeros(len(curve[1]))
                    else:
                        y2=refcurve[1]
                        maxval = y2.max()
                    if title!='NOSCAT':
                        ymax_non_NOSCAT[0] = max(ymax_non_NOSCAT[0],
                                                 maxval,
                                                 curve[1].max())
                    lbl = title if title!='NOSCAT' else 'Transmitted'
                    if fraction is not None:
                        lbl = f'{lbl} {_fractionval_fmt(fraction)}'
                    axis.fill_between(*curve,y2,
                                      label=lbl,
                                      edgecolor="none",
                                      facecolor=colors[title])



                def calcfrac(h):
                    ti = tot_integral
                    if tot_integral_with_absorption is not None:
                        ti = tot_integral_with_absorption
                    if ti is not None and ti > 0.0:
                        return h.stats['integral'] / ti

                _addplot( h.title, curve, fraction = calcfrac(h) )

                for i in range(1,len(hists)):
                    hi = hists[i]
                    h.add_contents(hi)
                    newcurve = h._hist_curve()
                    assert _np.array_equal(curve[0],newcurve[0])
                    _addplot(hi.title, newcurve, curve,
                             fraction = calcfrac(hi))
                    curve = newcurve

                plt.gca().errorbar(**hist_main.errorbar_args())
                if not logy:
                    plt.ylim(0.0,ymax_non_NOSCAT[0]*1.3 or None)
                plt.xlim(0.0,180.0)
            else:
                if rebin_factor != 1:
                    hist = hist.clone( rebin_factor = rebin_factor )
                hist.plot_hist(plt=plt,do_show = False)
                if not title:
                    plt.title(hist.title or '<untitled histogram>')
                if not logy:
                    plt.ylim(0.0)

            plt.xticks(_np_linspace(0.0,180.0,180//30+1))
            plt.xticks(_np_linspace(0.0,180.0,180//15+1),minor=True)
            plt.xlabel('Exit Angle (degrees)')
            plt.ylabel('Intensity (arbitrary units)')
            if title:
                plt.title(title)

            suptitle_fs = 'medium'
            if len(self.cfgstr)>40:
                suptitle_fs = 'small'
            if len(self.cfgstr)>80:
                suptitle_fs = 'x-small'
            if len(self.cfgstr)>101:
                suptitle_fs = 'xx-small'
            plt.suptitle(self.cfgstr,fontsize=suptitle_fs)

            absfracstr = 'unknown fraction'
            if tot_integral_with_absorption and tot_integral:
                _absfrac =  1.0 - tot_integral/tot_integral_with_absorption
                absfracstr = _fractionval_fmt(_absfrac)
            plt.plot([], [], ' ', label="Absorbed %s"%absfracstr)

            legargs = {}

            _plt_final(do_grid = do_grid,
                       do_legend = do_legend,
                       do_show = do_show,
                       logx = False,
                       logy = logy,
                       plt = plt,
                       extra_legend_kwargs = legargs )

    return Results( main_hist_content = ct,
                    main_hist_errsq = errsq,
                    json_details = res_json,
                    cfgstr = cfgstr,
                    setup_info = setup_info)

_length_units = {'km':1000.0,
                 'm':1.0,
                 'meter':1.0,
                 'cm':0.01,
                 'mm':0.001,
                 'mfp': None,#special
                 'nm':1e-9,
                 'aa':1e-10,
                 'Aa':1e-10,
                 'AA':1e-10,
                 'angstrom':1e-10}

_energy_units = {'eV':1.0,
                 'keV':1e3,
                 'MeV':1e6,
                 'GeV':1e9,
                 'meV':0.001,
                 'neV':1e-9,
                 'aa':None,
                 'Aa':None,
                 'AA':None,
                 'angstrom':None}

def _tofloat(s):
    try:
        return float(s)
    except ValueError:
        return None

def _parse_unit(valstr,unitmap):
    v = _tofloat(valstr)
    if v is not None:
        return v, None, None
    valstr=valstr.strip()
    for unit,unitvalue in sorted(unitmap.items(),key=lambda x : (-len(x),x)):
        if valstr.endswith(unit):
            v = _tofloat(valstr[:-len(unit)])
            if v is not None:
                return v, unit, unitvalue
    return None,None,None

def _parse_energy( valstr ):
    v,u,uv = _parse_unit( valstr, _energy_units )
    if v is not None and u is None and uv is None:
        from .exceptions import NCBadInput
        raise NCBadInput('Invalid energy specification (missing unit'
                         f' like Aa or eV): "{valstr}"')
    if v is None:
        from .exceptions import NCBadInput
        raise NCBadInput(f'Invalid energy specification: "{valstr}"')
    if u is not None and u.lower() in ('aa','angstrom'):
        from .constants import wl2ekin
        return wl2ekin(v)
    v *= uv
    return v

def _parse_length( valstr, mfp = None ):
    v,u,uv = _parse_unit( valstr, _length_units )
    if v is not None and u is None and uv is None:
        from .exceptions import NCBadInput
        _ex0="mfp" if mfp is not None else "mm"
        raise NCBadInput('Invalid length specification (missing unit like '
                         f'{_ex0} or cm): "{valstr}"')
    if v is None:
        from .exceptions import NCBadInput
        raise NCBadInput(f'Invalid length specification: "{valstr}"')
    if u=='mfp':
        if mfp is None:
            raise ValueError('Invalid length specification ("mfp" '
                             f'not supported for this parameter): "{valstr}"')
        v *= mfp
    else:
        v *= uv
    return v

def _macroxs_if_isotropic( mat, **xsect_kwargs ):
    #macroxs_scatter in units of [1/m]
    return ( None if mat.scatter.isOriented()
             else ( mat.info.factor_macroscopic_xs
                    * mat.scatter.xsect( **xsect_kwargs )
                    / _parse_length('1cm') ) )

def quick_diffraction_pattern_autoparams( cfgstr ):
    from .core import load as ncload
    mat = ncload( cfgstr )

    neutron_wl = 1.8 # todo: depend on e.g. Bragg threshold?

    unit_mm = _parse_length('1mm')
    unit_cm = _parse_length('1cm')
    unit_m = _parse_length('1m')
    assert unit_m == 1.0
    macroxs_scatter = _macroxs_if_isotropic( mat, wl=neutron_wl )
    if not macroxs_scatter:
        material_thickness = '1cm'
    else:
        mfp_scatter = 1.0 / macroxs_scatter
        def _round2digits(x):
            return int(x*100+0.5) * 0.01
        if mfp_scatter <= unit_cm:
            material_thickness = f'{_round2digits(mfp_scatter/unit_mm):g}mm'
        elif mfp_scatter <= unit_m:
            material_thickness = f'{_round2digits(mfp_scatter/unit_cm):g}cm'
        else:
            material_thickness = f'{_round2digits(mfp_scatter/unit_m):g}m'
    return dict( neutron_energy_str = f'{neutron_wl}Aa',
                 material_thickness_str = material_thickness )

def quick_diffraction_pattern( cfgstr, *,
                               neutron_energy,
                               material_thickness,
                               nstat = 'auto',
                               nthreads = 'auto' ):

    ### TODO: change c-api and use:
    ###    from .misc import MaterialSource
    ###    matsrc = MaterialSource(material)
    ###    if matsrc.is_preloaded():
    ###        #from .exceptions import NCBadInput
    ###        raise NCBadInput( 'Diffraction patterns can not be produced for'
    ###                          ' preloaded materials (for instance simply passing'
    ###                          ' a cfgstring will work).' )
    ###

    neutron_energy_eV = _parse_energy( neutron_energy )

    from .core import load as ncload
    mat = ncload( cfgstr )
    #unit_cm = _parse_length('1cm')
    unit_m = _parse_length('1m')

    macroxs_scatter = _macroxs_if_isotropic( mat, ekin=neutron_energy_eV )
    mfp_scatter = 1/macroxs_scatter if macroxs_scatter else float('inf')

    sphere_diameter = _parse_length(material_thickness, mfp=mfp_scatter)

    def simfct( n, cfgstr ):
        import time
        t0 = time.time()
        r_meter = 0.5*sphere_diameter / unit_m
        res = runsim_diffraction_pattern( cfgstr,
                                          geomcfg = f'sphere;r={r_meter}',
                                          srccfg = (f'constant;ekin={neutron_energy_eV}'
                                                    f';n={n};z={-r_meter*(1-1e-13)}'),
                                          tally_detail_lvl = 2,
                                          nthreads = nthreads
                                         )
        t1 = time.time()
        return t1-t0, res

    if nstat is None or nstat=='auto':
        #simfct(1,cfgstr)#build mat cache
        for nstat in [1e4,1e5,1e6,1e7]:
            t,res = simfct(nstat,cfgstr)
            #Usually, end within a second in total, but in worst cases, up to
            #10seconds:
            if t>0.1 and nstat >= 1e6:
                break
            if t>1.0:
                break
    else:
        t,res=simfct(nstat,cfgstr)

    res.setup_info['nstat'] = nstat

    return res
