
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

from ._numpy import _np, _np_linspace
from .exceptions import NCBadInput

_breakdown_colors = {
    'SINGLESCAT_ELAS' : 'blue',
    'SINGLESCAT_INELAS' : 'orangered',
    'MULTISCAT_PUREELAS' : 'cornflowerblue',
    'MULTISCAT_OTHER' : 'orange',
    'NOSCAT' : 'green',
}

_cache_tally_info = [None]
def tally_info():
    if _cache_tally_info[0] is not None:
        return _cache_tally_info[0]
    from .misc import evaluate_query as evalquery
    _cache_tally_info[0] =  evalquery(['mmc','cfgdoc','tally'], readonly = True)
    return _cache_tally_info[0]

def run( *, resclass, unpack,
         cfgstr, geomcfg, srccfg, scenario, enginecfg,
         callback, callback_options ):
    def check( v, sn):
        if not ( v is None or isinstance(v,str) ):
            raise NCBadInput(f'The {sn} parameter must be a string.')
    if cfgstr is None:
            raise NCBadInput('Missing required parameter: cfgstr.')
    enginecfg = '' if enginecfg is None else enginecfg
    check(cfgstr,'cfgstr')
    check(geomcfg,'geomcfg')
    check(srccfg,'srccfg')
    check(scenario,'scenario')
    check(enginecfg,'engiornecfg')
    check(callback_options,'callback_options')
    if unpack not in ('dict', 'json', 'dict_jsoncompat', 'object'):
        raise NCBadInput('Invalid value of unpack (must be "dict",'
                         ' "json", "dict_jsoncompat", or "object"):'
                         f' {repr(unpack)}')
    query = ['mmc','run', cfgstr]#, geomcfg, srccfg, enginecfg]
    n_geomsrc = ( ( 1 if geomcfg is not None else 0 )
                  + ( 1 if srccfg is not None else 0 ) )
    if n_geomsrc == 0 and scenario is None:
        raise NCBadInput('Missing required parameters for geometry and source'
                         '. Please supply either a scenario string,'
                         ' or both of geomcfg + srccfg strings.')
    if n_geomsrc > 0 and scenario is not None:
        raise NCBadInput('Inconsistent parameters. Do not supply geomcfg or'
                         ' srccfg when also supplying a scenario string.')
    if n_geomsrc == 2 and scenario is None:
        query += [ geomcfg, srccfg ]
    elif n_geomsrc == 0 and scenario is not None:
        query += [ scenario ]
    else:
        #Should have been caught above, but just as a safety we throw also here:
        raise NCBadInput('Inconsistent parameters.')
    query += [ enginecfg ]
    if callback:
        from ._chooks import _get_raw_cfcts
        _rawfct = _get_raw_cfcts()
        res = _rawfct['flexmmcrun']( query, callback, callback_options )
    else:
        if callback_options is not None:
            raise NCBadInput('Inconsistent parameters. Do not supply'
                             ' callback_options without a callback function.')
        from .misc import evaluate_query
        res = evaluate_query( query, unpack = False )

    if unpack == 'json':
        return res
    import json
    res = json.loads(res)
    if unpack == 'dict_jsoncompat':
        return res
    from .hist import Hist1D
    res = Hist1D.objectify_data(res)
    if unpack == 'dict':
        return res
    return resclass( res )

def results_check_compat_impl( _self, other, threshold, errfct ):
    if _self is other:
        return
    #Allow things like seed, nthreads, and size of statistics to fluctuate:
    volatile=[
        ('engine','cfgstr'),
        ('engine','decoded','seed'),
        ('engine','decoded','nthreads'),
        ('src','cfgstr'),
        ('src','decoded','n'),
    ]
    volatile = set(volatile)
    def cmp(d1,d2,keylist):
        k1=d1.keys()
        if k1!=d2.keys():
            return False
        for k in k1:
            keylist.append(k)
            block = tuple(keylist) in volatile
            if block:
                keylist.pop()
                continue
            v1 = d1[k]
            v2 = d2[k]
            if isinstance(v1,dict) and isinstance(v2,dict):
                if not cmp(v1,v2,keylist):
                    return False
            elif v1 != v2:
                return False
            keylist.pop()
        return True
    keylist=[]
    if not cmp( _self._raw_data()['input'],
                other._raw_data()['input'],
                keylist ):
        k='/'.join(keylist)
        return errfct(f'incompatible input values for "{k}"')

    tallies_s = dict( (t.name, t) for t in _self.tallies )
    tallies_o = dict( (t.name, t) for t in other.tallies )
    tally_names = sorted(tallies_s.keys())

    if ( tally_names != sorted(tallies_o.keys())
         or len(tally_names)!=len(tallies_s)
         or len(tally_names)!=len(tallies_o) ):
        #should already have been caught unless data format changed
        return errfct('incompatible tally list')

    nhists=sum( t._nhists for t in tallies_s.values() )
    actual_threshold = ( 0.0 if threshold is None else threshold ) / nhists

    if actual_threshold <= 0.0:
        return

    for tn in tally_names:
        t_s = tallies_s.get( tn )
        t_o = tallies_o.get( tn )
        if not t_s or not t_o:
            return errfct(f'unexpectedly tally missing: {tn}')
        hists_s = t_s.histograms
        hists_o = t_o.histograms
        hk = hists_s.keys()
        if hists_o.keys() != hk:
            return errfct(f'incompatible histograms available for tally {tn}')
        for hname in hk:
            h_s = hists_s[hname]
            h_o = hists_o[hname]
            if not h_s.check_compat( h_o,
                                     threshold = actual_threshold,
                                     force_norm = True ):
                return (f'incompatible histograms: {tn}/{hk}')
    return

def _validate_mmcresults_dict(data):
    #Very brief high-level validation:
    if not data.get('datatype') == 'NCrystalMiniMCResults_v1':
        return False
    return set(data.keys()) == set(['datatype','input','output'])

def _determine_rebin_factor( current_nbins,
                             max_nbins = None,
                             rebin_factor = None ):
    assert current_nbins >= 1
    if max_nbins is None and rebin_factor is None:
        return 1
    if max_nbins is not None and rebin_factor is not None:
        raise NCBadInput('Can not set both max_nbins and rebin_factor.')
    if rebin_factor is not None:
        if current_nbins%rebin_factor != 0:
            raise NCBadInput(f'Invalid rebin factor {rebin_factor} '
                             f'is not a divisor of nbins={current_nbins}.')
        return rebin_factor
    assert max_nbins is not None
    assert max_nbins >= 1
    if max_nbins >= current_nbins:
        return 1
    for n in list(range(max_nbins,0,-1)):
        if current_nbins % n == 0:
            return current_nbins//n
    return 1

def _plot_tally( minimcresults_dict, tallyname,
                 do_legend, breakdown, max_nbins, rebin_factor,
                 do_show, do_newfig, do_grid, logy, title, plt, axis ):

    from .plot import _import_matplotlib_plt, _plt_final
    from .hist import Hist1D
    assert isinstance(minimcresults_dict,dict)
    assert minimcresults_dict.get('datatype') == 'NCrystalMiniMCResults_v1'
    assert 'input' in minimcresults_dict
    assert 'output' in minimcresults_dict

    avail_tallies = sorted(minimcresults_dict['output']['tally'].keys())
    if not avail_tallies:
        raise NCBadInput('MiniMC results has no tallies!')

    do_title = not ( not title or title=='none' )

    if tallyname not in avail_tallies:
        raise NCBadInput('Requested tally "%s" not available.'%tallyname)


    tally_dict = minimcresults_dict['output']['tally'][tallyname]
    assert isinstance(tally_dict,dict)
    assert 'total' in tally_dict
    assert 'tallyname' in tally_dict
    assert tally_dict.get('datatype')=='NCrystalMiniMCTallyHistBreakdown_v1'
    breakdown_dict = tally_dict.get('breakdown') or None
    if breakdown == 'auto':
        breakdown = bool(breakdown_dict)
    if breakdown and not breakdown_dict:
        raise NCBadInput('breakdown=True but provided tally does not'
                         ' have breakdown histograms')
    def ensure_hist( e ):
        if not isinstance(e,Hist1D):
            e = Hist1D.objectify_data(e)
        if not isinstance(e,Hist1D):
            raise NCBadInput('Invalid input data: not a histogram')
        return e
    if breakdown:
        breakdown =  dict( (k,ensure_hist(v))
                           for k,v in breakdown_dict.items() )
    else:
        breakdown = None
    mainhist = ensure_hist(tally_dict['total'])

    if not breakdown:
        do_legend = False

    if breakdown:
        nbins = [h.nbins for h in breakdown.values()]
        nbins += [mainhist.nbins]
        if len(set(nbins))!=1:
            raise NCBadInput('Inconsistent binning of breakdown histograms')
        nbins = nbins[0]
        _ = set(breakdown.keys()) - set(_breakdown_colors.keys())
        if _:
            raise NCBadInput('Unexpected breakdown histogram key: %s'%_.pop())
    else:
        nbins = mainhist.nbins
    rebin_factor = _determine_rebin_factor( nbins,
                                            max_nbins = max_nbins,
                                            rebin_factor = rebin_factor )
    if rebin_factor != 1:
        mainhist = mainhist.clone(rebin_factor=rebin_factor)
        if breakdown:
            breakdown = dict( (k,v.clone(rebin_factor=rebin_factor))
                              for k,v in breakdown.items() )

    assert bool(axis) == (axis is not None)
    if do_newfig == 'auto':
        do_newfig = not axis
    if do_newfig and axis:
        raise NCBadInput('incompatible: axis is not '
                         'None and do_newfig=True')
    if not plt and not axis:
        from .plot import _import_matplotlib_plt
        plt = _import_matplotlib_plt()
    if do_newfig:
        if not plt:
            raise NCBadInput('unsupported combination: '
                             'axis=None, do_newfig=True, plt=None')
        plt.figure()
        assert axis is None
        axis = plt.gca()
    if not axis:
        assert plt is not None
        axis = plt.gca()

    if not mainhist.integral:
        from ._common import warn
        warn('Aborting plotting of empty histogram!')
        return plt

    #Look at config to determine total weights:
    out_md = minimcresults_dict['output']['metadata']
    tot_weight_incoming = out_md['provided']['weight']
    nrays_incoming = out_md['provided']['count']
    eopts = minimcresults_dict['input']['engine']
    if eopts['decoded']['ignoremiss']:
        tot_weight_incoming -= out_md['miss']['weight']
        nrays_incoming -= out_md['miss']['count']

    def find_nonabsfrac():
        #TODO: move logic to MMCResults method
        if eopts['decoded']['absorption']:
            return mainhist.integral/tot_weight_incoming
        #no absorption, must be 0. We simply verify that the result is
        #consistent with that.
        totw, totw_err = mainhist.integrate_bins()
        nonabsfrac = totw/tot_weight_incoming
        nonabsfrac_err = totw_err/tot_weight_incoming
        assert (nonabsfrac-1.0) < 0.2
        assert (nonabsfrac-1.0) < 5.0*nonabsfrac_err
        return 1.0

    nonabsfrac = find_nonabsfrac()
    absfrac = 1.0 - nonabsfrac

    def _fractionval_fmt(x):
        return f'({x*100.0:.3g}%)'

    label_order = []

    if breakdown:
        hists = [ (k,v) for k,v in breakdown.items() if v.integral > 0.0 ]
        assert hists
        #Sort histograms for (possibly) more meaningful plot order:
        hists.sort( key = lambda kv: ['MULTISCAT_OTHER',
                                      'SINGLESCAT_INELAS',
                                      'MULTISCAT_PUREELAS',
                                      'SINGLESCAT_ELAS',
                                      'NOSCAT'].index(kv[0]))

        h = hists[0][1].clone()
        curve = h._hist_curve()
        ymax_non_NOSCAT = [0.0]
        def _addplot( ptitle, curve, refcurve = None, fraction = None):
            maxval = 0.0
            if refcurve is None:
                y2 = _np.zeros(len(curve[1]))
            else:
                y2=refcurve[1]
                maxval = y2.max()
            if ptitle!='NOSCAT':
                ymax_non_NOSCAT[0] = max(ymax_non_NOSCAT[0],
                                         maxval,
                                         curve[1].max())
            lbl = ptitle if ptitle!='NOSCAT' else 'Transmitted'
            if fraction is not None:
                lbl = f'{lbl} {_fractionval_fmt(fraction)}'
            label_order.append(lbl)
            axis.fill_between(*curve,y2,
                              label=lbl,
                              edgecolor="none",
                              facecolor=_breakdown_colors[ptitle])

        def calcfrac(hh):
            return hh.integral / tot_weight_incoming

        _addplot( hists[0][0], curve, fraction = calcfrac(h) )

        for i in range(1,len(hists)):
            ki, hi = hists[i]
            h.add_contents(hi)
            newcurve = h._hist_curve()
            assert _np.array_equal(curve[0],newcurve[0])
            _addplot(ki, newcurve, curve,
                     fraction = calcfrac(hi))
            curve = newcurve

        axis.errorbar(**mainhist.errorbar_args())
        axis.set_xlim(mainhist.xmin,mainhist.xmax)
    else:
        lbl = 'All outgoing %s'%_fractionval_fmt(nonabsfrac)
        label_order.append(lbl)
        mainhist.plot_hist( plt=plt, do_show = False, label=lbl)
        if do_title and not title and mainhist.title:
            axis.set_title(mainhist.title)
        if not logy:
            axis.set_ylim(0.0)

    from ._mmc_impl import tally_info
    t_info = tally_info()['tallyhistinfo'][tallyname]
    xlbl = t_info['short_descr'].capitalize()
    if t_info['unit']:
        xlbl += ' (%s)'%t_info['unit']
    axis.set_xlabel(xlbl)

    if tallyname=='theta' and mainhist.xmin==0 and mainhist.xmax==180:
        axis.set_xticks(_np_linspace(0.0,180.0,180//30+1))
        axis.set_xticks(_np_linspace(0.0,180.0,180//15+1),minor=True)
    axis.set_ylabel('Intensity (arbitrary units)')
    if do_title:
        axis.set_title(title)

    if absfrac > 0.0:
        lbl="Absorbed %s"%_fractionval_fmt(absfrac)
        label_order.append(lbl)
        axis.plot([], [], ' ', label=lbl)

    if do_legend:
        #Enforce ordering!

        #Instead of handles, labels = axis.get_legend_handles_labels(), we do
        #the following because of our unit tests:
        is_fake_plt = hasattr(axis,'the_real_inspected_object')
        if is_fake_plt:
            #Not real calls, just fake it.
            class FakeHandle:
                def __repr__(self):
                    return 'FakeHandle'
                pass
            handles = [FakeHandle() for lbl in label_order]
            labels = [lbl for lbl in label_order]
        else:
            handles,labels = axis.get_legend_handles_labels()

        _ = sorted([hl for hl in zip(handles,labels)],
                   key = lambda hl : label_order.index(hl[1]))
        handles, labels = [hl[0] for hl in _], [hl[1] for hl in _]
        plt.legend(handles,labels)

    return _plt_final(do_grid = do_grid,
                      do_legend = False,#plt.legend was already called above
                      do_show = do_show,
                      logx = False,
                      logy = logy,
                      plt = plt )
