
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

from ._numpy import _np, _ensure_numpy, _np_linspace
from .exceptions import NCBadInput

_breakdown_colors = {
    'SINGLESCAT_ELAS' : 'blue',
    'SINGLESCAT_INELAS' : 'orangered',
    'MULTISCAT_PUREELAS' : 'cornflowerblue',
    'MULTISCAT_OTHER' : 'orange',
    'NOSCAT' : 'green',
}

def validate_mmcresults_dict(data):
    #FIXME!!!
    #Check that histograms have been objectified.
    pass


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

def _plot_tally( minimcresults_dict,
                 tallyname,
                 breakdown = 'auto',
                 max_nbins = None,
                 rebin_factor = None,
                 do_show = True,
                 do_newfig = 'auto',
                 do_grid = False,
                 logy = True,
                 title = None,
                 plt = None,
                 axis = None ):

    from .plot import _import_matplotlib_plt, _plt_final
    from ._hist import Hist1D
    assert isinstance(minimcresults_dict,dict)
    assert minimcresults_dict.get('datatype') == 'NCrystalMiniMCResults_v1'
    assert 'input' in minimcresults_dict
    assert 'output' in minimcresults_dict

    avail_tallies = sorted(minimcresults_dict['output']['tally'].keys())
    if not avail_tallies:
        raise NCBadInput('MiniMC results has no tallies!')

    do_title = not ( not title or title=='none' )


    #if tallyname == 'auto':
    #    if len(avail_tallies)==1:
    #        tallyname = avail_tallies[0]
    #    else:
    #        raise NCBadInput('MiniMC provided more than one tally. Please '
    #                         'select one with the tallyname keyword (available'
    #                         ' tallies: "%s")'%(", ".join(avail_tallies)))
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
    if minimcresults_dict['input']['engine']['decoded']['ignoremiss']:
        tot_weight_incoming -= out_md['miss']['weight']
        nrays_incoming -= out_md['miss']['count']
    nonabsfrac = mainhist.integral/tot_weight_incoming
    absfrac = 1.0 - nonabsfrac

    #FIXME: If enginecfg has absorption=False, we should simply verify that
    #absfrac is super small (i.e. roulette fluctuations), taking statistics into
    #account, and then enforce sum of remaining parts to be 100%. It could be
    #that we should move this logic onto the MMCResults object.

    def _fractionval_fmt(x):
        return f'({x*100.0:.3g}%)'

    label_order = []

    do_legend = False
    if breakdown:
        do_legend = True
        #hists = [h for h in self.histograms if h.title!='MAIN']
        #hist_main = [h for h in self.histograms if h.title=='MAIN'][0]
        #if not hists:
        #    from .exceptions import NCCalcError
        #    raise NCCalcError('Can not produce breakdown plots'
        #                      ' without more detailed tally.')
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
        #if not logy:
        #    axis.set_ylim(0.0,ymax_non_NOSCAT[0]*1.3 or None)#Fixme: hiding transmitted!
        axis.set_xlim(mainhist.xmin,mainhist.xmax)
    else:
        lbl = 'All outgoing %s'%_fractionval_fmt(nonabsfrac)
        label_order.append(lbl)
        mainhist.plot_hist(plt=plt,do_show = False,
                           label=lbl)
        if do_title and not title and mainhist.title:
            axis.set_title(mainhist.title)
        if not logy:
            axis.set_ylim(0.0)

    if tallyname=='mu':
        axis.set_xticks(_np_linspace(0.0,180.0,180//30+1))
        axis.set_xticks(_np_linspace(0.0,180.0,180//15+1),minor=True)
        #FIXME: xlabel needs to be embedded in Tally?
        axis.set_xlabel('Exit Angle (degrees)')
    else:
        axis.set_xlabel(tallyname)
    axis.set_ylabel('Intensity (arbitrary units)')
    if do_title:
        axis.set_title(title)

    #FIXME suptitle_fs = 'medium'
    #FIXME if len(self.cfgstr)>40:
    #FIXME     suptitle_fs = 'small'
    #FIXME if len(self.cfgstr)>80:
    #FIXME     suptitle_fs = 'x-small'
    #FIXME if len(self.cfgstr)>101:
    #FIXME     suptitle_fs = 'xx-small'
    #FIXME figure = axis.get_figure()
    #FIXME if figure:
    #FIXME     figure.suptitle(self.cfgstr,fontsize=suptitle_fs)

    if absfrac > 0.0:
        lbl="Absorbed %s"%_fractionval_fmt(absfrac)
        label_order.append(lbl)
        axis.plot([], [], ' ', label=lbl)

    #legargs = {}#fixme
    do_legend=True#fixme
    if do_legend:
        #Enforce ordering!
        handles, labels = axis.get_legend_handles_labels()
        _ = sorted([hl for hl in zip(handles,labels)],
                   key = lambda hl : label_order.index(hl[1]))
        handles, labels = [hl[0] for hl in _], [hl[1] for hl in _]
        plt.legend(handles,labels)

    _plt_final(do_grid = do_grid,
               do_legend = False,
               do_show = do_show,
               logx = False,
               logy = logy,
               plt = plt )


#FIXME: Missing doc-strings
#FIXME: Option to check_compat on Results objects, and thus check all comps?

class MMCTallyView:

    """Convenience class to handle a given MiniMC tally. These objects should
    not be constructed manually by users, but rather created from an MMCResult
    object's .tally() method or .tallies property."""

    def __new__(cls, *args, **kwargs):
        raise TypeError("Do not create MMCTallyView objects directly")

    @classmethod
    def _internal_create( cls, mmcresults, td ):
        assert td['datatype']=='NCrystalMiniMCTallyHistBreakdown_v1'
        o = super().__new__(cls)
        #Keeps ref to MMCResult mother object, so mother should never keep refs
        #to MMCTallyView objects!
        o.__mmcresults = mmcresults
        o.__data = td
        return o

    @property
    def name( self ):
        return self.__data['tallyname']

    @property
    def mother( self ):
        """Return the MMCResults object which this MMCTallyView object belongs
        to."""
        return self.__mmcresults

    @property
    def hist_total( self ):
        """Return main histogram with all tallied events."""
        h = self.__data.get('total')
        assert h
        return h

    @property
    def hist_breakdown( self ):
        """Return histograms broken down by component, as a dictionary of
        (component_name,histogram). Returns None if such breakdown histograms
        were not enabled for the tally.
        """
        #Fixme: Make it possible to "lock" a histogram, making it effectively
        #read-only (it can always be cloned in case mutations are needed.
        return self.__o.get('breakdown') or None

    def plot( self,
              breakdown = 'auto',
              max_nbins = None,
              rebin_factor = None,
              do_show = True,
              do_newfig = 'auto',
              do_grid = False,
              logy = True,
              title = None,
              plt = None,
              axis = None ):
        """Fixme more here.

        If title is None, "auto" or "short", a short title will be auto
        generated. If it is "long", a longer title with full configuration
        strings will be used. If title is "none" or False, no title will be
        shown.

        do_newfig='auto' will cause a new figure to be created, unless axis is
        provided.

        returns plt.

        """

        if title in (None,'auto','short'):
            title = self.__mmcresults.short_title(latex=True)
        elif title=='long':
            title = self.__mmcresults.long_title(latex=True)

        _plot_tally( self.__mmcresults._raw_data(),
                     tallyname = self.name,
                     breakdown = breakdown,
                     max_nbins = max_nbins,
                     rebin_factor = rebin_factor,
                     do_show = do_show,
                     do_newfig = do_newfig,
                     do_grid = do_grid,
                     logy = logy,
                     title = title,
                     plt = plt,
                     axis = axis )

class MMCResults:

    #fixme: docstrings (everywhere in file). And move to other files as
    #appropriate.

    def __init__(self, data):
        _ensure_numpy()
        import copy
        needsval = True
        if isinstance(data,MMCResults):
            needsval = False
            data = copy.deepcopy(data.__data)
        elif isinstance( data, bytes ) or isinstance( data, str ):
            from ._common import flex_load_json
            from ._hist import Hist1D
            data = Hist1D.objectify_data( flex_load_json( data ) )
        elif isinstance(data,dict):
            from ._hist import Hist1D
            data = Hist1D.objectify_data(data)
        else:
            raise NCBadInput('Unsupported data format')
        if needsval:
            validate_mmcresults_dict(data)
        self.__data = data

    def to_dict( self, json_compat = False ):
        if json_compat:
            from ._common import copy_and_deobjectify_data
            return copy_and_deobjectify_data( self.__data )
        else:
            import copy
            return copy.deepcopy( self.__data )

    def to_json( self ):
        import json
        return json.dumps(self.to_dict(json_compat=True))

    def _raw_data( self ):
        return self.__data

    def tally( self, tallyname ):
        o = self.__data['output']
        t = o['tally'].get(tallyname)
        if t is None:
            msg=f'Tally not available in MiniMC dataset: "{tallyname}"'
            tn = self.tally_names
            if not tn:
                msg += ' (not tallies were enabled!).'
            else:
                msg += ' (available tallies are "%s")'%('", "'.join(tn))
            raise NCBadInput(msg)
        return MMCTallyView._internal_create( self, t )

    @property
    def tallies( self ):
        o = self.__data['output']
        t = o['tally']
        return [ MMCTallyView._internal_create( self, t[tn] )
                 for tn in self.tally_names ]

    @property
    def tally_names( self ):
        return sorted(self.__data['output']['tally'].keys())

    @property
    def src_miss_counts( self ):
        return self.__data['output']['metadata']['miss']['count']

    @property
    def src_miss_totalweight( self ):
        return self.__data['output']['metadata']['miss']['count']

    @property
    def src_provided_counts( self ):
        return self.__data['output']['metadata']['provided']['count']

    @property
    def src_provided_totalweight( self ):
        return self.__data['output']['metadata']['provided']['count']

    @property
    def setup_material_cfgstr( self ):
        return self.__data['input']['material']['cfgstr']

    @property
    def setup_material_cfgstr_decoded( self ):
        return self.__data['input']['material']['decoded']

    @property
    def setup_source_cfgstr( self ):
        return self.__data['input']['source']['cfgstr']

    @property
    def setup_source_cfgstr_decoded( self ):
        return self.__data['input']['source']['decoded']

    @property
    def setup_geometry_cfgstr( self ):
        return self.__data['input']['geometry']['cfgstr']

    @property
    def setup_geometry_cfgstr_decoded( self ):
        return self.__data['input']['geometry']['decoded']

    @property
    def setup_engine_cfgstr( self ):
        return self.__data['input']['engine']['cfgstr']

    @property
    def setup_engine_cfgstr_decoded( self ):#fixme engine_cfgstr -> enginecfg, etc.
        return self.__data['input']['engine']['decoded']

    def short_title( self, latex = False ):
        n = self.src_provided_counts#fixme: weights?
        s = self.setup_source_cfgstr_decoded['short_neutron_description']
        g = self.setup_geometry_cfgstr_decoded['short_description']
        if latex:
            from ._common import _latex_format
            n = '$%s$'%_latex_format(n)
        return f'{n} {s} neutrons through {g}'

    def long_title( self, latex = False ):
        #NB: latex parameter currently does nothing.
        s = self.setup_source_cfgstr
        g = self.setup_geometry_cfgstr
        e = self.setup_engine_cfgstr
        r = f'"{s}" on "{g}"'
        return f'{r} ("{e}")' if e else r

    def plot_breakdown( self, *args, **kwargs ):
        """Obsolete function provided for backwards compatibility purposes."""
        if __cache_plotbdwarn[0]:
            from ._common import warn
            __cache_plotbdwarn[0] = False
            warn('The .plot_breakdown(..) method is obsolete and only works '
                 'when there is just a single tallied quantity. Please migrate '
                 'your code to use .tally("mu").plot() instead (assuming "mu" '
                 'is the tallied quantity). Or access all tallies with the '
                 '.tallies property, and invoke .plot() on them individually.')
        tns = self.tally_names
        if len(tns)!=1:
            raise NCBadInput('The .plot_breakdown is obsolete and only works'
                             ' when there is exactly one tallied quantity.')
        self.tally(tns[0]).plot(*args,**kwargs)

    def dump( self, do_print = True, prefix = '' ):
        o = []
        o.append('%sNCrystal MiniMC results:'%prefix)
        o.append('  inputs cfg:')
        o.append('    material : "%s"'%self.setup_material_cfgstr)
        o.append('    engine   : "%s"'%self.setup_engine_cfgstr)
        o.append('    source   : "%s"'%self.setup_source_cfgstr)
        o.append('    geometry : "%s"'%self.setup_geometry_cfgstr)
        o.append('  output:')
        outmd = self.__data['output']['metadata']
        def fmti( x ):
            xs = '%.2g'%x
            return xs if int(float(xs))==x else str(x)
        o.append('    src ray count: %s particles (weight sum: %g)'%(
            fmti(outmd['provided']['count']),outmd['provided']['weight']))
        o.append('    src rays missing geometry: %s particles (weight sum: %g)'%(
            fmti(outmd['miss']['count']),outmd['miss']['weight']))
        f_c = outmd['miss']['count']*100.0/outmd['provided']['count']
        f_w = outmd['miss']['weight']*100.0/outmd['provided']['weight']
        o.append('    src rays miss fraction:'
                 ' %g%% (by count) %g%% (by weight)'%(f_c,f_w))
        for t in self.tallies:
            o.append('    tally "%s":'%t.name)
            o += t.hist_total.dump( prefix = '      ',
                                    contents = False,
                                    do_print=False ).splitlines()
        #finish up:
        o.append('')
        o = ('\n%s'%prefix).join(o)
        if do_print:
            print(o)
        return o

__cache_plotbdwarn = [True]
