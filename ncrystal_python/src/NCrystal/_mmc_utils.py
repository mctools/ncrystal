
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

from ._numpy import _np, _ensure_numpy, _np_linspace
from .exceptions import NCBadInput

_breakdown_colors = {
    'SINGLESCAT_ELAS' : 'blue',
    'SINGLESCAT_INELAS' : 'orangered',
    'MULTISCAT_PUREELAS' : 'cornflowerblue',
    'MULTISCAT_OTHER' : 'orange',
    'NOSCAT' : 'green',
}

#FIXME: Check all doc-strings in entire file, and possibly move
#       MMCResults/MMCTallyView objects to minimc.y


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
        return self.__data.get('breakdown') or None

    @property
    def histograms( self ):
        """Returns histograms as a dictionary of (key,hist). The keys are either
        'total' or one of the breakdown component names (if available).
        """
        d=dict(total=self.hist_total)
        hd = self.hist_breakdown
        if hd is not None:
            assert 'total' not in hd
            d.update(hd)
        return d

    def __eq__(self, o ):
        if id(self) == id(o):
            return True
        if self.name != o.name:
            return False
        if not ( self.hist_total == o.hist_total ):
            return False
        hd = self.hist_breakdown
        ohd = o.hist_breakdown
        if (hd is None) != (ohd is None):
            return False
        if hd is None:
            return True
        if len(hd) != len(ohd):
            return False
        for k,h in hd.items():
            oh = ohd.get(k)
            if oh is None:
                return False
            if ( h != oh ):
                return False
        return True

    @property
    def _nhists( self ):
        return 1 + (len(self.hist_breakdown) if self.hist_breakdown else 0)

    def histogram_sum( self, *, select=None, exclude=None ):
        if isinstance(exclude,str):
            exclude=[exclude]
        if isinstance(select,str):
            select=[select]
        if not exclude and not select:
            return self.hist_total
        histmap = self.hist_breakdown
        if select:
            histmap = dict( (hn,h) for hn,h in histmap.items()
                            if hn in select )
        if exclude:
            histmap = dict( (hn,h) for hn,h in histmap.items()
                            if hn not in exclude )
        hl = [h for hn,h in sorted(histmap.items()) ]
        if len(hl) <= 1:
            return hl[0] if hl else None
        h = hl[0].clone()
        for o in hl[1:]:
            h.add_contents( o )
        return h

    def unit( self ):
        from ._mmc_impl import tally_info
        return tally_info()['tallyhistinfo'][self.name]['unit']

    def short_description( self ):
        from ._mmc_impl import tally_info
        return tally_info()['tallyhistinfo'][self.name]['short_descr']

    def dump( self, *args, **kwargs):
        return self.hist_total.dump(*args,**kwargs)

    def plot( self,
              breakdown = 'auto',
              max_nbins = None,
              rebin_factor = None,
              do_show = True,
              do_newfig = 'auto',
              do_grid = False,
              do_legend = 'auto',
              logy = True,
              title = None,
              plt = None,
              axis = None ):
        """todo: more here.

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
                     do_legend = do_legend,
                     logy = logy,
                     title = title,
                     plt = plt,
                     axis = axis )

class MMCResults:

    def __init__(self, data):
        _ensure_numpy()
        import copy
        needsval = True
        if isinstance(data,MMCResults):
            needsval = False
            data = copy.deepcopy(data.__data)
        elif ( isinstance( data, bytes )
               or isinstance( data, str )
               or hasattr( data, '__fspath__' ) ):
            from ._common import flex_load_json
            from .hist import Hist1D
            data = Hist1D.objectify_data( flex_load_json( data ) )
        elif isinstance(data,dict):
            from .hist import Hist1D
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
    def setup( self ):
        return self.__data['input']

    @property
    def output_metadata( self ):
        return self.__data['output']['metadata']

    def short_title( self, latex = False ):
        n = self.output_metadata['provided']['count']
        s = self.setup['src']['metadata']['energy_description']
        g = self.setup['geom']['decoded']['short_description']
        if latex:
            from ._common import _latex_format
            n = '$%s$'%_latex_format(n)
        return f'{n} {s} neutrons through {g}'

    def long_title( self, latex = False ):
        #NB: latex parameter currently does nothing
        s = self.setup['src']['cfgstr']
        g = self.setup['geom']['cfgstr']
        e = self.setup['engine']['cfgstr']
        r = f'"{s}" on "{g}"'
        n = self.output_metadata['tallied']['count']
        r = f'{r} ({n} fills)'
        return f'{r} ("{e}")' if e else r

    def dump( self, do_print = True, prefix = '',
              tally_filter_fct = None ):
        o = []
        o.append('%sNCrystal MiniMC results:'%prefix)
        o.append('  inputs cfg:')
        o.append('    material : "%s"'%self.setup['material']['cfgstr'])
        o.append('    engine   : "%s"'%self.setup['engine']['cfgstr'])
        o.append('    source   : "%s"'%self.setup['src']['cfgstr'])
        o.append('    geometry : "%s"'%self.setup['geom']['cfgstr'])
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
        o.append('    tallied ray count: %s particles (weight sum: %g)'%(
            fmti(outmd['tallied']['count']),outmd['tallied']['weight']))
        for t in self.tallies:
            if tally_filter_fct and not tally_filter_fct(t.name):
                continue
            o.append('    tally "%s":'%t.name)
            o += t.hist_total.dump( prefix = '      ',
                                    contents = False,
                                    do_print = False ).splitlines()
        #finish up:
        o.append('')
        o = ('\n%s'%prefix).join(o)
        if do_print:
            from ._common import print as ncprint
            ncprint(o)
        return o

    def check_compat( self, other, threshold = 0.05, check=False ):
        """Checks compatiblity between this and other MMCResults object. Two
        objects are deemed compatible if they contain the same input settings,
        and contents of output tally histograms are all compatible (if
        threshold>0).

        The threshold parameter is the nominal p-value threshold required to be
        exceeded in the histogram comparisons to achieve compatibility. However,
        to reduce false positives, note that the individual histogram
        comparisons will merely have to exceed a reduced p-value of threshold/N,
        where N is the total number of histograms in the MMCResults object.

        Returns True in case of compatibility, False otherwise.

        If check is True, a CalcError exception is raised in case of
        incompatibilities.
        """

        from ._mmc_impl import results_check_compat_impl
        assert isinstance(other,MMCResults)
        assert isinstance(threshold,float)
        if check:
            def errfct( errmsg ):
                from .exceptions import NCCalcError
                raise NCCalcError(f'Incompatible MMCResults ({errmsg})')
        else:
            def errfct( errmsg ):
                return errmsg
        errmsg = results_check_compat_impl(self,other,threshold,errfct)
        if errmsg:
            return False
        return True
