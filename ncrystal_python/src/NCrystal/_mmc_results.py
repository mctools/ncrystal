
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

class MMCDiffractionResults:

    def __new__(cls, *args, **kwargs):
        raise TypeError("Use named constructors instead of direct instantiation")

    @classmethod
    def from_dict( cls, data ):
        from ._hist import Hist1D
        import copy
        assert isinstance( data, dict )
        o = super().__new__(cls)
        assert set(data.keys()) == set(['setup_info','hists','cfgstr'])
        o.__hists = [Hist1D(e) for e in data['hists']]
        o.__setup_info = copy.deepcopy(data['setup_info'])
        o.__cfgstr = str(data['cfgstr'])
        return o

    def nbins( self ):
        nb = set([h.nbins for h in self.__hists])
        assert len(nb)==1, "Inconsistent binnings"
        return nb.pop()

    @classmethod
    def from_json( cls, data ):
        import json
        return cls.from_dict( json.loads(data) )

    def to_dict( self, json_compat = False ):
        return dict(
            setup_info = self.__setup_info,
            hists = [ h.to_dict( json_compat = json_compat)
                      for h in self.__hists ],
            cfgstr = self.__cfgstr,
        )

    def to_json( self ):
        import json
        return json.dumps(self.to_dict(json_compat=True))

    @classmethod
    def _from_C( cls,
                 main_hist_content,
                 main_hist_errsq,
                 json_details,
                 cfgstr,
                 setup_info ):
        _ensure_numpy()
        o = super().__new__(cls)
        o.__setup_info = setup_info

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
        o.__hists = hists
        o.__cfgstr = cfgstr
        return o


    def rebin( self, rebin_factor ):
        """Reduce granularity of binnings in histogram by provided factor, which
        must be a divisor in the current number of bins. Returns self."""
        nbins = self.nbins()
        if not rebin_factor>=1 or not nbins % rebin_factor == 0:
            raise NCBadInput( f'Can not rebin nbins={nbins} with '
                              f'rebin_factor={rebin_factor} (rebin_factor'
                              ' must be a positive divisor of nbins).' )
        if rebin_factor != 1:
            for h in self.__hists:
                h.rebin(rebin_factor)
        return self

    def clone( self, rebin_factor = 1 ):
        """Return clone of object, while possibly rebinning all histograms."""
        import copy
        c = super().__new__(MMCDiffractionResults)
        c.__hists = [h.clone(rebin_factor) for h in self.__hists]
        c.__setup_info = copy.deepcopy( self.__setup_info )
        c.__cfgstr = copy.deepcopy( self.__cfgstr )
        return c

    @property
    def histograms( self ):
        return self.__hists

    @property
    def histogram_main( self ):
        """Return main histogram with totals"""
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
        print
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
    def histogram_titles( self ):
        return [h.title for h in self.__hists]

    def histogram( self, title ):
        hh = [ h for h in self.__hists if h.title == title ]
        if not hh:
            raise NCBadInput(f'No histogram found with title {title}'
                             f' (possibilities are {self.histogram_titles}).')
        assert len(hh)==1, "multiple histograms with same title"
        return hh[0]

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
            #axis = plt.gca()
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
                plt.fill_between(*curve,y2,
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

            plt.errorbar(**hist_main.errorbar_args())
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
#FIXME: Option to check_compat on Results objects, and thus check all components.
