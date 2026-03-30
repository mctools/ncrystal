
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

__all__ = ['MMCResults','MMCTallyView']

from ._numpy import _ensure_numpy
from .exceptions import NCBadInput

class MMCResults:

    """Convenience class to handle the results of a given MiniMC run, providing
    access not only to tally results, but also configuration metadata and other
    statistics collected during the run.

    Note that these result objects can be converted to and from dictionaries and
    JSON data. This means for instance that it is possible to store the data in
    a file in JSON format, and load it again later.
    """

    def __init__(self, data):
        """Initialise object from JSON data, a result dictionary, or an existing
        MMCResult object.
        """
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
            from ._mmc_impl import _validate_mmcresults_dict
            if not _validate_mmcresults_dict(data):
                raise NCBadInput('Data seems to be in an unsupported format')
        self.__data = data

    def to_dict( self, json_compat = False ):
        """Serialise as a dictionary. This can optionally be JSON compatible, if
        json_compat is set to true.
        """

        if json_compat:
            from ._common import copy_and_deobjectify_data
            return copy_and_deobjectify_data( self.__data )
        else:
            import copy
            return copy.deepcopy( self.__data )

    def to_json( self ):
        """Serialise as JSON data (technically a string)."""
        import json
        return json.dumps(self.to_dict(json_compat=True))

    def _raw_data( self ):
        return self.__data

    @property
    def tally_names( self ):
        """Get a list of names of all tallied quantities available."""
        return sorted(self.__data['output']['tally'].keys())

    def __repr__( self ):
        return str(self)

    def __str__( self ):
        s = self.setup['src']['cfgstr']
        g = self.setup['geom']['cfgstr']
        e = self.setup['engine']['cfgstr']
        return f'MMCResults(srccfg="{s}",geomcfg="{g}",enginecfg="{e}")'

    @property
    def tallies( self ):
        """List of all available tallies, in the form of MMCTallyView
        objects.
        """
        o = self.__data['output']
        t = o['tally']
        return [ MMCTallyView._internal_create( self, t[tn] )
                 for tn in self.tally_names ]

    def tally( self, tallyname ):
        """Return an MMCTallyView object for a particular named tally (use the
        .tally_names property to see which are available).
        """
        o = self.__data['output']
        t = o['tally'].get(tallyname)
        if t is None:
            msg=f'Tally not available in MiniMC dataset: "{tallyname}"'
            tn = self.tally_names
            if not tn:
                msg += ' (no tallies were enabled!).'
            else:
                msg += ' (available tallies are "%s")'%('", "'.join(tn))
            raise NCBadInput(msg)
        return MMCTallyView._internal_create( self, t )

    @property
    def setup( self ):
        """Return a dictionary with information about the configuration of the
        MiniMC run.
        """
        return self.__data['input']

    @property
    def output_metadata( self ):
        """Return a dictionary with information about the configuration of the
        MiniMC run. This contains summation (counts or weight sums) about
        neutrons provided by the source, neutrons from the source missing the
        geometry, and neutrons entering the tallies.
        """
        return self.__data['output']['metadata']

    def short_title( self, latex = False ):
        """Returns a short string with a very brief description of the
        simulation setup, which is for instance suitable for quick plot titles.

        The latex parameter can be used to optionally change how neutron count
        is formatted.

        """
        n = self.output_metadata['provided']['count']
        s = self.setup['src']['metadata']['energy_description']
        g = self.setup['geom']['decoded']['short_description']
        if latex:
            from ._common import _latex_format
            n = '$%s$'%_latex_format(n)
        return f'{n} {s} neutrons through {g}'

    def long_title( self ):
        """Returns a longer title string. This is less likely to fit nicely on
        top of a plot, but includes contents of the various configuration
        strings used to setup the simulation run.
        """
        s = self.setup['src']['cfgstr']
        g = self.setup['geom']['cfgstr']
        e = self.setup['engine']['cfgstr']
        r = f'"{s}" on "{g}"'
        n = self.output_metadata['tallied']['count']
        r = f'{r} ({n} fills)'
        return f'{r} ("{e}")' if e else r

    def plot(self, **kwargs):
        """Loops over all tallies and plots them, passing along any kwargs. The
        same functionality can be obtained by (if the current results object is
        "results"):

        for t in results.tallies:
            t.plot(**kwargs)

        """
        p = False
        for t in self.tallies:
            p = True
            t.plot(**kwargs)
        if not p:
            from ._common import ncwarn
            ncwarn('No tallies were enabled.')

    def plot_xsect(self, **kwargs):
        """Plots the material cross sections with the plot_xsect function from
        the NCrystal.plot module (passing along any kwargs). This is the same
        as (if the current results object is "results"):

        from NCrystal.plot import plot_xsect
        plot_xsect(results.setup['material']['cfgstr'],**kwargs)
        """
        from .plot import plot_xsect
        return plot_xsect(self.setup['material']['cfgstr'])

    def dump( self, do_print = True, prefix = '',
              tally_filter_fct = None ):
        """Prints a short summary of the simulation results, and also returns it
        as a string. If do_print is False, nothing will be printed. Optionally,
        all printed lines can be prefixed by setting the prefix parameter.

        The summary will also include a dump of the histogram of all the
        available tallies. To reduce the verbosity, one can optionally select
        which tallies to include in the dump by the tally_filter_fct
        parameter. For instance to only include the "mu" tally, one would use:

          tally_filter_fct = lambda tname : tname=="mu"

        """

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

    def __eq__(self, o ):
        if not isinstance(o,MMCResults):
            return False
        ds, do = self.__data, o._raw_data()
        return( id(self) == id(o) or id(ds)==id(do) or ds == do )

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
        """The tally name (e.g. "q", "mu", "theta", ...)."""
        return self.__data['tallyname']

    @property
    def unit( self ):
        """Return unit of tallied quantity like "eV", "Aa", or "1/Aa". Will be
        an empty string for unit-less quantities.
        """
        from ._mmc_impl import tally_info
        return tally_info()['hists'][self.name]['unit']

    @property
    def short_description( self ):
        """Return a string with a short description of the tallied quantity,
        like "wavelength" or "cosine scattering angle".
        """
        from ._mmc_impl import tally_info
        return tally_info()['hists'][self.name]['short_descr']

    def _raw_data( self ):
        return self.__data

    @property
    def mother( self ):
        """Return the MMCResults object which this MMCTallyView object belongs
        to."""
        return self.__mmcresults

    @property
    def hist_total( self ):
        """Return main histogram with all tallied events."""
        h = self.__data['total']
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
        if not isinstance(o,MMCTallyView):
            return False
        ds, do = self.__data, o._raw_data()
        return( id(self) == id(o) or id(ds)==id(do) or ds == do )

    @property
    def _nhists( self ):
        return 1 + (len(self.hist_breakdown) if self.hist_breakdown else 0)

    def histogram_sum( self, *, select=None, exclude=None ):
        """Creates a new histogram by adding up selected breakdown histograms.

        Both "select" and "exclude" parameters can either be one of the
        component names as a string, or a list of such strings.

        Examples:

        h = tally.histogram_sum(select=['NOSCAT','MULTISCAT_PUREELAS'])
        h = tally.histogram_sum(exclude='SINGLESCAT_ELAS')

        """

        if isinstance(exclude,str):
            exclude=[exclude]
        if isinstance(select,str):
            select=[select]
        if not exclude and not select:
            return self.hist_total.clone()
        histmap = self.hist_breakdown
        if select:
            histmap = dict( (hn,h) for hn,h in histmap.items()
                            if hn in select )
        if exclude:
            histmap = dict( (hn,h) for hn,h in histmap.items()
                            if hn not in exclude )
        hl = [h for hn,h in sorted(histmap.items()) ]
        if len(hl) <= 1:
            return hl[0].clone() if hl else None
        h = hl[0].clone()
        for o in hl[1:]:
            h.add_contents( o )
        return h

    def dump( self, *args, **kwargs):
        """Shorthand for .hist_total.dump(*args,**kwargs)."""
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
        """Launch a plot (via matplotlib) of a histogram of the tallied
        quantity.

        By default this will include a breakdown into the components
        (single/multi scattering, elastic/inelastic) contributing, if
        available. Setting the "breakdown" parameter to True or False overrides
        this.

        If title is None, "auto" or "short", a short title will be auto
        generated. If it is "long", a longer title with full configuration
        strings will be used. If title is "none" or False, no title will be
        shown. Finally, any other non-empty string provided will simply become
        the title.

        The max_nbins or rebin_factor parameters can be used to reduce the
        binning granularity. For instance, setting max_nbins=100 will ensure
        that at most 100 bins are shown.

        The parameters do_grid, do_legend, and logy are hopefully
        self-explanatory.

        Finally, for advanced users, follows a series of parameters which can be
        used by advanced users who wish to customise the plotting further: If
        plt is None, matplotlib.pyplot will be used, and if do_newfig is True, a
        call to plt.figure() then follows. Setting do_newfig='auto' will cause a
        new figure to be created, unless axis is provided. Finally, if axis is
        None, it will default to plt.gca(). After all plotting calls have been
        done, plt.show() will be called unless do_show is False.

        The return value of the function is the plt object actually used
        (normally matplotlib.pyplot).

        """

        if title in (None,'auto','short'):
            title = self.__mmcresults.short_title(latex=True)
        elif title == 'long':
            title = self.__mmcresults.long_title(latex=True)
        title = (title or '').strip() or False

        from ._mmc_impl import _plot_tally
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

