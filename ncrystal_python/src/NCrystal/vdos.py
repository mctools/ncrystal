
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

"""Module with VDOS-related utilities"""

from .constants import constant_planck as _constant_planck
from .constants import constant_boltzmann as _constant_boltzmann
from .constants import constant_c as _constant_c

_is_unit_test = False

vdos_units_2_eV = {
    'eV' : 1.0,
    'meV' : 1e-3,
    'keV' : 1e3,
    'MeV' : 1e6,
    'THz' : _constant_planck*1e12,
    'GHz' : _constant_planck*1e9,
    'MHz' : _constant_planck*1e6,
    'kHz' : _constant_planck*1e3,
    'Hz' : _constant_planck,
    '1/cm' : _constant_planck * _constant_c * 1e-8,# 1e-8 to get c in [cm/s] instead of [Aa/s]
}

def createVDOSDebye( debye_temperature ):
    """Create simplified VDOS according to the Debye model"""
    from ._numpy import _np, _ensure_numpy, _np_linspace
    _ensure_numpy()
    #NB: Must keep function exactly synchronised with createVDOSDebye function
    #in .cc src (although leaving out temperature,boundXS,elementMassAMU args
    #here):
    debye_energy = _constant_boltzmann*debye_temperature
    if 0.5*debye_energy < 1.001e-5:
        from .exceptions import NCBadInput
        raise NCBadInput('Too low Debye temperature')
    vdos_egrid = _np_linspace(0.5*debye_energy,debye_energy,20)
    scale = 1.0 / (debye_energy*debye_energy)
    vdos_density = scale * (vdos_egrid**2)
    #Actual returned egrid should contain only first and last value:
    return (_np.asarray([vdos_egrid[0],vdos_egrid[-1]]) ,vdos_density)

def debyeIsotropicMSD( *, debye_temperature, temperature, mass ):
    """Estimate (isotropic, harmonic) atomic mean-squared-displacement (a.k.a. "U_iso") using the
       Debye Model (eq. 11+12 in R.J. Glauber, Phys. Rev. Vol98 num 6,
       1955). Unit of returned MSD value is Aa^2. Input temperatures should be
       in Kelvin, and input atomic mass should be in amu.
    """
    from ._chooks import _get_raw_cfcts
    return float(_get_raw_cfcts()['ncrystal_debyetemp2msd'](debye_temperature, temperature, mass))

def debyeTempFromIsotropicMSD( *, msd, temperature, mass ):
    """The inverse of debyeIsotropicMSD (implemented via root-finding), allowing to
       get the Debye temperature which will give rise to a given
       mean-squared-displacement (a.k.a. "U_iso").
    """
    from ._chooks import _get_raw_cfcts
    return float(_get_raw_cfcts()['ncrystal_msd2debyetemp'](msd, temperature, mass))

def analyseVDOS(emin,emax,density,temperature,atom_mass_amu):
    """Analyse VDOS curve to extract mean-squared-displacements, Debye temperature,
    effective temperature, gamma0 and integral. Input VDOS must be defined via
    an array of density values, over an equidistant energy grid over [emin,emax]
    (in eV). Additionally, it is required that emin>0, and a parabolic trend
    towards (0,0) will be assumed for energies in [0,emin]. Units are kelvin and
    eV where appropriate.
    """
    from ._chooks import _get_raw_cfcts
    from ._numpy import _np, _ensure_numpy
    _ensure_numpy()
    density = _np.asarray(density,dtype=float)
    return _get_raw_cfcts()['nc_vdoseval'](emin,emax,density,temperature,atom_mass_amu)

def extractGn( vdos, n, mass_amu, temperature, scatxs = 1.0, expand_egrid = True ):
    """Extract Sjolander's Gn function of order n."""
    assert 1 <= n <= 99999
    from .misc import AnyVDOS
    v = AnyVDOS(vdos)
    from ._chooks import _get_raw_cfcts
    emin, emax, Gn =  _get_raw_cfcts()['raw_vdos2gn'](v.egrid(),v.dos(),scatxs, mass_amu, temperature, int(n) )
    if not expand_egrid:
        return (emin,emax),Gn
    else:
        from ._numpy import _ensure_numpy, _np_linspace
        _ensure_numpy()
        return _np_linspace( emin, emax, len(Gn) ), Gn

def extractKnl( vdos, mass_amu, temperature, vdoslux = 3, scatxs = 1.0,
                order_weight_fct = None, target_emax = None,
                plot = False, **plotkwargs ):
    """Expand the VDOS to a scattering kernel, based on the provided
    parameters.

    If provided, the order_weight_function can be used to assign a weight to
    each order of the expansion (e.g. 0.0 to remove the contribution of that
    order). It must be a function taking a single parameter n (the order), and
    return a floating point value (the weight).

    The target_emax parameter can be used to override the Emax value (eV)
    targeted by the expansion, which will otherwise be based on the vdoslux
    value.

    The results will include alpha, beta, and sab arrays, as well as a
    suggested_emax value (eV) which provides the actual Emax value of the
    kernel. The latter will be None if an order_weight_function is provided.

    If plot=True, the extracted kernel will be plotted with the
    NCrystal.plot.plot_knl function.

    """
    from .misc import AnyVDOS
    v = AnyVDOS(vdos)
    from ._chooks import _get_raw_cfcts
    f = _get_raw_cfcts()['raw_vdos2knl']
    a, b, sab, suggested_emax = f( v.egrid(),
                                   v.dos(),
                                   scatxs,
                                   mass_amu,
                                   temperature,
                                   vdoslux,
                                   order_weight_fct,
                                   target_emax )
    k = dict( alpha = a,
              beta = b,
              sab = sab,
              mass_amu = mass_amu,
              temperature = temperature,
              scatxs = scatxs,
              suggested_emax = suggested_emax )
    if plot:
        from .plot import plot_knl
        plot_knl( k, **plotkwargs )
    return k


class PhononDOSAnalyser:

    """Immutable class which reads and facilitates interpretation of phonon DOS information.

    The input data format (fmt) can either be "raw", which is a list of DOS
    (label,egrid,density). It can also be an NCrystal.misc.AnyVDOS,
    NCrystal.Info.DI_VDOS or NCrystal.Info.DI_VDOSDebye object, or a list of
    such. Or it can be "quantumespresso", pointing to a matdyn.dos file produced
    by QuantumEspresso, expecting a format like (leftmost column is the
    frequency, the next column is a combined DOS and the rest are partial DOS's
    for the atoms in the same order as input into QE):

     # Frequency[cm^-1] DOS PDOS
     -2.2774053913E+02  0.0000000000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
     -2.2674053913E+02  3.3672834745E-05  1.0609E-05  1.0446E-05  6.2610E-06  6.3576E-06
     -2.2574053913E+02  7.5113401337E-05  2.3699E-05  2.3292E-05  1.3952E-05  1.4171E-05
     -2.2474053913E+02  1.1552338186E-04  3.6486E-05  3.5824E-05  2.1434E-05  2.1780E-05
     ...

    The partial DOS columns (the 4 rightmost columns in the example above) will
    be labelled "pdos_0", "pdos_1", etc. After construction, the .plot(..)
    method is then typically used to inspect the data, and the labels can be
    merged, dropped or changed using the methods .drop(..), .merge(..) or
    .update_label(..).

    Next, one might way to apply a lower threshold on the curves if needed with
    .apply_cutoff(..), after choosing an appropriate threshold value by look at
    the plots produced by plot_cutoff_effects and plot_cutoff_effects_on_xsects.

    A word warning: The PhononDOSAnalyser is immutable, hence this code does not
    do anything (mydosana is an PhononDOSAnalyser object):

      mydosana.apply_cutoff( 0.015 )

    Instead one must do:

      mydosana_with_cutoff = mydosana.apply_cutoff()

    Or simply (if you don't care to keep the old un-edited object around):

      mydosana = mydosana.apply_cutoff()

    In addition to setting a threshold with .apply_cutoff, for output files
    intended to have a high quality but smaller footprint (e.g. files destined
    for the NCrystal stdlib) one might also consider applying a regularisation,
    i.e. ensuring that the energy grid is linearly spaced and all grid points
    being a multiple of the binwidth (so that 0.0 would be a grid point, if the
    grid would have been extended downwards).

    Once the phonon curves are as desired, one most likely want to use the
    .apply_to(..) method to apply all or selected VDOS curves to an
    NCMATComposer instance (which one might then for instance use to write the
    data to an NCMAT file).

    Note that many methods on this class can optionally take a list of labels or
    indices. If so, the operation will only affect those DOS curves.

    """

    def __init__(self, data, fmt = None ):
        """Initialise from data (see class description)."""
        from .exceptions import NCBadInput

        if isinstance(data,tuple) and len(data)==2 and data[0]=='__internal_state__' and isinstance(data[1],dict):
            import copy
            self.__d = copy.deepcopy( data[1] )
            return

        if fmt is None:
            from .misc import AnyVDOS
            from .core import Info
            def _is_anyvdos(x):
                return any( isinstance(data,e)
                            for e in (AnyVDOS,Info.DI_VDOS,Info.DI_VDOSDebye) )
            if _is_anyvdos( data ):
                fmt = 'anyvdos'
            elif data and hasattr(data,'__len__') and all( _is_anyvdos(e) for e in data ):
                fmt = 'anyvdos_list'
            elif data and hasattr(data,'__len__') and all( ( hasattr(e,'__len__') and len(e)==3 ) for e in data ):
                fmt = 'raw'
            else:
                fmt = 'quantumespresso'

        def _extract_from_anyvdos(x):
            from .misc import AnyVDOS
            v = AnyVDOS( x )
            return v.label,v.egrid(),v.dos()

        _srcname = None
        if fmt == 'anyvdos':
            doslist = [_extract_from_anyvdos(data)]
        elif fmt == 'anyvdos_list':
            doslist = [_extract_from_anyvdos(e) for e in data]
        elif fmt == 'quantumespresso':
            from .misc import AnyTextData
            td = AnyTextData( data )
            doslist = _read_quantumespresso( td.content )
            _srcname = td.name
        elif fmt == 'raw':
            import copy
            doslist = copy.deepcopy( data )
            _srcname = 'raw data'
        else:
            raise NCBadInput('Invalid input format. Currently only "raw", "anyvdos", "anyvdos_list", and "quantumespresso" are supported.')

        _dl = []
        from ._numpy import _ensure_numpy, _np
        _ensure_numpy()

        nminpts = 5
        for lbl,egrid,density in doslist:
            egrid = _np.asarray(egrid,dtype=float).copy()
            if len(egrid) < nminpts:
                raise NCBadInput(f'DOS egrid has too few points (at least {nminpts} required).')
            if len(density) != len(egrid):
                raise NCBadInput('DOS egrid and density arrays have different lengths.')
            def _is_grid( a ):
                return _np.all(a[:-1] < a[1:])#stackoverflow question 47004506
                                              #but with < instead of <=
            if not _is_grid(egrid):
                raise NCBadInput('DOS egrid does not consist of increasing unique values')
            if _np.isinf(egrid[-1]) or not egrid[-1] > 0.0:
                raise NCBadInput(f'DOS egrid has invalid upper edge value: {egrid[-1]}')

            density = _np.asarray(density,dtype=float).copy()
            _densmin, _densmax = density.min(), density.max()
            if not _densmax > 0.0:
                raise NCBadInput('Maximum density value is not > 0.0.')
            if _densmin < 0.0:
                if abs(_densmin) < 1e-3*_densmax:
                    density.clip( min = 0.0, out = density )
                    from . import _common as nc_common
                    nc_common.warn('Clipping tiny negative DOS density values to 0.0')
                else:
                    raise NCBadInput('Negative density values observed.')

            #peel off excess zeroes:
            assert nminpts>3
            while len(egrid) > nminpts and egrid[-3]>0.0 and density[-2]==0.0 and density[-1]==0.0:
                egrid = egrid[:-1]
                density = density[:-1]
            while len(egrid) > 5 and density[0]==0.0 and density[1]==0.0:
                egrid = egrid[1:]
                density = density[1:]
            if not density.max() > 0.0:
                raise NCBadInput('Invalid input DOS (non-positive everywhere)')
            _dl.append( ( str(lbl), egrid, density ) )
        doslist = _dl

        self.__d = dict( doslist = doslist, srcname = _srcname )

    @property
    def source_name( self ):
        "Name of source data (such as the file name). Might be absent (None)."
        return self.__d['srcname']

    @property
    def ndos( self ):
        """Number of DOS curves."""
        return len( self.__d['doslist'])

    @property
    def labels( self ):
        """Available DOS labels."""
        return list(lbl for lbl,_,_ in self.__d['doslist'])

    def dos( self, label_or_idx ):
        """Get DOS by index or label. Returns tuple of two arrays: (egrid,dos)."""
        lbl,egrid,dos = self.__d['doslist'][self.__dosidx( label_or_idx )]
        return egrid,dos

    def drop( self, *dos_labels_or_indices ):
        """
        Return a new instance from which the indicated DOS has been removed.
        DOS can be identified via labels or indices, and it is possible to
        provide a list or tuple of them, to drop multiple
        """
        droplist = self.__dosidxlist( *dos_labels_or_indices )
        if not droplist:
            return self
        o = self.__clone()
        o.__d['doslist'] = [ e for i,e in enumerate(o.__d['doslist'])
                             if i not in droplist ]
        return o

    def update_label( self, old_label_or_idx, newlabel ):
        """Return a new instance in which the label of the indicated DOS curve
        has been updated."""
        idx = self.__dosidx( old_label_or_idx )
        o = self.__clone()
        dl = list(o.__d['doslist'])
        ll = list(dl[idx])
        ll[0] = str(newlabel)
        dl[idx] = tuple(ll)
        o.__d['doslist'] = tuple(dl)
        return o

    def merge( self, *dos_labels_or_indices, weights = None, newlabel = None ):
        """
        Return a new instance from which the indicated DOS curves have been merged.
        If weights is given, it must be a list of weights to to use in this merging.

        If newlabel is not provided, the label of the resulting curve will be
        automatically generated.
        """
        mergelist = self.__dosidxlist( *dos_labels_or_indices )
        if weights and len(weights) != len(mergelist):
            from .exceptions import NCBadInput
            raise NCBadInput('Invalid weights provided (must be list of same length as the number of DOS curves to merge')
        if len(mergelist)==1 and newlabel:
            #simply a label update:
            return self.update_label( mergelist[0], newlabel )
        if len(mergelist)<=1:
            #do nothing:
            return self

        egrids = [(lbl,egrid) for i,(lbl,egrid,dos) in enumerate(self.__d['doslist']) if i in mergelist]
        from ._numpy import _ensure_numpy, _np
        _ensure_numpy()
        newegrid = egrids[0][1].copy()
        for lbl,eg in egrids[1:]:
            if not _np.array_equal(newegrid,eg):
                from .exceptions import NCBadInput
                raise NCBadInput('Can not merge DOS curves with incompatible egrids (problems merging labels "%s" and "%s")'%(egrids[0][0],lbl))

        if weights is None:
            weights = [ 1.0 ] * len(mergelist)
        assert all( w >= 0.0 for w in weights )
        import math
        _wsum = math.fsum(w for w in weights)
        assert _wsum > 0.0
        weights = [ w/_wsum for w in weights ]

        c = None
        for idx,w in zip(mergelist,weights):
            _ = self.dos(idx)[1].copy()
            _ *= w
            if c is None:
                c = _
            else:
                c += _
        assert c is not None

        if newlabel is None:
            newlabel = 'merged(%s)'%(','.join( [lbl for i,(lbl,_,_) in
                                                enumerate(self.__d['doslist'])
                                                if i in mergelist] ))
            while newlabel in self.labels:
                newlabel += '(uniquelabel)'
        o = self.drop( mergelist )
        o.__d['doslist'].append( (newlabel,newegrid,c) )
        return o

    def apply_regularisation( self, target_n, *dos_labels_or_indices, quiet = True ):
        """Returns a new instance where the indicated lower threshold value have been
        "regularised", i.e. put on a linearly spaced grid which if it had been
        extended downwards, would eventually cross the point E=0. It is an error
        to apply it to a curve whose first egrid point is not positive (use the
        .apply_cutoff(..) method first for such curves)."""
        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )
        if not selected:
            return self
        target_n = int(target_n)
        assert target_n >= 10
        def doregfct( _e, _d ):
            return _do_regularise( egrid=_e, density=_d,
                                   n=target_n, quiet = quiet)
        o = self.__clone()
        oldlist = o.__d['doslist']
        newlist = []
        for idx,(lbl,egrid,dos) in enumerate(oldlist):
            if idx not in selected:
                newlist.append( (lbl,egrid,dos) )
                continue
            assert len(egrid)==len(dos)
            if not egrid[0] > 0.0:
                from .exceptions import NCBadInput
                raise NCBadInput('Can not regularise DOS whose first egrid value '
                                 'is not positive. Use the .apply_cutoff method first.')
            egrid, dos = doregfct( egrid, dos )
            assert len(egrid)==len(dos)
            newlist.append( (lbl,egrid,dos) )
        o.__d['doslist'] = newlist
        return o

    def apply_cutoff( self, threshold, *dos_labels_or_indices ):
        """Returns a new instance where the indicated lower threshold value has been
        applied to all (selected) DOS curves."""
        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )

        if not selected:
            return self

        threshold = self._parse_threshold( threshold )

        from ._numpy import _ensure_numpy, _np
        _ensure_numpy()

        o = self.__clone()
        oldlist = o.__d['doslist']
        newlist = []
        for idx,(lbl,egrid,dos) in enumerate(oldlist):
            if idx not in selected:
                newlist.append( (lbl,egrid,dos) )
                continue
            assert len(egrid)==len(dos)
            j = _np.argmax( egrid >= threshold )
            egrid,dos = egrid[j:],dos[j:]
            assert egrid[0] >= threshold
            assert len(egrid)==len(dos)
            newlist.append( (lbl,egrid,dos) )
        o.__d['doslist'] = newlist
        return o

    def plot( self, *dos_labels_or_indices,
              do_newfig = True,
              do_show = True,
              do_legend = True,
              logy=None,
              unit = 'eV',
              ymin = None,
              ymax = None,
              xmin = None,
              xmax = None,
             ):
        """Plot contained DOS curves for the selected DOS labels (or indices),
        defaulting to all curves. Set do_show to false to avoid plt.show(), logy
        to default to a semilogy plot, ymin to set a minimum plotrange of the
        y-axis, and finally the unit argument can be used to show the DOS in
        other units (e.g. "THz", "1/cm", "meV", etc.
        """
        self.__plot( *dos_labels_or_indices, do_newfig = do_newfig,
                     do_show = do_show, logy=logy, unit = unit,
                     do_legend=do_legend,
                     ymin = ymin, ymax = ymax,
                     xmin = xmin, xmax = xmax )

    def plot_gn( self, *dos_labels_or_indices, n=1,
                 temperature = 293.15, masses = None,
                 do_newfig = True, do_show = True, do_legend = True,
                 logy=None, unit = 'eV' ):
        """Similar to .plot() but showing the Sjolander Gn function instead. If
           labels are not elements or isotopes, masses must be provided in a
           list (in daltons). It is an error to try to plot a Gn function for a
           vdos curve whose first egrid point is not positive (i.e. you must
           first .apply_cutoff(..)).
        """
        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )
        self.__plot( *selected, do_newfig = do_newfig,
                     do_show = do_show, logy=logy, unit = unit,
                     do_legend=do_legend,
                     sjolanderGn = self.__sjolanderGn_args( selected = selected,
                                                            n=n,
                                                            masses=masses,
                                                            temperature=temperature) )

    def plot_cutoff_effects( self, thresholds, *dos_labels_or_indices, gn = None,
                             temperature = 293.15, masses = None,
                             regularise_n_values = None, **plot_kwargs ):
        """Plot the effect of the listed thresholds on the DOS curves for the
        selected DOS labels (or indices), defaulting to all curves. Any
        plot_kwargs will be used as on the .plot() method.

        To see the effect on Sjolander's Gn curves instead, supply the gn
        keyword (i.e. gn=1 to see G1) -- optionally along with temperature and
        masses keywords (cf. .plot_gn(..)).

        """
        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )

        if gn is not None:
            sjolanderGn = self.__sjolanderGn_args( selected = selected,
                                                   n=gn,
                                                   masses=masses,
                                                   temperature=temperature)
            plot_kwargs['sjolanderGn'] = sjolanderGn

        unitname,unitfactor = _parsevdosunit( plot_kwargs.get('unit','eV') )
        do_show = plot_kwargs.get('do_show',True)
        do_newfig = plot_kwargs.get('do_newfig',True)
        color_offset = plot_kwargs.get('color_offset',0)
        plot_kwargs['do_show'] = False
        plot_kwargs['do_newfig'] = False
        do_legend = plot_kwargs.get('do_legend',True)
        do_grid = plot_kwargs.get('do_grid',True)
        plot_kwargs['do_grid'] = False
        plot_kwargs['do_legend'] = False

        from .plot import _import_matplotlib_plt
        plt = _import_matplotlib_plt()
        if do_newfig:
            plt.figure()

        plot_thresholds = [e for e in thresholds]
        if gn is None:
            #Can not plot Gn without a cutoff
            plot_thresholds = [None] + plot_thresholds
        reg_n_vals = [None] + ( [ e for e in regularise_n_values] if regularise_n_values else [] )
        for t in plot_thresholds:
            for regn in reg_n_vals:
                plot_kwargs['color_offset'] = color_offset
                if t is None:
                    if regn is None:
                        lblcomments=['orig']
                    else:
                        lblcomments = [f'npts={regn}']
                else:
                    t_parsed = self._parse_threshold( t )
                    lblcomments = [f'cut@{t_parsed/unitfactor:g}{unitname}']
                    if regn is not None:
                        lblcomments+=[f'npts={regn}']
                lblcomment = ','.join(lblcomments)
                plot_kwargs['labelfct'] = lambda lbl : f'{lbl} ({lblcomment})'
                color_offset += len(selected)
                o = self if t is None else self.apply_cutoff( t, *selected )
                if regn is not None:
                    if t is None:
                        #Only go ahead in this case, if all selected egrids already have a cutoff:
                        if not all( self.dos(idx)[0][0]>0.0 for idx in selected ):
                            continue
                    o = o.apply_regularisation( regn, *selected )
                o.__plot(*selected,**plot_kwargs)

        from .plot import _plt_final
        _plt_final(do_grid,do_legend,do_show)

    def plot_cutoff_effects_on_xsects( self, ncmatcomposer, thresholds, cfg_params = None,
                                       *dos_labels_or_indices, lblmap = None, **plot_kwargs ):
        from .ncmat import NCMATComposer
        assert isinstance(ncmatcomposer,NCMATComposer), ( "First argument in call to .plot_cutoff_"
                                                          "effects_on_xsects(..) must be an NCMATComposer object" )
        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )
        unitname,unitfactor = _parsevdosunit( plot_kwargs.get('unit','eV') )
        do_show = plot_kwargs.get('do_show',True)
        do_newfig = plot_kwargs.get('do_newfig',True)
        do_grid = plot_kwargs.get('do_grid',True)
        do_legend = plot_kwargs.get('do_legend',True)
        plot_kwargs['scatter_breakdown'] = False
        plot_kwargs['show_scattering'] = True
        plot_kwargs['show_absorption'] = False
        plot_kwargs['do_newfig'] = False
        plot_kwargs['do_show'] = False
        plot_kwargs['do_grid'] = False
        plot_kwargs['do_legend'] = False
        plot_kwargs['cfg_params'] = cfg_params
        color = plot_kwargs.get('color')
        lblmap = self.__determine_lblmap( selected, ncmatcomposer, lblmap = lblmap, warn = True )
        lbls = list(sorted(lblmap.keys()))
        if not lbls:
            return
        from .plot import _import_matplotlib_plt
        plt = _import_matplotlib_plt()
        if do_newfig:
            plt.figure()
        colorder = self.__colorder()
        for iplot, t in enumerate([e for e in thresholds]):
            t = self._parse_threshold( t )
            c = ncmatcomposer.clone()
            o = self.apply_cutoff( t, lbls )
            thr_description = f'cut@{t/unitfactor:g}{unitname}'
            for lbl in lbls:
                if lbl not in lblmap:
                    from . import _common as nc_common
                    nc_common.warn('Not using PhononDOSAnalyser label "{lbl}" in plot.')
                    continue
                c.set_dyninfo_vdos( lblmap[lbl], comment = 'From PhononDOSAnalyser', **o.get_dyninfo_args(lbl) )
            if not color:
                plot_kwargs['color'] = colorder[iplot%len(colorder)]
            plot_kwargs['labelfct'] = lambda x : thr_description
            c.plot_xsect( **plot_kwargs )

        if do_legend:
            plt.legend()
        if do_grid:
            plt.grid()
        t = 'DOS cutoff effect'
        if cfg_params:
            t += ' (%s)'%cfg_params.strip()
        plt.title(t)
        if do_show:
            plt.show()

    def apply_to( self, ncmatcomposer, *dos_labels_or_indices, lblmap = None, warn = True, cutoff = None ):
        """Apply DOS curves to NCMATComposer objects, resulting in updates to
        the relevant dyninfo sections. If lblmap is not given, the
        determine_mapping_to_composer_labels() method is used to infer one
        automatically.

        If the cutoff value is provided, it will be applied to all DOS curves
        first. Otherwise, a warning will be emitted and an ad-hoc threshold will
        be supplied to any curves whose initial frequency point is not positive.

        """
        from .ncmat import NCMATComposer
        assert isinstance(ncmatcomposer,NCMATComposer), ( "First argument in call to .apply"
                                                          "(..) must be an NCMATComposer object" )

        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )
        if not selected:
            return

        from . import _common as nc_common
        warnfct = nc_common.warn if warn else (lambda s : None)

        lblmap = self.__determine_lblmap( selected, ncmatcomposer, lblmap = lblmap, warn = warn )
        cutoff = ( self._parse_threshold( cutoff ) if cutoff is not None else None ) or 0.0
        def _access_egrid( idx ):
            return self.__d['doslist'][idx][1]
        selected_needscutoff = [ idx for idx in selected if not _access_egrid(idx)[0] > cutoff ]
        if selected_needscutoff and not cutoff:
            #must autodetermine a suitable cutoff. We take it as 1% of the maximum egrid value:
            cutoff = 0.03 * max( _access_egrid(idx)[-1] for idx in selected_needscutoff )
            warnfct(f'Applying an ad-hoc lower DOS egrid cutoff of {cutoff:g}eV since the cutoff parameter'
                    ' was not provided. For important work, it is recommended to instead select one'
                    ' explicitly (probably after investigating - for instance with the .plot_*() methods).')

        o = self.apply_cutoff( cutoff, *selected_needscutoff ) if (cutoff and selected_needscutoff) else self
        if not lblmap:
            warnfct('Not applying any DOS curves to NCMATComposer object'
                    ' (perhaps try again with the lblmap argument).')
            return
        for srclbl, tgtlbl in sorted(lblmap.items()):
            diargs = o.get_dyninfo_args(srclbl)
            assert diargs['vdos_egrid'][0] > 0.0
            ncmatcomposer.set_dyninfo_vdos( tgtlbl,
                                            comment = f'From PhononDOSAnalyser ("{self.source_name}", atom with label "{srclbl}")',
                                            **diargs )

    def is_regular( self, *dos_labels_or_indices ):
        """A regular DOS has an equidistant energy grid whose first point, e0,
           is positive and a binwidth which divides e0 an even number of times."""
        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )
        if not selected:
            return True
        return all( _vdos_egrid_is_regular(self.__d['doslist'][idx][1]) for idx in selected )

    def requires_cutoff( self, *dos_labels_or_indices ):
        """Whether or not any of the selected curves has an initial energy grid
           point which is not positive."""
        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )
        return any( (not (self.__d['doslist'][idx][1][0]>0.0) ) for idx in selected )

    def get_dyninfo_args( self, label_or_idx ):
        """returns a dict with 'vdos_egrid' and 'vdos' keys, suitable for usage
           when calling NCMATComposer.set_dyninfo_vdos(..)"""
        lbl,egrid,dos = self.__d['doslist'][self.__dosidx( label_or_idx )]
        return dict( vdos_egrid = egrid, vdos = dos )

    def determine_mapping_to_composer_labels( self, ncmatcomposer, warn = True ):
        """Try to determine a mapping between labels in this object and an
        NCMATComposer object. The mapping will either be based on Z-values or
        the actual label names themselves.

        Unless warn=False, warnings will be emitted when a label could not be mapped.

        """
        from .ncmat import NCMATComposer
        assert isinstance(ncmatcomposer,NCMATComposer), ( "First argument in call to .determine"
                                                          "_mapping_to_composer_labels(..) must be an NCMATComposer object" )
        #First check if any labels can be mapped based on Z-values:
        from . import _common as nc_common
        from . import atomdata as nc_atomdata
        warnfct = nc_common.warn if warn else (lambda s : None)
        lbl2z = {}
        self_labels = self.labels
        for lbl in self_labels:
            elemiso = nc_common.check_elem_or_isotope_marker( lbl )
            if elemiso:
                lbl2z[lbl] = nc_atomdata.elementNameToZValue(elemiso,allow_isotopes=True)
        lbl2tgt_zbased = {}
        for z in set(lbl2z.values()):
            lbls = set( _lbl for _lbl,_z in lbl2z.items() if _z==z )
            if len(lbls) == 1:
                lbl = lbls.pop()
                tgtlbl = ncmatcomposer.find_label( z, allow_multi = False )
                if tgtlbl:
                    lbl2tgt_zbased[lbl] = tgtlbl

        #Now, see if we can map based on having the same label in the two objects:
        common_lbls = set(self_labels).intersection(set( ncmatcomposer.get_labels() ))

        res = {}
        for lbl in self_labels:
            tgt_z = lbl2tgt_zbased.get(lbl)
            tgt_lblname = lbl if lbl in common_lbls else None
            if tgt_z and tgt_lblname and tgt_z != tgt_lblname:
                warnfct( f'Could not map atom with PhononDOSAnalyser label "{lbl}" '
                         'unambiguously to NCMATComposer (Z value implies mapping to'
                         f' "{tgt_z}" but the label name implies the target label "{tgt_lblname}")')
            else:
                tgtlbl = tgt_z or tgt_lblname
                if tgtlbl:
                    res[lbl] = tgtlbl
                else:
                    warnfct( f'Could not map atom with PhononDOSAnalyser label "{lbl}" '
                             'unambiguously to NCMATComposer label' )
        return res

    def __dosidx( self, label_or_idx ):
        import numbers
        ll = self.__d['doslist']
        if isinstance(label_or_idx,numbers.Integral):
            i = int( label_or_idx )
            if not 0 <= i < len(ll):
                from .exceptions import NCBadInput
                raise NCBadInput('Index out of range (%i is not in range 0..%i'%(i,len(ll)-1))
            return i
        _ = [i for i,(lbl,_,_) in enumerate(ll) if lbl == label_or_idx ]
        if not _:
            from .exceptions import NCBadInput
            _ = '","'.join(lbl for lbl,_,_ in ll)
            raise NCBadInput('Invalid label "%s" (available labels are "%s")'%(label_or_idx,_))
        assert len(_) == 1
        return _[0]

    def __dosidxlist( self,  *dos_labels_or_indices, default_is_all = False ):
        if not dos_labels_or_indices:
            return list( range( self.ndos ) ) if default_is_all else []
        res = []
        for ll in dos_labels_or_indices:
            if not hasattr(ll,'__len__') or hasattr(ll,'startswith'):
                #single item
                res.append( self.__dosidx( ll ) )
            else:
                res += [ self.__dosidx( e ) for e in ll ]
        return res

    def __sjolanderGn_args( self, selected, n=1, masses = None, temperature = 293.15 ):
        assert 1<=n<=9999
        assert temperature >= 0.001
        from .exceptions import NCBadInput
        if masses is None:
            masses = []
            from .atomdata import atomDB
            for idx in selected:
                lbl = self.__d['doslist'][idx][0]
                ad = atomDB(lbl,throwOnErrors=False)
                if not ad:
                    raise NCBadInput( 'Can not plot Gn function for label "lbl" which does'
                                      ' not correspond to a known element or isotope. Either'
                                      ' change the label with .update_label(..), or directly'
                                      ' provide masses using the "masses" parameter' )
                masses.append( ad.averageMassAMU() )
        if len(masses) != len(selected):
            raise NCBadInput( 'Invalid number of masses provided' )
        return dict( n = n, masses = masses, temperature=temperature )

    def __colorder( self ):
        from ._common import _palette_Few as _palette
        return [ _palette[e] for e in ['red',
                                       'blue',
                                       'orange',
                                       'green',
                                       'purple',
                                       'yellow',
                                       'pink',
                                       'gray',
                                       ]]

    def __plot( self, *dos_labels_or_indices, do_newfig = True, do_show = True, logy=None, unit = 'eV',
                color_offset = 0, labelfct = None, do_legend=True,do_grid=True, sjolanderGn = None,
                ymin=None, ymax = None, xmin=None, xmax=None ):

        selected = self.__dosidxlist( *dos_labels_or_indices, default_is_all = True )
        gnfcts = {}
        if sjolanderGn is not None:
            for idx, mass in zip(selected,sjolanderGn['masses']):
                lbl, egrid, dos = self.__d['doslist'][idx]
                if not egrid[0]>0.0:
                    from .exceptions import NCBadInput
                    raise NCBadInput('Can not plot Gn functions for DOS curves whose initial grid'
                                     ' point is not at a positive value (use apply_cutoff(..) to rectify')
                gnfcts[idx] = extractGn( vdos = (egrid,dos),
                                         n = sjolanderGn['n'],
                                         mass_amu = mass,
                                         temperature = sjolanderGn['temperature'] )



        unitname,unitfactor = _parsevdosunit( unit )
        from .plot import _import_matplotlib_plt
        plt = _import_matplotlib_plt()
        if do_newfig:
            plt.figure()
        from ._numpy import _ensure_numpy, _np_linspace

        colorder = self.__colorder()
        for idx,(lbl, egrid, dos) in enumerate( self.__d['doslist'] ):
            if idx not in selected:
                continue
            gnfct = gnfcts.get(idx)
            color = colorder[(color_offset+idx)%len(colorder)]
            actual_lbl = lbl if labelfct is None else labelfct(lbl)
            if gnfct:
                _x,_y = gnfct[0]/unitfactor, gnfct[1]
                if _is_unit_test:
                    def _fixup_y(y):
                        from ._numpy import _ensure_numpy, _np
                        _ensure_numpy()
                        assert y.min() >= 0.0
                        return _np.clip(y,y.max()*1e-13,None)
                    _y = _fixup_y(_y)#discard tiny values
                plt.plot( _x,_y,
                          label = actual_lbl,
                          color = color )
            else:
                plt.plot( egrid/unitfactor, dos,
                          label = actual_lbl,
                          color = color )
                if egrid[0]>0.0:
                    _k = dos[0] / egrid[0]**2
                    _ensure_numpy()
                    _x = _np_linspace(0.0, egrid[0], 2000+2)[1:-1]
                    plt.plot( _x/unitfactor, _k*(_x**2),ls=':',color=color)

        plt.xlabel('Frequency (%s)'%unitname)
        if sjolanderGn is not None:
            plt.ylabel('G%i (arbitrary scale)'%sjolanderGn['n'])
        else:
            plt.ylabel('DOS (arbitrary scale)')
        if ymin is not None or ymax is not None:
            plt.ylim(ymin,ymax)
        if xmin is not None or xmax is not None:
            plt.xlim(xmin,xmax)
        from .plot import _plt_final
        _plt_final(do_grid,do_legend,do_show,logy=logy)

    def __clone( self ):
        return PhononDOSAnalyser( ('__internal_state__',self.__d) )

    def _parse_threshold( self, value ):
        import numbers
        from ._common import _decodeflt
        from .exceptions import NCBadInput
        def impl(x):
            if isinstance(x,numbers.Real):
                return float(x)
            if isinstance(x,str):
                x=x.strip()
                v = None
                for un,unval in vdos_units_2_eV.items():
                    if x.endswith(un):
                        v = _decodeflt(x[:-len(un)].strip())
                        if v is not None:
                            return v *unval
                if v is None:
                    raise NCBadInput(f'Invalid threshold string: {x}')
            if hasattr(x,'__len__') and len(x)==2 and isinstance(x[1],str):
                _,unitfactor = _parsevdosunit( x[1] )
                v = _decodeflt( x[0] )
                if v is None:
                    raise NCBadInput(f'Invalid threshold value: {x[0]}')
                return v*unitfactor
            raise NCBadInput(f'Invalid threshold: {x}')

        v = impl(value)
        if not v>0.0 or not v < 1e6:
            raise NCBadInput(f'Invalid threshold value (out of range): {v:g}')
        return v

    def __determine_lblmap( self, selected, ncmatcomposer, lblmap = None, warn = True ):
        from .exceptions import NCBadInput
        if lblmap is None:
            lblmap = self.determine_mapping_to_composer_labels( ncmatcomposer, warn = warn )
        else:
            #use provided lblmap, but with a few sanity checks:
            composer_lbls = ncmatcomposer.get_labels()
            missing = set(lblmap.values()).difference(composer_lbls)
            def fmtlabellist(lbls):
                return ( '"%s"'%('", "'.join(lbls)) if lbls else '' )
            if missing:
                raise NCBadInput('Some values in lblmap are not present in provided'
                                 f' NCMATComposer object: {fmtlabellist(missing)} (the'
                                 f' following labels are available: {fmtlabellist(composer_lbls)})')
            missing = set(lblmap.keys()).difference(self.labels)
            if missing:
                raise NCBadInput('Some keys in lblmap are not actual labels in PhononDOSAnalyser'
                                 f' object: {fmtlabellist(missing)} (the'
                                 f' following labels are available: {fmtlabellist(self.labels)})')
            lblmap = dict( lblmap.items() )
        res = {}
        for k,v in lblmap.items():
            idx = self.__dosidx(k)
            if idx in selected:
                res[k] = v
        return res


def _read_quantumespresso( raw_text_data ):
    assert isinstance(raw_text_data,str)
    _data_lines = raw_text_data.splitlines()

    def _get_header():
        out=[]
        for ll in _data_lines:
            ll = ll.strip()
            if not ll:
                continue
            if not ll.startswith('#'):
                break
            out.append( ' '.join(ll[1:].split()) )
        return out

    hdr = _get_header()
    _expected_hdr = 'Frequency[cm^-1] DOS PDOS'
    if not hdr or _expected_hdr not in hdr:
        from .exceptions import NCBadInput
        raise NCBadInput('Invalid input format. Did not find expected header line "# %s"'%_expected_hdr)
    from ._numpy import _ensure_numpy, _np
    _ensure_numpy()

    rawdata = _np.loadtxt( _data_lines )
    assert len( rawdata.shape ) == 2
    nrows, ncols = rawdata.shape if len( rawdata.shape ) == 2 else ( 0, 0 )
    if not nrows >= 10 or not ncols >= 3:
        from .exceptions import NCBadInput
        raise NCBadInput('Invalid input format. Expected at least 3 columns and 10 rows')
    npdos = ncols - 2

    def _get_col( icol ):
        return rawdata[:,icol]

    egrid = vdos_units_2_eV['1/cm'] * _get_col( 0 )

    doslist  = [ ( 'combined_dos', egrid, _get_col( 1 ) ) ]
    doslist += [ ('pdos%i'%ipdos, egrid, _get_col(2+ipdos) ) for ipdos in range(npdos)]
    return doslist

def _vdos_egrid_is_regular( egrid ):
    if not ( egrid[0] > 0.0 ):
        return False
    from ._common import _grid_is_linspace
    if not _grid_is_linspace( egrid ):
        return False
    bw = ( egrid[-1]-egrid[0] ) / ( len(egrid) - 1 )
    x = egrid[0] / bw
    #regular if x is very near an integer >= 1
    nx = round(x)
    return nx >= 1 and abs(x-nx) < 1e-4

def _do_regularise( egrid, density, n, quiet = False):
    #Regularisation function adapted from ncrystal_vdos2ncmat

    if n >= len(egrid) and _vdos_egrid_is_regular(egrid):
        #Already fine, won't get better just by wasting more points!
        return egrid,density

    from ._numpy import _ensure_numpy, _np, _np_linspace
    _ensure_numpy()

    egrid = _np.asarray(egrid,dtype=float)
    density = _np.asarray(density,dtype=float)
    assert len(egrid)==len(density)
    assert len(egrid) >= 2
    assert egrid[0] > 0.0

    emin,emax=egrid[0],egrid[-1]
    if quiet:
        def nc_print( *a,**kw ):
            pass
    else:
        from ._common import print as nc_print
    nc_print('old range',emin,emax)
    THZ = vdos_units_2_eV['THz']
    nc_print('old range [THz]',emin/THZ,emax/THZ)

    for k in range(1,1000000000):
        #k is number of bins below emin, an integral number by definition in a regularised grid.
        binwidth = emin/k
        nbins=int(_np.floor((emax-emin)/binwidth))+1
        eps = (emin+nbins*binwidth)-emax
        assert eps>=0.0
        if nbins+1 >= n:
            break
    n=nbins+1
    binwidth = emin/k
    new_emax = emin + (n-1) * binwidth
    if abs( (new_emax-binwidth) - emax ) < 1e-3*binwidth:
        nbins -= 1
        n -= 1
        new_emax -= binwidth
    nc_print(f" ==> Choosing regular grid with n={n} pts from emin={emin} to emax={new_emax} ({new_emax-emax} beyond old emax)")
    def retry():
        from ._common import warn
        nnew = n+100
        warn('Something went wrong in DOS regularisation with n={n}. Retrying with n={nnew}.')
        return  _do_regularise(egrid,density,n = nnew, quiet = quiet )
    if not ( new_emax >= emax-binwidth*1.001e-3 ):
        return retry()
    if not ( new_emax >= emax-binwidth*1.001e-3 ):
        return retry()
    new_egrid = _np_linspace(emin,new_emax,n)
    test=new_egrid[0] / ( (new_egrid[-1]-new_egrid[0])/(len(new_egrid)-1) )
    if not abs(round(test)-test)<1e-6:
        return retry()
    new_density = _np.interp(new_egrid,egrid,density, left=0.0, right=0.0)
    nc_print('last density values in new grid:',new_density[-5:])
    return new_egrid,new_density

def _parsevdosunit( name ):
    unitfactor = vdos_units_2_eV.get(name or 'eV')
    if unitfactor:
        return name, unitfactor
    from .exceptions import NCBadInput
    _ = '","'.join(u for u in sorted(vdos_units_2_eV.keys()))
    raise NCBadInput('Invalid frequency unit "%s" (must be one of "%s")'%(name,_))
