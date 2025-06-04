
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

Internal implementation of ncmat2endf.py

"""

from ._numpy import _np
from . import core as nc_core
from . import constants as nc_constants
from . import vdos as nc_vdos
from ._common import print as ncprint
from ._common import warn as ncwarn
from ._common import write_text as ncwrite_text
print = ncprint

#Materials with other processes than the following must be validated by experts
#before conversion to ENDF can be supported:
allowed_scat_proc_names = set( ['ElIncScatter', 'PowderBragg', 'SABScatter'] )

mass_neutron = (nc_constants.const_neutron_mass_amu*
               nc_constants.constant_dalton2eVc2/
               ((nc_constants.constant_c*1e-12)**2)) # eV*ps^2*Angstrom^-2

hbar = nc_constants.constant_planck/nc_constants.k2Pi*1e12 # eV*ps
T0 = 293.6 # K - Reference temperature for LAT=1 in ENDF-6 MF=7/MT=4
ENDF_DESCR_MAXW = 66


_cacheimport=[None]
def import_endfparserpy():
    if _cacheimport[0] is not None:
        return _cacheimport[0]
    try:
        import endf_parserpy
        from endf_parserpy.interpreter.fortran_utils import read_fort_floats
        from endf_parserpy.interpreter.fortran_utils import write_fort_floats
    except ImportError as e:
        ncprint('ERROR: Could not import endf_parserpy.\n\n'
                '       Please check that the endf-parserpy package was'
                ' correctly installed.\n'
                '\n'
                '       Examples of how to do this:\n'
                '\n'
                '          $> conda install -c conda-forge endf_parserpy\n'
                '             (always prefer this if using Conda)\n'
                '\n'
                '          $> pip install ncrystal[endf]\n'
                '             (best way if not using Conda and NCrystal was itself installed via pip)\n'
                '\n'
                '          $> pip install ncrystal[all]\n'
                '             (same as previous but provides even more optional dependencies)\n'
                '\n'
                '          $> pip install endf_parserpy\n'
                '             (might work in other cases if you at least can use pip for dependencies)\n'
                '\n'
                )
        raise ImportError('Could not import endf_parserpy') from e
    _cacheimport[0] = ( endf_parserpy,read_fort_floats,write_fort_floats)
    return _cacheimport[0]

def _endf_roundoff(x):
    #
    # Limit the precision of a float to what can be represented
    # in an ENDF-6 file
    #
    # Receives an  Iterable of float
    #
    _,read_fort_floats,write_fort_floats = import_endfparserpy()
    return _np.array(read_fort_floats( write_fort_floats(x, {'width':11}),
                                       n=len(x),
                                       read_opts={'width':11}),
                     dtype = float)

def _endf_clean(x):
    #
    # Return an array of unique floats that can be represented
    #   in an ENDF-6 file
    #
    # Receives an  Iterable of float
    #
    return _np.unique(_endf_roundoff(x) )

class ElementData():
    #
    # Container for nuclear data for a single element or isotope.
    #
    def __init__(self, ad):
        self._sigma_i = ad.incoherentXS()
        self._sigma_free = ad.freeScatteringXS()
        self._awr = ad.averageMassAMU() / nc_constants.const_neutron_mass_amu
        self._dwi = []
        self._alpha = None
        self._beta = None
        self._beta_total = None
        self._sab_total = []
        self._teff = []
        self._elastic = None
        self._sym = ad.elementName()
        self._zsymam = '{:3d}-'.format(ad.Z()) + self._sym.ljust(2)+'    '
        self._za = ad.Z()*1000

    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, x):
        self._alpha = x

    @property
    def beta(self):
        return self._beta
    @beta.setter
    def beta(self, x):
        self._beta = x

    @property
    def beta_total(self):
        return self._beta_total
    @beta_total.setter
    def beta_total(self, x):
        self._beta_total = x

    @property
    def sab_total(self):
        return self._sab_total

    @property
    def dwi(self):
        return self._dwi
    @dwi.setter
    def dwi(self, x):
        self._dwi = x

    @property
    def teff(self):
        return self._teff

    @property
    def elastic(self):
        return self._elastic
    @elastic.setter
    def elastic(self, x):
        self._elastic = x

    @property
    def awr(self):
        return self._awr

    @property
    def za(self):
        return self._za

    @property
    def zsymam(self):
        return self._zsymam

    @property
    def sigma_free(self):
        return self._sigma_free
    @property
    def sigma_i(self):
        return self._sigma_i
    @sigma_i.setter
    def sigma_i(self, x):
        self._sigma_i = x

class NuclearData():
    #
    # Container for nuclear data for a material.
    #

    # NOTE: The terminology "coh elas" in this class is ENDF-parlance, refering
    #       to Bragg diffraction. Similarly "incoh elas" is used to indicate
    #       incoherent-elastic (as sum of Debye-Waller factors) AND/OR
    #       coherent-elastic approximated via the incoherent approximation to
    #       the shape of incoherent-elastic.

    def __init__(self, *, ncmat_cfg, temperatures, elastic_mode,
                          requested_emax, verbosity ):

        self.__loaded = _decodecfg_and_loadobjs( ncmat_cfg )
        del ncmat_cfg
        self.__di2knlcache = {}
        self._temperatures = tuple(temperatures)
        self._elems = {}
        self._requested_emax = requested_emax
        self._verbosity = verbosity
        self._elastic_mode = elastic_mode
        scattering_components = self.__loaded['scat_comps']
        self._enable_coh_elas = ( 'coh_elas' in scattering_components and
                                  self.__loaded['info_obj'].hasAtomInfo() )
        if not self._enable_coh_elas:
            ncwarn('ENDF output will not contain a Bragg diffraction'
                   ' component')
        self._enable_incoh_elas = ( 'incoh_elas' in scattering_components )
        if not self._enable_incoh_elas:
            ncwarn('ENDF output will not contain any incoherent elastic '
                   'component')
        self._enable_inelas = ( 'inelas' in scattering_components )
        assert self._enable_inelas, ( 'Inelastic component always'
                                     ' must be present for ENDF output')
        # _combine_temperatures:
        # False: use (alpha, beta) grid for lowest temperature
        # True: combine all temperatures
        self._combine_temperatures = False
        for frac, ad in self.composition:
            if not ad.isNaturalElement():
                # TODO: properly handle isolated isotopes and enriched
                #       elements
                from .exceptions import NCBadInput
                raise NCBadInput('Conversion to ENDF is currently supported'
                                 ' only for natural elements')
            sym = ad.elementName()
            if sym in self._elems:
                from .exceptions import NCBadInput
                raise NCBadInput('Conversion to ENDF is not currently supported'
                                 ' for materials where the same element has'
                                 ' multiple dynamic roles (problematic'
                                 f' element: {sym})')

            self._elems[sym] = ElementData(ad)
        if self._enable_coh_elas:
            self._edges = []
            self._sigmaE = []
        else:
            self._edges = None
            self._sigmaE = None
        self._incoherent_fraction = -1
        if len(self.composition) > 1 and elastic_mode == 'scaled':
            #
            # Find element with minimum incoherent contribution.
            #
            self._designated_coherent_atom = None
            for frac, ad in self.composition:
                sym = ad.elementName()
                if (frac/(1.0-frac)*ad.incoherentXS() <
                    self._incoherent_fraction or
                    self._incoherent_fraction == -1):
                    self._incoherent_fraction = (frac/
                                                 (1.0-frac)*ad.incoherentXS())
                    self._designated_coherent_atom = sym
            if self._verbosity > 1:
                ncprint(f'Designated incoherent: {sym}')
        if self._enable_inelas:
            self._get_alpha_beta_grid()
            self._get_inelastic_data()

        self._get_elastic_data()
        self._extract_ncrystal_comments()

    @property
    def comments(self):
        return self._comments

    @property
    def temperatures(self):
        return self._temperatures

    @property
    def ncmat_cfg(self):
        return self.__loaded['cfgstr']

    @property
    def composition(self):
        return self.__loaded['info_obj'].composition

    @property
    def elements(self):
        return self._elems

    @property
    def edges(self):
        return self._edges

    @property
    def sigmaE(self):
        return self._sigmaE

    @property
    def elastic_mode(self):
        return self._elastic_mode

    def _loadKernel( self, di, infoobj_uid ):
        uid = ( infoobj_uid, di.atomIndex )
        cacheobj = self.__di2knlcache.get(uid)
        if cacheobj is not None:
            return cacheobj

        from .vdos import extractKnl
        #Note: using the extractKnl function rather than di.loadKernel means
        #that there will be no reduction of vdoslux for VDOSDebye objects. This
        #is why we use this function here, since ENDF files need the higher
        #energy range of the resulting sab.

        kwargs = dict( vdos = di,
                       mass_amu = di.atomData.averageMassAMU(),
                       temperature = di.temperature,
                       scatxs = 1.0,# seems to be the right thing
                       vdoslux = self.__loaded['vdoslux'] )

        k = extractKnl( target_emax = None, **kwargs )
        emax0 = k.get('suggested_emax',0.0)

        req_emax = self._requested_emax
        if abs(emax0-req_emax) > 1e-9*abs(emax0+req_emax):
            k = extractKnl( target_emax = req_emax, **kwargs )
            lbl = di.atomData.displayLabel()
            ncwarn(f'The extracted kernel for "{lbl}" would normally cover'
                   f' energies up to Emax={emax0:g}eV with the provided'
                   ' cfg-string but a value of'
                   f' Emax={req_emax:g}eV was enforced.')
        self.__di2knlcache[uid] = k
        return k

    def _combine_alpha_beta_grids(self):
        #
        # Combine (alpha, beta) grids from different temperatures
        # into a single grid. This usually results in a huge grid and
        # it is only kept as an option to debug libraries.
        #
        ncwarn('Combining (alpha, beta) grids from different temperatures.'
               ' This usually results in a huge grid.')
        for T in self._temperatures[1:]:
            cfg = self.ncmat_cfg+f';temp={T}K'
            info_obj = nc_core.createInfo(cfg)
            for di in info_obj.dyninfos:
                sym = di.atomData.elementName()
                sctknl = self._loadKernel(di, info_obj.uid )
                self._elems[sym].alpha = _np.unique(_np.concatenate((
                                         self._elems[sym].alpha,
                                         sctknl['alpha']*T/T0)))
                self._elems[sym].beta_total = _np.unique(
                    _np.concatenate(
                        ( self._elems[sym].beta_total, sctknl['beta']*T/T0 )
                    )
                )

    def _get_alpha_beta_grid(self):
        T = self._temperatures[0]
        cfg = self.ncmat_cfg+f';temp={T}K'
        info_obj = nc_core.createInfo(cfg)
        for di in info_obj.dyninfos:
            sym = di.atomData.elementName()
            sctknl = self._loadKernel(di,info_obj.uid)
            self._elems[sym].alpha = sctknl['alpha']*T/T0
            self._elems[sym].beta_total = sctknl['beta']*T/T0
        if self._combine_temperatures:
            self._combine_alpha_beta_grids()

        for frac, ad in self.composition:
            sym = ad.elementName()
            #
            # Remove points that cannot be represented in ENDF data
            #
            self._elems[sym].alpha = _endf_clean(self._elems[sym].alpha)
            self._elems[sym].beta_total = _endf_clean(self._elems[sym].
                                                      beta_total)
            _ = self._elems[sym].beta_total[
                _np.where(self._elems[sym].beta_total<=0)] # get negative beta
            self._elems[sym].beta = -_[::-1] # Invert beta and change sign
            self._elems[sym].beta[0] = 0.0
            if self._verbosity > 2:
                def fmta(x):
                    return '%.4g'%x
                def fmtb(x):
                    return '%.4g'%x
                if unit_test_chop_vals[0]:
                    def fmta(x):
                        return '%.3g'%x
                    def fmtb(x):
                        if x > 50:
                            return '%.2g'%x
                        return '%.3g'%x
                a,b = self._elems[sym].alpha, self._elems[sym].beta
                ncprint(f'>>> alpha points: {len(a)}, alpha range: '
                        f'({fmta(_np.min(a*T0/T))}, {fmta(_np.max(a*T0/T))})')
                ncprint(f'>>> beta points: {len(b)}, beta range: '
                        f'({fmtb(_np.min(b*T0/T))}, {fmtb(_np.max(b*T0/T))})')

    def _get_coherent_elastic(self, T):
        #
        # Load coherent elastic data (Bragg diffraction)
        #
        cfg = (self.ncmat_cfg+f';temp={T}K;comp=bragg')
        m = nc_core.load(cfg)
        if T == self._temperatures[0]:
            # Find unique Bragg edges, as represented in ENDF-6 floats
            edges = _np.array([nc_constants.wl2ekin(2.0*e.dspacing)
                              for e in m.info.hklObjects()])
            self._edges = _endf_clean(edges)
            # Coherent scattering XS is evaluated between edges
            eps = 1e-3
            self._evalpoints = _np.concatenate((self._edges[:-1]**(1-eps)*
                               self._edges[1:]**eps, [self._edges[-1]]))
        sigmaE = _endf_roundoff(m.scatter.xsect(self._evalpoints)
                                *self._evalpoints)
        assert _np.all(sigmaE[:-1] <= sigmaE[1:]),\
               'Sigma*E in Bragg edges not cummulative'
        self._sigmaE.append(sigmaE)

    def _get_elastic_data(self):
        for T in self._temperatures:
            if self._enable_coh_elas:
                self._get_coherent_elastic(T)
            if self._enable_incoh_elas:
                cfg = self.ncmat_cfg+f';temp={T}K'
                info_obj = nc_core.createInfo(cfg)
                for di in info_obj.dyninfos:
                    sym = di.atomData.elementName()
                    emin = di.vdosData()[0][0]
                    emax = di.vdosData()[0][1]
                    rho = di.vdosData()[1]
                    res = nc_vdos.analyseVDOS(emin, emax, rho, di.temperature,
                                              di.atomData.averageMassAMU())
                    #
                    # Load incoherent elastic data
                    #
                    msd = res['msd']
                    self._elems[sym].dwi.append(
                        _tidy_teffwp( msd*2*mass_neutron/hbar**2 )
                    )
            else:
                for frac, ad in self.composition:
                    sym = ad.elementName()
                    self._elems[sym].sigma_i =  None
                    self._elems[sym].dwi =  None
        if self._verbosity > 1:
            ncprint('>> Prepare elastic approximations')
        if ( self._elastic_mode == 'scaled' ):
            if ( len(self.composition) > 1 and
                 self._incoherent_fraction < 1e-6 ):
                    self._elastic_mode = 'greater'
                    ncwarn('Scaled elastic mode requested '
                           'but all elements are coherent. '
                           '"greater" option will be used instead.')
        for frac, ad in self.composition:
            sym = ad.elementName()
            if ( self._sigmaE is None and
                 self._elems[sym].dwi is None ):
                ncwarn('No coherent elastic data or '
                       'incoherent elastic data found: '
                       'only inelastic data will be written')
                self._elems[sym].elastic = None
                continue
            if self._elastic_mode == 'mixed': # iel = 100
                if (self._sigmaE is None):
                    # mixed elastic requested but only incoherent available
                    ncwarn(f'Mixed elastic mode for {sym} but no '
                            'Bragg edges found: incoherent approximation')
                    self._elems[sym].sigma_i = (ad.incoherentXS() +
                                                ad.coherentXS())
                    self._elems[sym].elastic = 'incoherent'
                elif (self._elems[sym].dwi is None):
                    ncwarn(f'Mixed elastic mode for {sym} but no '
                            'incoherent elastic data found: '
                            'coherent approximation')
                    self._elems[sym].elastic = 'coherent'
                else:
                    if self._verbosity>1:
                        ncprint(f'>> Mixed elastic mode for {sym}')
                    self._elems[sym].elastic = 'mixed'
            elif self._elastic_mode == 'greater': # iel = 98
                if ((ad.incoherentXS() > ad.coherentXS()) or
                    (self._sigmaE is None)):
                    if self._verbosity>1:
                        ncprint( '>> Principal elastic mode '
                                f'for {sym}: incoherent')
                    self._edges = None
                    self._sigmaE = None
                    self._elems[sym].elastic = 'incoherent'
                else:
                    if self._verbosity>1:
                        ncprint(f'>> Principal elastic mode for {sym}: coherent')
                    self._elems[sym].sigma_i =  None
                    self._elems[sym].dwi =  None
                    self._elems[sym].elastic = 'coherent'
            elif self._elastic_mode == 'scaled':
                if len(self.composition) == 1: # iel = 99, single atomic case
                    if ((ad.incoherentXS() > ad.coherentXS()) or
                        (self._sigmaE is None)):
                        if self._verbosity>1:
                            ncprint( '>> Scaled elastic mode for '
                                    f'single atom {sym}: incoherent')
                        self._edges = None
                        self._sigmaE = None
                        self._elems[sym].sigma_i = (ad.incoherentXS() +
                                                    ad.coherentXS())
                        self._elems[sym].elastic = 'incoherent'
                    else:
                        if self._verbosity>1:
                            ncprint( '>> Scaled elastic mode for '
                                    f'single atom {sym}: coherent')
                        self._elems[sym].sigma_i =  None
                        self._elems[sym].dwi =  None
                        self._elems[sym].elastic = 'coherent'
                        self._sigmaE = [x*(ad.incoherentXS() +
                                        ad.coherentXS())/ad.coherentXS()
                                        for x in self._sigmaE]
                elif (self._sigmaE is None):
                    # iel = 99, incoherent approximation
                    if self._verbosity>1:
                        ncprint(f'>> Scaled elastic mode for {sym}: '
                               'incoherent approximation')
                    self._elems[sym].sigma_i = (ad.incoherentXS() +
                                                ad.coherentXS())
                    self._elems[sym].elastic = 'incoherent'
                else:
                    # iel = 99, multi atomic case
                    if sym == self._designated_coherent_atom:
                        if self._verbosity>1:
                            ncprint(f'>> Scaled elastic mode for {sym} in '
                                     'compound: designated coherent atom, '
                                    f'dividing by frac^2={frac**2}')
                        self._elems[sym].sigma_i =  None
                        self._elems[sym].dwi =  None
                        self._elems[sym].elastic = 'coherent'
                        self._sigmaE = [x/frac for x in self._sigmaE]
                    else:
                        if self._verbosity>1:
                            ncprint(f'>> Scaled elastic mode for {sym} '
                                     'in compound: incoherent')
                        self._elems[sym].elastic = 'incoherent'
                        self._elems[sym].sigma_i = ((1.0+
                                                     self._incoherent_fraction/
                                                     ad.incoherentXS())*
                                                    self._elems[sym].sigma_i)

    def _get_inelastic_data(self):
        for T in self._temperatures:
            cfg = self.ncmat_cfg+f';temp={T}K'
            info_obj = nc_core.createInfo(cfg)
            for di in info_obj.dyninfos:
                sym = di.atomData.elementName()
                #
                # Load incoherent inelastic data
                #
                sctknl = self._loadKernel(di,info_obj.uid)
                if self._verbosity > 2:
                    ncprint(f'>>> Interpolating T={T}K for {sym}')

                alpha = sctknl['alpha']
                beta = sctknl['beta']
                sab = sctknl['sab']
                #
                # Interpolate S(a,b) because the NCrystal grid might
                # contain numbers that cannot be represented
                # as distinct FORTRAN reals in the ENDF-6 file
                #
                sab.shape = (beta.size, alpha.size)
                sab_int = _interp2d(self._elems[sym].alpha*T0/T,
                          self._elems[sym].beta_total*T0/T,
                          alpha, beta, sab.transpose()).transpose()
                self._elems[sym].sab_total.append(sab_int)
                emin = di.vdosData()[0][0]
                emax = di.vdosData()[0][1]
                rho = di.vdosData()[1]
                res = nc_vdos.analyseVDOS(emin, emax, rho, di.temperature,
                                          di.atomData.averageMassAMU())
                self._elems[sym].teff.append( _tidy_teffwp( res['teff'] ) )

    def _extract_ncrystal_comments( self ):
        # TODO: handle multi phase materials
        ncmat_fn = self.__loaded['cfgstr_decoded']['data_name']
        td = nc_core.createTextData(ncmat_fn)
        from ._ncmatimpl import _extractInitialHeaderCommentsFromNCMATData
        from textwrap import wrap
        raw = _extractInitialHeaderCommentsFromNCMATData(td)

        #We need to rewrap for ENDF_DESCR_MAXW. The following gymnastics allow
        #us to identify each "paragraph" of text:

        paragraphs = []
        current_pg = []
        for line in raw:
            line = line.strip()
            if not line:
                if not current_pg:
                    continue#ignore leading or multiple empty lines
                paragraphs.append( ' '.join(current_pg) )
                current_pg = []
                continue

            # Detect and shorten horisontal ascii art lines like the following,
            # and give them their own paragraph entry:
            #
            #  --------------------------------------------
            #  ++++++++++++ Stuff here ++++++++++++++++++++
            #
            is_hr, line =  _detect_and_short_horisontal_ruler_in_line( line )

            if is_hr:
                #HR. Flush previous paragraph and add hr line as its own:
                if current_pg:
                    paragraphs.append( ' '.join(current_pg) )
                    current_pg = []
                paragraphs.append( line )

            else:
                #Just add the line to the current paragraph:
                current_pg.append( line )

        if current_pg:
            paragraphs.append( ' '.join(current_pg) )

        wrapped_pgs = []
        for pg in paragraphs:
            if len(pg) > ENDF_DESCR_MAXW:
                pg = '\n'.join(wrap( pg,
                                     width=ENDF_DESCR_MAXW,
                                     break_long_words=False ))
            assert isinstance(pg,str)
            wrapped_pgs.append( pg )

        res = []
        for pg in wrapped_pgs:
            if res:
                res.append('')
            res += pg.splitlines()
        self._comments = res



class EndfFile():
    #
    # Container for data for a therma ENDF file.
    # Includes a write() method to create the file using endf-parserpy
    #
    def __init__(self, element, data, mat, endf_metadata, *,
                 include_gif=False, isotopic_expansion=False,
                 smin=None, emax=None, lasym=None, verbosity=1):
        endf_parserpy,_,_ = import_endfparserpy()

        self._endf_parserpy_version = endf_parserpy.__version__
        self._sym = element
        self._data = data
        self._endf_metadata = endf_metadata
        self._mat = mat
        self._include_gif = include_gif
        assert not isotopic_expansion, "isotopic_expansion not supported yet"
        self._isotopic_expansion = isotopic_expansion
        self._verbosity = verbosity
        assert smin is not None, 'smin not set'
        self._smin = smin
        assert emax is not None, 'emax not set'
        self._emax = emax
        assert lasym is not None, 'lasym not set'
        self._lasym = lasym
        self._endf_dict = endf_parserpy.EndfDict()
        self._endf_dict['0/0'] = {}
        self._endf_dict['0/0']['MAT'] = self._mat
        self._endf_dict['0/0']['TAPEDESCR'] = 'Created with ncmat2endf'
        self._createMF1()
        self._createMF7()
        self._parser = None

    def _createMF7MT2(self, elastic):
        #
        # Creates endf-parserpy dictionary for MF=7 file,
        # MT=2 reaction (thermal elastic) of a thermal ENDF file
        #
        mat = self._mat
        data = self._data
        awr = data.elements[self._sym].awr
        za = data.elements[self._sym].za
        temperatures = data.temperatures
        self._endf_dict['7/2'] = {}
        d = self._endf_dict['7/2']
        d['MAT'] = mat
        d['ZA'] = za
        d['AWR'] = awr
        if elastic in ['coherent', 'mixed']:
            edges  = data.edges
            sigmaE = data.sigmaE[0]
            d['T0'] = temperatures[0]
            d['LT'] = len(temperatures)-1
            d['S_T0_table'] = {}
            d['S_T0_table']['NBT'] = [len(edges)]
            d['S_T0_table']['INT'] = [1]
            d['S_T0_table']['Eint'] = edges.tolist()
            d['S_T0_table']['S'] = sigmaE.tolist()
            d['T'] = {k: v for k, v in
                      enumerate(temperatures[1:], start=1)}
            S = {}
            for q, v in enumerate(edges, start=1):
                S[q] = {}
                for i,v in enumerate(temperatures[1:], start=1):
                    S[q][i] = float( data.sigmaE[i][q-1] )
            d['S'] = S
            d['LI'] = 2
            d['NP'] = len(edges)
        if elastic in ['incoherent', 'mixed']:
            d['SB'] = data.elements[self._sym].sigma_i
            d['Wp']  = data.elements[self._sym].dwi
            d['Tint'] = list(temperatures)
            d['INT'] = [2]
            d['NBT'] = [len(temperatures)]
        lthr_values = {'coherent':1, 'incoherent':2, 'mixed':3}
        d['LTHR'] = lthr_values[data.elements[self._sym].elastic]

    def _createMF7MT4(self):
        #
        # Creates endf-parserpy dictionary for MF=7 file,
        # MT=4 reaction (thermal inelastic) of a thermal ENDF file
        #
        mat = self._mat
        data = self._data
        awr = data.elements[self._sym].awr
        za = data.elements[self._sym].za
        temperatures = data.temperatures
        self._endf_dict['7/4'] = {}
        d = self._endf_dict['7/4']
        d['MAT'] = mat
        d['ZA'] = za
        d['AWR'] = awr
        d['LAT'] =  1 # (alpha, beta) grid written for 0.0253 eV
        d['LASYM'] = self._lasym # symmetric/asymmetric S(a,b)
        d['LLN'] = 0 # linear S is stored
        d['NI'] = 6
        d['NS'] = 0
        d['B'] = {1:data.elements[self._sym].sigma_free,
                  2:self._emax,
                  3:awr,
                  4:self._emax,
                  5:0,                                # unused
                  6:1                                 # natom
                 }
        alpha = data.elements[self._sym].alpha
        beta = data.elements[self._sym].beta
        beta_total = data.elements[self._sym].beta_total
        sab_data = []
        for sab_total, T in zip(data.elements[self._sym].sab_total,
                                temperatures):
            _, beta_grid = _np.meshgrid(alpha*T0/T, beta_total*T0/T)
            if (self._lasym == 0) or (self._lasym == 1):
                detailed_balance_factor = _np.exp(beta_grid/2)
            if self._lasym == 3:
                # S(a,b) for all beta
                sab_data.append(sab_total.transpose())
                continue
            if self._lasym == 2:
                # S(a,b) for negative beta
                # get negative branch of S(a,b)
                sab_sym2 = sab_total[_np.where(beta_grid<=0)]
                sab_sym2.shape = (len(beta), len(alpha))
                sab_sym3 = sab_sym2[::-1,:]  # Invert S(a,b) for negative beta
                sab_data.append(sab_sym3.transpose())
                continue
            if self._lasym == 1:
                # S(a,b)*exp(-b/2) for all beta
                sab_sym = sab_total*detailed_balance_factor
                sab_data.append(sab_sym.transpose())
                continue
            if self._lasym == 0:
                # S(a,b)*exp(-b/2) for negative beta
                sab_sym = sab_total*detailed_balance_factor
                # get negative branch of S(a,b)
                sab_sym2 = sab_sym[_np.where(beta_grid<=0)]
                sab_sym2.shape = (len(beta), len(alpha))
                # Invert S(a,b) for negative beta
                sab_sym3 = sab_sym2[::-1,:]
                sab_data.append(sab_sym3.transpose())
                continue

        if (self._lasym == 0) or (self._lasym == 2):
            # Save S(a,b) or S(a,b)*exp(-b/2) for negative beta
            d['NB'] = len(beta)
            d['beta_interp/NBT'] = [len(beta)]
            d['beta_interp/INT'] = [4]

            d['T0'] = temperatures[0]
            d['beta'] = {k:_tidy_beta(v) for k, v in enumerate(beta,start=1)}
            d['LT'] = {k:len(temperatures)-1
                       for k, v in enumerate(beta,start=1)}
            d['T'] = {k:v for k,v in enumerate(temperatures[1:], start=1)}
            d['LI'] = {k:4 for k,v in enumerate(temperatures[1:], start=1)}
            d['NP'] = len(alpha)
            S1 = {}
            sab = sab_data[0]
            sab[sab < self._smin] = 0.0
            for j,v in enumerate(beta, start=1):
                S1[j] = {}
                S1[j]['NBT'] = [len(alpha)]
                S1[j]['INT'] = [4]
                S1[j]['alpha'] = _tidy_alpha_list( (alpha/awr).tolist() )
                S1[j]['S'] = _tidy_sab_list(sab[:,j-1].tolist())
            d['S_table'] = S1

            S2 = {}
            if len(temperatures) > 1:
                for q, v in enumerate(alpha, start=1):
                    S2[q] = {}
                    for j,v in enumerate(beta, start=1):
                        sab = []
                        for i,v in enumerate(temperatures[1:], start=1):
                                sval = sab_data[i][q-1,j-1]
                                if sval < self._smin:
                                    sval = 0.0
                                sab.append(sval)
                        sab = _tidy_sab_list(sab)
                        S2[q][j] = {k: float(v) for k,v in enumerate(sab, start =1)}
            d['S'] = S2
        elif (self._lasym == 1) or (self._lasym == 3):
            # Save S(a,b) or S(a,b)*exp(-b/2) for all beta
            alpha = data.elements[self._sym].alpha
            beta = data.elements[self._sym].beta_total

            d['NB'] = len(beta)
            d['beta_interp/NBT'] = [len(beta)]
            d['beta_interp/INT'] = [4]

            d['T0'] = temperatures[0]
            d['beta'] = {k:_tidy_beta(v, allow_negative=True)
                         for k, v in enumerate(beta,start=1)}
            d['LT'] = {k:len(temperatures)-1
                       for k, v in enumerate(beta,start=1)}
            d['T'] = {k:v for k,v in enumerate(temperatures[1:], start=1)}
            d['LI'] = {k:4 for k,v in enumerate(temperatures[1:], start=1)}
            d['NP'] = len(alpha)
            S1 = {}
            sab = sab_data[0]
            sab[sab < self._smin] = 0.0
            for j,v in enumerate(beta, start=1):
                S1[j] = {}
                S1[j]['NBT'] = [len(alpha)]
                S1[j]['INT'] = [4]
                S1[j]['alpha'] = _tidy_alpha_list((alpha/awr).tolist())
                S1[j]['S'] = _tidy_sab_list(sab[:,j-1].tolist())
            d['S_table'] = S1

            S2 = {}
            if len(temperatures) > 1:
                for q, v in enumerate(alpha, start=1):
                    S2[q] = {}
                    for j,v in enumerate(beta, start=1):
                        sab = []
                        for i,v in enumerate(temperatures[1:], start=1):
                                sval = sab_data[i][q-1,j-1]
                                if sval < self._smin:
                                    sval = 0.0
                                sab.append(sval)
                        sab = _tidy_sab_list(sab)
                        S2[q][j] = {k: float(v) for k,v in enumerate(sab, start =1)}
            d['S'] = S2

        d['teff0_table/Teff0'] = data.elements[self._sym].teff
        d['teff0_table/Tint'] = list(temperatures)
        d['teff0_table/NBT'] = [len(temperatures)]
        d['teff0_table/INT'] = [2]

    def _createMF7(self):
        #
        # Creates MF=7 file of a thermal ENDF file.
        #  See ENDF-102, sect. 7.
        #  https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf
        #
        if self._verbosity > 1:
            ncprint('> Generate MF7')
        mat = self._mat
        data = self._data
        awr = data.elements[self._sym].awr
        za = data.elements[self._sym].za
        #
        # Prepare dictionary for thermal elastic reaction
        #
        elastic = data.elements[self._sym].elastic
        if elastic is not None:
            self._createMF7MT2(elastic)
        #
        # Prepare dictionary for thermal inelastic reaction
        #
        self._createMF7MT4()

        if self._include_gif:
            assert not self._isotopic_expansion
            self._endf_dict['7/451'] = {}
            d = self._endf_dict['7/451']
            d['MAT'] = mat
            d['ZA'] = za
            d['AWR'] = awr
            d['NA'] = 1
            d['NAS'] = 1
            d['NI'] = {1:1}
            d['ZAI'] = {1:{1:za}}
            d['LISI'] = {1:{1:0}}
            d['AFI'] = {1:{1:1.0}}
            d['SFI'] = {1:{1:data.elements[self._sym].sigma_free}}
            d['AWRI'] = {1:{1:awr}}

    def _createMF1MT451description(self, data, endf_metadata):
        desc = []
        desc.append(ENDF_DESCR_MAXW*'*')
        desc.append('')
        desc.append(' This file was converted from the following NCrystal [1,2]')
        desc.append(' cfg-string:')
        desc.append('')

        if len(data.ncmat_cfg) < ENDF_DESCR_MAXW-6:
            desc.append(f'"{data.ncmat_cfg}"'.center(ENDF_DESCR_MAXW))
        else:
            cfg_parts = data.ncmat_cfg.split(';')
            cfg_lines = []
            from textwrap import wrap
            for part in cfg_parts:
                cfg_lines += wrap(part.lstrip(),
                                  width=60,
                                  break_long_words=True)
                cfg_lines[-1] += ';'
            desc.append(f'  "{cfg_lines[0]}')
            desc.extend(f'   {_}' for _ in cfg_lines[1:-1])
            desc.append(f'   {cfg_lines[-1][:-1]}"')

        desc.append('')
        nc_version = nc_core.get_version()
        ep_version = self._endf_parserpy_version
        if unit_test_chop_vals[0] or unit_test_not_write_version[0]:
            nc_version = 'NCVERSION'
            ep_version = 'EPVERSION'
        desc.append(f' using NCrystal {nc_version} ')
        desc.append(f' and endf-parserpy [3] {ep_version} ')
        desc.append(' with the following options:')
        desc.append('')
        desc.append(f'  smin:{self._smin}')
        desc.append(f'  emax:{self._emax}')
        desc.append(f'  lasym:{self._lasym}')
        desc.append(f'  include_gif:{self._include_gif}')
        desc.append(f'  isotopic_expansion:{self._isotopic_expansion}')
        desc.append(f'  elastic_mode:{data.elastic_mode}')
        desc.append('')
        desc.append(' Temperatures:')
        for T in data.temperatures:
            desc.append(f'       {T:.2f} K')
        desc.append('')
        desc.append('References:')
        desc.append('[1] https://github.com/mctools/ncrystal')
        desc.append('[2] https://doi.org/10.1016/j.cpc.2019.07.015')
        desc.append('[3] https://endf-parserpy.readthedocs.io/en/latest/')
        desc.append('')
        desc.append(ENDF_DESCR_MAXW*'*')
        desc.append('')
        desc.append('Comments from NCMAT file:')
        desc.append('')
        for line in data.comments:
            desc.append(line)
        desc.append(ENDF_DESCR_MAXW*'*')

        if any( len(line) > ENDF_DESCR_MAXW for line in desc ):
            longmarker = '+LINECONT+'
            ncwarn( 'Desciption contains lines too long for ENDF. Will break '
                    f'lines and use "{longmarker}" to mark line continuations' )
            newdesc = []
            for line in desc:
                if len(line)<=ENDF_DESCR_MAXW:
                    newdesc.append( line )
                    continue
                newdesc.append( line[:ENDF_DESCR_MAXW] )
                s = line[ENDF_DESCR_MAXW:]
                while s:
                    s = longmarker + s
                    newdesc.append( s[:ENDF_DESCR_MAXW] )
                    s = s[ENDF_DESCR_MAXW:]
            desc = newdesc
        assert max( len(line) for line in desc) <= ENDF_DESCR_MAXW
        return [ line.ljust(ENDF_DESCR_MAXW) for line in desc ]

    def _createMF1(self):
        #
        # Creates MF=1 file of a thermal ENDF file.
        #    See ENDF-102, sect. 1.
        #    https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf
        #
        mat = self._mat
        data = self._data
        awr = data.elements[self._sym].awr
        endf_metadata = self._endf_metadata
        za = data.elements[self._sym].za
        zsymam = data.elements[self._sym].zsymam
        self._endf_dict['1/451'] = {}
        d = self._endf_dict['1/451']
        d['MAT'] = mat
        d['ZA'] = za
        d['AWR'] = awr
        d['LRP'] = -1
        d['LFI'] = 0
        d['NLIB'] = endf_metadata.nlib
        d['NMOD'] = 0
        d['ELIS'] = 0
        d['LIS'] = 0
        d['LISO'] = 0
        d['STA'] = 0
        d['NFOR'] = 6
        d['AWI'] = 1.0
        d['EMAX'] = self._emax
        d['LREL'] = endf_metadata.lrel
        d['NSUB'] = 12
        d['NVER'] = endf_metadata.nver
        d['TEMP'] = 0.0
        d['LDRV'] = 0
        d['HSUB/1'] = f'----{endf_metadata.libname:18s}MATERIAL {mat:4d}'
        d['HSUB/1'].ljust(ENDF_DESCR_MAXW)
        d['HSUB/2'] =  '-----THERMAL NEUTRON SCATTERING DATA'.ljust(ENDF_DESCR_MAXW)
        d['HSUB/3'] =  '------ENDF-6 FORMAT'.ljust(ENDF_DESCR_MAXW)
        d['NXC'] = 1
        d['ZSYMAM'] = zsymam.ljust(11)
        d['ALAB'] = endf_metadata.alab.ljust(11)
        d['AUTH'] = endf_metadata.auth.ljust(33)
        d['REF'] = endf_metadata.reference.ljust(21)
        d['EDATE'] = ( 'EVAL-MMMYY' if endf_metadata.edate is None else
                                      'EVAL-'+endf_metadata.edate )
        d['DDATE'] = ( 'DIST-MMMYY' if endf_metadata.ddate is None else
                                      'DIST-'+endf_metadata.ddate )
        d['RDATE'] = ( 'REV0-MMMYY' if endf_metadata.rdate is None else
                                      f'REV{endf_metadata.lrel:1d}-'+
                                      endf_metadata.rdate )
        d['ENDATE'] = endf_metadata.endate.ljust(8)
        desc = self._createMF1MT451description(data, endf_metadata)
        d['DESCRIPTION'] = {k:v for k, v in enumerate(desc, start=1)}
        d['NWD'] = 5+len(desc)
        d['MFx/1'] = 1
        d['MTx/1'] = 451
        d['NCx/1'] = 5
        d['MOD/1'] = d['NMOD']

    def write(self, endf_fn, *, force, outdir ):
        import pathlib
        outfile = pathlib.Path(outdir).joinpath( pathlib.Path(endf_fn) )
        outfile.parent.mkdir(exist_ok = True, parents = True)

        if self._verbosity > 0:
            ncprint(f'Write ENDF file {outfile.name} ...')
        if outfile.exists() and not force:
            from .exceptions import NCBadInput
            raise NCBadInput('Error: output file already exists'
                             ' (run with force=True or --force to overwrite)')
        assert outfile.parent.is_dir()
        if unit_test_abort_write[0]:
            if unit_test_abort_write[0] == 'dump':
                self.dump_endf_dict()
            return

        if self._parser is None:
            endf_parserpy,_,_ = import_endfparserpy()
            self._parser = endf_parserpy.EndfParser(
                explain_missing_variable=True,
                cache_dir=False
            )
            endf_parserpy.update_directory(self._endf_dict, self._parser)

        text = '\n'.join(self._parser.write(self._endf_dict,
                                            zero_as_blank=True))
        ncwrite_text(outfile,text)

    def dump_endf_dict( self ):
        ncprint('DUMPING endf dict begin')
        _dump_dict(self._endf_dict,prefix='  ')
        ncprint('DUMPING endf dict end')


def _decodecfg_and_loadobjs( cfgstr ):
    #We want to load Info and Scatter objects, but also normalise the cfg-string
    #and extract related items like vdoslux and scattering components. However,
    #with the caveat that we want the temperature to be contained in the final
    #cfg-string, and emit a warning if it was absent in the initial cfgstr. Will
    #throw exceptions in case of multiphase materials.

    def is_already_loaded_obj( x ):
        if not isinstance(x,dict):
            return False
        expected = set(['info_obj','scat_obj','cfgstr','cfgstr_decoded','temp',
                        'scat_comps','vdoslux'])
        return set(x.keys()) == expected

    if is_already_loaded_obj( cfgstr ):
        return cfgstr #nothing to do, already loaded

    from .cfgstr import normaliseCfg, decodeCfg, decodecfg_vdoslux
    from .misc import detect_scattering_components

    def cfg_has_explicit_temp( cstr ):
        # Due to the possibility of cfg-strings embedded in NCMAT data, we use
        # the following manner of detection:
        assert '>' not in cstr, "logic error (unecpected multiphase syntax)"
        return any ( (''.join(e.strip().split())).startswith('temp=')
                    for e in cstr.split(';')[1:] )

    def _decode_cfg( cstr ):
        dc = decodeCfg(cstr)
        assert dc['format'] == 'NCrystal-MatCfg-v1'
        return dc

    def fmtfp( v ):
        s = '%.14g'%v
        return s if float(s) == v else '%.17g'%v

    multiphase_errmsg = 'Only single phase materials supported'
    if _decode_cfg(cfgstr)['ismultiphase']:
        from .exceptions import NCBadInput
        raise NCBadInput(multiphase_errmsg)

    info0 = nc_core.createInfo(cfgstr)
    temp = info0.getTemperature()
    del info0

    warn_missing_explicit_temp = not cfg_has_explicit_temp( cfgstr )

    norm_cfg = normaliseCfg( cfgstr )
    del cfgstr
    if not cfg_has_explicit_temp( norm_cfg ):
        norm_cfg += f';temp={fmtfp(temp)}'

    d = {}
    d['info_obj'] = nc_core.createInfo(norm_cfg)
    if not d['info_obj'].isSinglePhase():
        from .exceptions import NCBadInput
        raise NCBadInput(multiphase_errmsg)
    d['scat_obj'] = nc_core.createScatter(norm_cfg)
    d['cfgstr'] = norm_cfg
    d['cfgstr_decoded'] = _decode_cfg(norm_cfg)
    d['temp'] = temp
    d['scat_comps'] = detect_scattering_components( norm_cfg )
    d['vdoslux'] = decodecfg_vdoslux( norm_cfg )

    if warn_missing_explicit_temp:
        ncwarn( 'Temperature not explicitly given in the cfg-string,'
                f' using T={fmtfp(temp)}K')

    assert is_already_loaded_obj( d )
    return d

def _impl_ncmat2endf( *,
                      ncmat_cfg,
                      material_name,
                      endf_metadata,
                      othertemps,
                      elastic_mode,
                      include_gif,
                      isotopic_expansion,
                      force,
                      smin,
                      emax,
                      lasym,
                      outdir,
                      verbosity ):
    from .exceptions import NCBadInput
    from . import core as nc_core
    from ._common import warn as ncwarn
    from ._common import print as ncprint
    from ._numpy import _ensure_numpy
    from .ncmat2endf import ( EndfMetaData,
                              available_elastic_modes )
    _ensure_numpy()

    if not isinstance(verbosity,int) or not ( 0<=verbosity<=3):
        raise NCBadInput('Invalid value of verbosity parameter (expects '
                         f'value 0, 1, 2, or 3): {verbosity}')

    loaded = _decodecfg_and_loadobjs( ncmat_cfg )
    del ncmat_cfg


    if not endf_metadata:
        endf_metadata = EndfMetaData()
    elif not isinstance(endf_metadata,EndfMetaData):
        _ = EndfMetaData()
        _.update_from_dict(endf_metadata)
        endf_metadata = _

    if elastic_mode not in available_elastic_modes:
        raise NCBadInput(f'Elastic mode {repr(elastic_mode)}'
                         f' not in ({available_elastic_modes})')
    if lasym > 0:
        ncwarn( 'Creating non standard S(a,b)'
               f' with LASYM = {lasym}')
    if loaded['scat_obj'].isOriented():
        raise NCBadInput('Oriented materials cannot be represented in the'
                         ' ENDF format and are not supported' )
    supported_comps = ['inelas','coh_elas','incoh_elas']
    unsupported_comps = set(loaded['scat_comps']) - set(supported_comps)
    if unsupported_comps:
        c = list(unsupported_comps)[0]
        raise NCBadInput(f'Materials with scattering component "{c}" can not'
                         ' be represented in the ENDF format')

    if 'inelas' not in loaded['scat_comps']:
        raise NCBadInput('MF7/MT4 is mandatory in an ENDF file '
                         'but no inelastic data found' )

    base_temp = loaded['temp']
    if othertemps is None:
        othertemps = tuple()
    else:
        if type(othertemps) in (int, float):
            othertemps = (othertemps,)
        elif type(othertemps) in (list, tuple):
            if any(type(T) not in (int, float) for T in othertemps):
                raise NCBadInput('Something wrong with the othertemps'
                                 f' parameter: ({othertemps})')
            else:
                othertemps = tuple(float( T ) for T in othertemps )
        else:
            raise NCBadInput('othertemps parameter: should be a list or tuple '
                             'of float or int')
    if base_temp in othertemps:
        raise NCBadInput('Repeated temperatures: othertemps parameter must not '
                         'include the temperature defined in the cfg string')
    temperatures = sorted(othertemps + (base_temp,))
    if len(temperatures) > 1:
        ncwarn('Multiple temperatures requested. Although this is supported, '
               'it is not recommended because NCrystal generates '
               'a custom (alpha,beta) grid for each temperature. '
               'Thus the (alpha,beta) grid for the first temperature will '
               'be used, and S(alpha, beta) for other temperatures '
               'will be interpolated onto it.')

    if any( T<=0 for T in temperatures ):
        raise NCBadInput('Non positive temperatures')

    for di in loaded['info_obj'].dyninfos:
        if type(di) not in (nc_core.Info.DI_VDOS, nc_core.Info.DI_VDOSDebye):
            raise NCBadInput('Conversion to ENDF supported only '
                             'for VDOS and VDOSDebye dyninfos')

    if (isotopic_expansion and not include_gif):
        raise NCBadInput( 'Isotopic expansion requires generalized information'
                          ' file, use --gif' )
    if isotopic_expansion:
        raise NCBadInput('Isotopic expansion in conversion to ENDF is not'
                         ' yet supported')

    if loaded['scat_obj'].isNull():
        raise NCBadInput('Material configuration indicates no active'
                         ' scattering processes to convert')

    scat_proc_names = set( _get_scat_proc_names(loaded['scat_obj']) )
    unsupported_scat_procs = scat_proc_names - allowed_scat_proc_names
    if unsupported_scat_procs:
        raise NCBadInput('Material configuration indicates scattering'
                         ' processes which has not been vetted for conversion'
                         ' to the ENDF format'
                         ': "%s"'%('", "'.join(sorted(unsupported_scat_procs))))

    if verbosity > 0:
        ncprint('Initialise nuclear data...')

    data = NuclearData(ncmat_cfg=loaded,
                       temperatures=temperatures,
                       elastic_mode=elastic_mode,
                       requested_emax=emax,
                       verbosity=verbosity)

    if endf_metadata.matnum is not None:
        n = len(endf_metadata.matnum)
        for frac, ad in data.composition:
            if ad.elementName() in endf_metadata.matnum.keys():
                n = n - 1
        if n != 0:
            raise NCBadInput('Incorrect material number assignment')

    if material_name is None:
        #Default depends on whether or not it is a monoatomic material:
        material_name = ( None
                          if len(data.composition) == 1
                          else 'UnknownCompound' )

    output_composition = []
    for frac, ad in data.composition:
        sym = ad.elementName()
        mat = ( 999 if not endf_metadata.matnum
                else endf_metadata.matnum.get(sym))
        assert mat is not None, ('Incorrect material number '
                                 f'assignment for symbol "{sym}"')
        endf_fn = ( f'tsl_{sym}.endf'
                    if not material_name
                    else f'tsl_{sym}_in_{material_name}.endf' )
        assert data.elements[sym].sab_total is not None, ('Scattering kernel'
                                            f' not available for: {endf_fn}')
        endf_file = EndfFile(sym, data, mat, endf_metadata,
                             include_gif=include_gif,
                             isotopic_expansion=isotopic_expansion,
                             smin=smin, emax=emax, lasym=lasym,
                             verbosity=verbosity)
        endf_file.write( endf_fn, force = force, outdir = outdir )
        output_composition.append( (endf_fn, frac, sym) )

    if verbosity > 0:
        ncprint('Files created:')
        from ._common import prettyFmtValue as pfmt
        for fn, frac, sym in output_composition:
            ncprint(f'  {fn} : {sym} with fraction {pfmt(frac)}')
        if len(temperatures)==1:
            ncprint('Suggested material density: '
                    '%.10g g/cm^3'%loaded['info_obj'].density)

    return [ (fn,frac) for fn,frac,_ in output_composition ]


_metadata_definitions = dict(
    ALAB = dict( defval = 'MyLAB' ),
    AUTH = dict( defval = 'NCrystal' ),
    LIBNAME = dict( defval = 'MyLib' ),
    NLIB = dict( datatype = int, defval = 0 ),
    REFERENCE = dict( defval = 'REFERENCE' ),
    LREL = dict( datatype = int, defval = 0 ),
    NVER = dict( datatype = int, defval = 1 ),
    MATNUM = dict( datatype = 'matnumbers', defval = {} ),
    ENDATE = dict ( defval = '' ),
    EDATE = dict( datatype = 'datestr', defval = 'MMMYY' ),
    DDATE = dict( datatype = 'datestr', defval = 'MMMYY' ),
    RDATE = dict( datatype = 'datestr',
                  defval = 'MMMYY' ),
)

def _impl_get_metadata_params_and_docs():
    d = {}
    from .ncmat2endf import EndfMetaData
    for k, v in _metadata_definitions.items():
        doc = getattr(EndfMetaData,k.lower()).__doc__
        assert doc is not None
        d[k] = ' '.join(doc.strip().split())
    return d

def _impl_emd_set( now_MMMYY, data, param, value,  ):
    from .exceptions import NCBadInput
    k, v = param.upper(), value
    md = _metadata_definitions.get(k)
    if not md:
        raise NCBadInput(f'Invalid EndfMetaData parameter "{k}"')
    if v is None:
        v = md['defval']
        assert v is not None
        data[k] = v
        return

    if isinstance(v,str):
        for e in ['"',"'",'`']:
            if e in v:
                raise NCBadInput(f'Forbidden character {e} in value '
                                 f'of EndfMetaData parameter "{k}"')

    datatype = md.get('datatype',str)

    if isinstance(datatype,str) and datatype == 'datestr':
        if not isinstance( v, str ):
            raise NCBadInput('ENDF date value must be a string')
        if v.lower()=='now':
            v = now_MMMYY
        if len('MMMYY') != len(v):
            raise NCBadInput('ENDF date value is not in expected'
                             ' format, which is either special value'
                             ' "NOW" or a date in the format "MMMYY"'
                             ' (e.g. "Jun25").')
        data[k] = v
        return

    if isinstance(datatype,str) and datatype == 'matnumbers':
        if not hasattr(v,'items') and isinstance(v,str):
            d = {}
            for e in v.split(','):
                p = [_.strip() for _ in e.split(':')]
                if ( not len(p)==2 or not p[0] or not p[1]
                     or not p[1].isdigit() ):
                    raise NCBadInput(f'Invalid MATNUM value {repr(v)}:'
                                     ' must be a dict or string with a'
                                     ' format like "Zn:101,O:102".')
                d[p[0]]=int(p[1])
            v = d
        if not hasattr(v,'items'):
            raise NCBadInput('MATNUM must be a dict')
        for kk,vv in v.items():
            if type(kk) is not str or type(vv) is not int:
                raise NCBadInput('MATNUM must be a dict from element',
                                 ' labels (str) to material values (int)' )
        data[k] = v
        return

    if isinstance(v,str) and datatype is int and v.isdigit():
        v = int(v)

    if not isinstance( v, datatype ):
        raise NCBadInput(f'EndfMetaData parameter "{k}" data '
                         f'must be of type {datatype.__name__}')
    data[k] = v

def _dump_dict( d, prefix, lvl = 1, snip_output=True ):
    if not hasattr(d,'items'):
        s = repr(d)
        if len(s) > 80:
            s = s[0:35]+'<<SNIP>>'+s[-35:]
        ncprint(f'{prefix}{s}')
        return

    ld = list(d.items())
    keys_all_digits = all( str(k).isdigit() for k,v in ld)
    nlim = 30 if lvl<2 else 15
    if not keys_all_digits:
        nlim = 100
    for i,(k,v) in enumerate(ld):
        if len(ld)>nlim and nlim//3<i<len(ld)-nlim//3 and snip_output:
            if i== (nlim//3) + 1:
                ncprint(f'{prefix}<..SNIP..>')
            continue
        vs = repr(v)
        if len(vs) < 80:
            ncprint(f'{prefix}{repr(k)} -> {vs}')
        else:
            ncprint(f'{prefix}{repr(k)} ->')
            snip_output = ( False if k == 'DESCRIPTION' else snip_output )
            _dump_dict(v,prefix+'    ',lvl=lvl+1, snip_output=snip_output)

def _interp2d(x, y, x0, y0, z0=None):
    #
    # Bilinear interpolation on irregular cartesian grids
    # Interpolates a 2D array z0, with coordinates given by 1D arrays x0, y0
    # into a grid defined by 1D arrays x and y
    #
    if isinstance(x, (int, float)):
        x =_np.array([x])
    if isinstance(y, (int, float)):
        y =_np.array([y])
    if isinstance(y0, list):
        y =_np.array(y)
    if isinstance(x, list):
        x =_np.array(x)
    if isinstance(y0, list):
        y =_np.array(y)
    if isinstance(x0, list):
        x0 =_np.array(x0)
    if isinstance(y0, list):
        y0 =_np.array(y0)

    Nx = x0.size
    Ny = y0.size
    assert isinstance(z0,_np.ndarray)
    assert _np.shape(z0) == (Nx, Ny)
    assert Nx>1 and Ny >1

    # Create matrix of x and y coordinates
    xx, yy =_np.meshgrid(x, y, indexing='ij')

    # find neighbour points
    i2 =_np.searchsorted(x0, x)
    outside =_np.where(i2 > Nx - 1)
    i2[outside] = Nx - 2
    i1 = i2 - 1
    i1[outside] = Nx - 2
    i1[_np.where(i1 < 0)] = 0

    j2 =_np.searchsorted(y0, y)
    outside =_np.where(j2 > Ny - 1)
    j2[outside] = Ny - 2
    j1 = j2 - 1
    j1[outside] = Ny - 2
    j1[_np.where(j1 < 0)] = 0

    # get corner values
    ii1, jj1 =_np.meshgrid(i1, j1, indexing='ij')
    xx1, yy1 =_np.meshgrid(x0[i1], y0[j1], indexing='ij')
    ii2, jj2 =_np.meshgrid(i2, j2, indexing='ij')
    xx2, yy2 =_np.meshgrid(x0[i2], y0[j2], indexing='ij')

    ii11, jj11 =_np.meshgrid(i1, j1, indexing='ij')
    ii12, jj12 =_np.meshgrid(i1, j2, indexing='ij')
    ii21, jj21 =_np.meshgrid(i2, j1, indexing='ij')
    ii22, jj22 =_np.meshgrid(i2, j2, indexing='ij')

    z11 = z0[ii11, jj11]
    z12 = z0[ii12, jj12]
    z21 = z0[ii21, jj21]
    z22 = z0[ii22, jj22]

    # get deltas
    dxx, dyy =_np.meshgrid(_np.diff(x0)[i1],_np.diff(y0)[j1], indexing='ij')
    assert _np.all(dxx>0)
    assert _np.all(dyy>0)

    z = 1.0/(dxx*dyy)*(z11*(yy2 - yy)*(xx2 - xx) +
                       z21*(yy2 - yy)*(xx - xx1) +
                       z12*(yy - yy1)*(xx2 - xx) +
                       z22*(yy - yy1)*(xx - xx1))
    return z

def _get_scat_proc_names( scat_obj ):
    d = scat_obj.getSummary()
    if d['name'] != 'ProcComposition':
        return [ d['name'] ]
    return sorted( dd['name']
                   for frac,dd in d['specific']['components'] )

def _tidy_beta( x, allow_negative=False):
    if not unit_test_chop_vals[0]:
        return x
    if allow_negative:
        assert -1e99 <= x <= 1e99
    else:
        assert 0.0 <= x <= 1e99
    if x > 8:
        return float('%.2g'%x)
    return float('%.3g'%x)

def _tidy_alpha_list( a_values ):
    if not unit_test_chop_vals[0]:
        return a_values
    def _chop(x):
        assert 0.0 <= x <= 1e99
        return float('%.1g'%x)
    return  [ _chop(x) for x in a_values ]

def _tidy_teffwp( x ):
    if not unit_test_chop_vals[0]:
        return x
    assert 0.0 < x <= 1e99
    return float('%.13g'%x)

def _tidy_sab_list( s_values ):
    s_values = [ float(e) for e in s_values ]
    if not unit_test_chop_vals[0]:
        return s_values
    def _chop(x):
        assert 0.0 <= x <= 1e99
        #Have to be rather harsh unfortunately:
        if x < 1e-2:
            return 0.0
        return float('%.1g'%x)
    return  [ _chop(x) for x in s_values ]

def _detect_and_short_horisontal_ruler_in_line( line, hr_chars = '+-=~*^#' ):
    # Not the greatest implementation:
    is_hr = False
    _hrtries = 1000
    while _hrtries > 0 and any( 16*c in line for c in hr_chars ):
        is_hr = True
        while _hrtries > 0 and len(line) > ENDF_DESCR_MAXW:
            _hrtries -= 1
            _norig = len(line)
            for c in hr_chars:
                line = line.replace(16*c,8*c,1)
                if len(line)<ENDF_DESCR_MAXW:
                    _hrtries = 0
                    break
            if _norig == len(line):
                _hrtries = 0
                break#did not help, give up
    return is_hr, line

unit_test_abort_write = [False]
unit_test_chop_vals = [False]
unit_test_not_write_version = [False]
