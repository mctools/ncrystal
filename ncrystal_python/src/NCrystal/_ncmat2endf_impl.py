
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
from . import misc as nc_misc
from . import cfgstr as nc_cfgstr
from ._common import print as ncprint
from ._common import warn as ncwarn
from ._common import write_text as ncwrite_text
print = ncprint

mass_neutron = (nc_constants.const_neutron_mass_amu*
               nc_constants.constant_dalton2eVc2/
               ((nc_constants.constant_c*1e-12)**2)) # eV*ps^2*Angstrom^-2

hbar = nc_constants.constant_planck/nc_constants.k2Pi*1e12 # eV*ps
T0 = 293.6 # K - Reference temperature for LAT=1 in ENDF-6 MF=7/MT=4

_cacheimport=[None]
def import_endfparserpy():
    if _cacheimport[0] is not None:
        return _cacheimport[0]
    try:
        # TODO: temporary fix to avoid syntax warning from endf-parserpy
        # https://github.com/IAEA-NDS/endf-parserpy/issues/10
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=SyntaxWarning)
            import endf_parserpy
            from endf_parserpy.interpreter.fortran_utils import read_fort_floats
            from endf_parserpy.interpreter.fortran_utils import write_fort_floats
    except ImportError:
        raise SystemExit('Could not import endf_parserpy. Check the package '
                         'was correctly installed. This can be done with:'
                         'pip install ncrystal[endf], '
                         'pip install ncrystal[all], '
                         'pip install endf_parserpy if NCrystal was installed'
                         'with pip, or'
                         'conda install -c conda-forge endf_parserpy if'
                         'NCrystal was installed with conda.')
    _cacheimport[0] = ( endf_parserpy,read_fort_floats,write_fort_floats)
    return _cacheimport[0]

def _endf_roundoff(x):
    """Limit the precision of a float to what can be represented
       in an ENDF-6 file

    Parameters
    ----------
    x : Iterable of float
        Array to process

    Returns
    -------
    numpy array
        Processed array

    """
    _,read_fort_floats,write_fort_floats = import_endfparserpy()
    return _np.array(read_fort_floats(write_fort_floats(x, {'width':11}),
                                     n=len(x),read_opts={'width':11}))

def _endf_clean(x):
    """Return an array of unique floats that can be represented
       in an ENDF-6 file

    Parameters
    ----------
    x : Iterable of float
        Array to process

    Returns
    -------
    numpy array
        Processed array

    """
    return _np.unique(_endf_roundoff(x))

class ElementData():
    r"""Container for nuclear data for a single element or isotope.

    Attributes
    ----------
    alpha : numpy array
        alpha grid
    beta : numpy array
        positive beta grid
    beta_total : numpy array
        asymmetric beta grid
    sab : list of numpy array
        symmetric S(alpha, beta) table
    dwi : numpy array
        Debye-Waller integral
    teff : numpy array
        Effective temperatures for short collision time approximation
    elastic : string
        Elastic approximation used in the element
        (coherent, incoherent or mixed)
    awr : float
        Atomic mass in neutron mass units
    za : int
        ZAID (Z*1000 + A)
    zsymam : float
        Text representation of the element or isotope
    sigma_i : float
        Incoherent bound atom cross section
    sigma_free : float
        Scattering free atom cross section
    """

    def __init__(self, ad):
        r"""
        Parameters
        ----------
        ad : NCrystal AtomData
        """
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
    @sab_total.setter
    def sab_total(self, x):
        self._sab_total = x

    @property
    def dwi(self):
        return self._dwi
    @dwi.setter
    def dwi(self, x):
        self._dwi = x

    @property
    def teff(self):
        return self._teff
    @teff.setter
    def teff(self, x):
        self._teff = x

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
    r"""Container for nuclear data for a material.

    Attributes
    ----------
    comments : list of string
        Comments about the origin of the data
    temperatures : list or tuple of float
        List of temperatures to process
    ncmat_cfg : string
        NCrystal cfg string to convert
    composition : list of (float, NCrystal AtomData)
        Composition of the material
    elements : dictionary
        Nuclear data for each of the elements of isotopes in the materal
    edges : list numpy array
        Energies for the Bragg edges for each temperature
    sigmaE : list numpy array
        XS*E for the Bragg edges for each temperature
    elastic_mode: string
        Elastic approximation used in the material
        (greater, scaled or mixed)
    """

    def __init__(self, *, ncmat_cfg, temperatures, elastic_mode,
                          requested_emax, verbosity=1):
        r"""
        Parameters
        ----------
        ncmat_cfg : string
            NCrystal cfg string to convert
        temperatures : list or tuple of float
            List of temperatures to process
        elastic_mode : string
            Elastic approximation used in the material
            (greater, scaled or mixed)
        verbosity : integer
            Level of verbosity for the output
        """

        self._temperatures = tuple(temperatures)
        self._ncmat_cfg = ncmat_cfg
        info_obj = nc_core.createInfo(ncmat_cfg)
        self._composition = info_obj.composition
        self._elems = {}
        self._comments = None
        self._vdoslux = nc_cfgstr.decodecfg_vdoslux(ncmat_cfg)
        self._requested_emax = requested_emax
        self._verbosity = verbosity
        self._elastic_mode = elastic_mode
        scattering_components = nc_misc.detect_scattering_components(ncmat_cfg)
        self._enable_coh_elas = ( 'coh_elas' in scattering_components and
                                  info_obj.hasAtomInfo() )
        if not self._enable_coh_elas:
            ncwarn('Coherent elastic component disabled')
        self._enable_incoh_elas = ( 'incoh_elas' in scattering_components )
        if not self._enable_incoh_elas:
            ncwarn('Incoherent elastic component disabled')
        self._enable_inelas = ( 'inelas' in scattering_components )
        if not self._enable_inelas:
            ncwarn('Inelastic component disabled')
        # _combine_temperatures:
        # False: use (alpha, beta) grid for lowest temperature
        # True: combine all temperatures
        self._combine_temperatures = False
        for frac, ad in self._composition:
            if not ad.isNaturalElement():
                # TODO: properly handle isolated isotopes and enriched
                #       elements
                raise NotImplementedError('Conversion supported only for'
                                          ' natural elements')
            sym = ad.elementName()
            self._elems[sym] = ElementData(ad)
        if self._enable_coh_elas:
            self._edges = []
            self._sigmaE = []
        else:
            self._edges = None
            self._sigmaE = None
        self._incoherent_fraction = -1
        if (len(self._composition) > 1) and (elastic_mode == 'scaled'):
            #
            # Find element with minimum incoherent contribution.
            #
            self._designated_coherent_atom = None
            for frac, ad in self._composition:
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
        self._get_ncrystal_comments()

    @property
    def comments(self):
        return self._comments

    @property
    def temperatures(self):
        return self._temperatures

    @property
    def ncmat_cfg(self):
        return self._ncmat_cfg

    @property
    def composition(self):
        return self._composition

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

    def _loadKernel( self, di ):
        from .vdos import extractKnl
        #Note: using the extractKnl function rather than di.loadKernel means
        #that there will be no reduction of vdoslux for VDOSDebye objects. This
        #is why we use this function here, since ENDF files need the higher
        #energy range of the resulting sab.
        return extractKnl( vdos = di,
                           mass_amu = di.atomData.averageMassAMU(),
                           temperature = di.temperature,
                           scatxs = 1.0,#fixme: double check if 1.0 is appropriate here
                           target_emax = self._requested_emax,
                           vdoslux = self._vdoslux )

    def _combine_alpha_beta_grids(self):
        #
        # Combine (alpha, beta) grids from different temperatures
        # into a single grid. This usually results in a huge grid and
        # it is only kept as an option to debug libraries.
        #
        ncwarn('Combining (alpha, beta) grids from different temperatures.'
               ' This usually results in a huge grid.')
        for T in self._temperatures[1:]:
            cfg = self._ncmat_cfg+f';temp={T}K'
            info_obj = nc_core.createInfo(cfg)
            for di in info_obj.dyninfos:
                sym = di.atomData.elementName()
                sctknl = self._loadKernel(di)
                self._elems[sym].alpha = _np.unique(_np.concatenate((
                                         self._elems[sym].alpha,
                                         sctknl['alpha']*T/T0)))
                self._elems[sym].beta_total = _np.unique(_np.concatenate((
                                         self._elems[sym].beta_total,
                                         sctknl['beta']*T/T0)))

    def _get_alpha_beta_grid(self):
        T = self._temperatures[0]
        cfg = self._ncmat_cfg+f';temp={T}K'
        info_obj = nc_core.createInfo(cfg)
        for di in info_obj.dyninfos:
            sym = di.atomData.elementName()
            sctknl = self._loadKernel(di)
            self._elems[sym].alpha = sctknl['alpha']*T/T0
            self._elems[sym].beta_total = sctknl['beta']*T/T0
        if self._combine_temperatures:
            self._combine_alpha_beta_grids()

        for frac, ad in self._composition:
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
                ncprint(f'>>> alpha points: {len(self._elems[sym].alpha)}, '
                         'alpha range: '
                        f'({_np.min(self._elems[sym].alpha*T0/T)}'
                        f', {_np.max(self._elems[sym].alpha*T0/T)})')
                ncprint(f'>>> beta points: {len(self._elems[sym].beta)}, '
                        f'beta range: ({_np.min(self._elems[sym].beta*T0/T)},'
                        f' {_np.max(self._elems[sym].beta*T0/T)})')

    def _get_coherent_elastic(self, T):
        #
        # Load coherent elastic data
        #
        cfg = (self._ncmat_cfg+f';temp={T}K;comp=bragg')
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
                cfg = self._ncmat_cfg+f';temp={T}K'
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
                    self._elems[sym].dwi.append(msd*2*mass_neutron/hbar**2)
            else:
                for frac, ad in self._composition:
                    sym = ad.elementName()
                    self._elems[sym].sigma_i =  None
                    self._elems[sym].dwi =  None
        if self._verbosity > 1:
            ncprint('>> Prepare elastic approximations')
        if ( self._elastic_mode == 'scaled' ):
            if ( len(self._composition) > 1 and
                 self._incoherent_fraction < 1e-6 ):
                    self._elastic_mode = 'greater'
                    ncwarn('Scaled elastic mode requested '
                           'but all elements are coherent. '
                           '"greater" option will be used instead.')
        for frac, ad in self._composition:
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
            if self._elastic_mode == 'greater': # iel = 98
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
                if len(self._composition) == 1: # iel = 99, single atomic case
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

    def _interpolate_sab(self, a, b, s, sym, T):
        #
        # Interpolate S(a,b) because the NCrystal grid might
        # contain numbers that cannot be represented
        # as distinct FORTRAN reals in the ENDF-6 file
        #
        # FIXME: replace scipy interpolation with other implementation
        import scipy.interpolate as scint
        s.shape = (len(b), len(a))
        aint, bint = _np.meshgrid(self._elems[sym].alpha*T0/T,
                                  self._elems[sym].beta_total*T0/T)
        sab_int = scint.interpn((a, b),
                                s.transpose(),
                                _np.column_stack((aint.ravel(), bint.ravel())),
                                bounds_error=False, fill_value=0.0,
                                method='linear')
        sab_int.shape = _np.shape(aint)
        return sab_int

    def _get_inelastic_data(self):
        for T in self._temperatures:
            cfg = self._ncmat_cfg+f';temp={T}K'
            info_obj = nc_core.createInfo(cfg)
            for di in info_obj.dyninfos:
                sym = di.atomData.elementName()
                #
                # Load incoherent inelastic data
                #
                sctknl = self._loadKernel(di)
                if self._verbosity > 2:
                    ncprint(f'>>> Interpolating T={T}K for {sym}')

                alpha = sctknl['alpha']
                beta = sctknl['beta']
                sab = sctknl['sab']
                sab_int = self._interpolate_sab(alpha, beta, sab, sym, T )
                self._elems[sym].sab_total.append(sab_int)
                emin = di.vdosData()[0][0]
                emax = di.vdosData()[0][1]
                rho = di.vdosData()[1]
                res = nc_vdos.analyseVDOS(emin, emax, rho, di.temperature,
                                          di.atomData.averageMassAMU())
                self._elems[sym].teff.append(res['teff'])

    def _get_ncrystal_comments(self):
        # TODO: handle multi phase materials
        ncmat_fn = nc_cfgstr.decodeCfg(self._ncmat_cfg)['data_name']
        td = nc_core.createTextData(ncmat_fn)
        from ._ncmatimpl import _extractInitialHeaderCommentsFromNCMATData
        from textwrap import wrap
        comments = _extractInitialHeaderCommentsFromNCMATData(td)
        self._comments = []
        for paragraph in '\n'.join(comments).split('\n'):
            if paragraph == '':
                self._comments += ['']
            else:
                self._comments += wrap(paragraph.lstrip(),
                                                width=66,
                                                break_long_words=False)

class EndfFile():
    r"""Creates thermal ENDF file.
       using endf-parserpy

    Methods
    ----------
    write(endf_fn)
        Write ENDF file.
    """
    def __init__(self, element, data, mat, endf_metadata, *,
                 include_gif=False, isotopic_expansion=False,
                 smin=None, emax=None, lasym=None, verbosity=1):
        r"""
        Parameters
        ----------
        element : string
            Element to be output

        data : NuclearData
            Nuclear data for the material

        mat: int
            ENDF material number

        endf_metadata : EndfMetaData
            Metadata for the ENDF-6 file

        include_gif: boolean
            Include the generalized information in MF=7/MT=451 in isotopes

        isotopic_expansion: boolean
            Expand the information in MF=7/MT=451 in isotopes

        verbosity : int
            Level of verbosity of the output (0: quiet)
        """
        endf_parserpy,_,_ = import_endfparserpy()

        self._endf_parserpy_version = endf_parserpy.__version__
        self._sym = element
        self._data = data
        self._endf_metadata = endf_metadata
        self._mat = mat
        self._include_gif = include_gif
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
        """
        Creates endf-parserpy dictionary for MF=7 file,
        MT=2 reaction (thermal elastic) of a thermal ENDF file
        """

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
                    S[q][i] = data.sigmaE[i][q-1]
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
        """
        Creates endf-parserpy dictionary for MF=7 file,
        MT=4 reaction (thermal inelastic) of a thermal ENDF file
        """

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
            d['beta'] = {k:v for k, v in enumerate(beta,start=1)}
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
                S1[j]['alpha'] = (alpha/awr).tolist()
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
                        S2[q][j] = {k: v for k,v in enumerate(sab, start =1)}
            d['S'] = S2
        elif (self._lasym == 1) or (self._lasym == 3):
            # Save S(a,b) or S(a,b)*exp(-b/2) for all beta
            alpha = data.elements[self._sym].alpha
            beta = data.elements[self._sym].beta_total

            d['NB'] = len(beta)
            d['beta_interp/NBT'] = [len(beta)]
            d['beta_interp/INT'] = [4]

            d['T0'] = temperatures[0]
            d['beta'] = {k:v for k, v in enumerate(beta,start=1)}
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
                S1[j]['alpha'] = (alpha/awr).tolist()
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
                        S2[q][j] = {k: v for k,v in enumerate(sab, start =1)}
            d['S'] = S2

        d['teff0_table/Teff0'] = data.elements[self._sym].teff
        d['teff0_table/Tint'] = list(temperatures)
        d['teff0_table/NBT'] = [len(temperatures)]
        d['teff0_table/INT'] = [2]

    def _createMF7(self):
        """Creates MF=7 file of a thermal ENDF file.
           See ENDF-102, sect. 7.
           https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf
        """
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
            if self._isotopic_expansion:
                # TODO: implement isotopic expansion
                pass
            else:
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
        desc.append(66*'*')
        desc.append('')
        desc.append(' This file was converted from the following NCrystal cfg')
        desc.append(' string [1]:')
        desc.append('')

        if len(data.ncmat_cfg) < 60:
            desc.append(f'"{data.ncmat_cfg}"'.center(66))
        else:
            cfg_parts = data.ncmat_cfg.split(';')
            cfg_lines = []
            from textwrap import wrap
            for part in cfg_parts:
                cfg_lines += wrap(part.lstrip(),
                                  width=60,
                                  break_long_words=True)
                cfg_lines[-1] += ';'
            if len(cfg_lines) == 1:
                desc.append(f'  "{cfg_lines[0]}"')
            elif len(cfg_lines) == 2:
                desc.append(f'  "{cfg_lines[0]}')
                desc.append(f'   {cfg_lines[-1]}"')
            elif len(cfg_lines) > 1:
                desc.append(f'  "{cfg_lines[0]}')
                desc.extend(f'   {_}' for _ in cfg_lines[1:-1])
                desc.append(f'   {cfg_lines[-1]}"')

        desc.append('')
        nc_version = nc_core.get_version()
        ep_version = self._endf_parserpy_version
        if unit_test_chop_svals[0]:
            nc_version = 'NCVERSION'
            ep_version = 'EPVERSION'
        desc.append(f' using NCrystal {nc_version} ')
        desc.append(f' and endf-parserpy [2] {ep_version} ')
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
        desc.append('[2] https://endf-parserpy.readthedocs.io/en/latest/')
        desc.append('')
        desc.append(66*'*')
        desc.append('')
        desc.append('Comments from NCMAT file:')
        desc.append('')
        for line in data.comments:
            desc.append(line)
        desc.append(66*'*')
        for _ in desc:
            assert len(_) <= 66
        return [_.ljust(66) for _ in desc]

    def _createMF1(self):
        """Creates MF=1 file of a thermal ENDF file.
           See ENDF-102, sect. 1.
           https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf
        """

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
        d['HSUB/1'].ljust(66)
        d['HSUB/2'] =  '-----THERMAL NEUTRON SCATTERING DATA'.ljust(66)
        d['HSUB/3'] =  '------ENDF-6 FORMAT'.ljust(66)
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

    def write(self, endf_fn, force_save):
        import pathlib
        if self._verbosity > 0:
            ncprint(f'Write ENDF file {endf_fn}...')
        outfile = pathlib.Path(endf_fn)
        if outfile.exists() and not force_save:
            raise SystemExit('Error: output file already exists'
                             ' (run with --force to overwrite)')
        if not outfile.parent.is_dir():
            raise SystemExit('Error: output directory does not exist:'
                             f' { outfile.parent }')

        if is_unit_test[0]:
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

def _dump_dict( d, prefix, lvl = 1 ):
    if not hasattr(d,'items'):
        s = repr(d)
        if len(s) > 80:
            s = s[0:35]+'<<SNIP>>'+s[-35:]
        ncprint(f'{prefix}{s}')
        return

    ld = list(d.items())#fixme: sorted(..) ?
    keys_all_digits = all( str(k).isdigit() for k,v in ld)
    nlim = 30 if lvl<2 else 15
    if not keys_all_digits:
        nlim = 100
    for i,(k,v) in enumerate(ld):
        if len(ld)>nlim and nlim//3<i<len(ld)-nlim//3:
            if i== (nlim//3) + 1:
                ncprint(f'{prefix}<..SNIP..>')
            continue
        vs = repr(v)
        if len(vs) < 80:
            ncprint(f'{prefix}{repr(k)} -> {vs}')
        else:
            ncprint(f'{prefix}{repr(k)} ->')
            _dump_dict(v,prefix+'    ',lvl=lvl+1)

def _tidy_sab_list( s_values ):
    if not unit_test_chop_svals[0]:
        return s_values
    _npf = _np.float64
    def _chop(x):
        assert 0.0 <= x <= 1e99
        if x > 1e-20:
            return x
        if x > 1e-20:
            return float('%.4g'%x)
        if x > 1e-60:
            return float('%.3g'%x)
        if x > 1e-80:
            return float('%.2g'%x)
        return float('%.1g'%x)
    return  [ _chop(x) for x in s_values ]

is_unit_test = [False]


unit_test_chop_svals = [False]
