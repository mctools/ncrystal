
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

Utility for creating a set of ENDF-6 thermal scattering files from a .ncmat
file. Parameters for the ENDF-6 file can be defined as arguments. It  allows to
handle multiple temperatures in one ENDF-6 file, but this is not
recommended, because NCrystal computes an optimal (alpha, beta) grid for each
material and temperature.

This utility uses the endf-parserpy package from IAEA to format and check the
syntaxis of the ENDF-6 file.

"""

__all__ = [ 'ncmat2endf']

import numpy as np
import scipy.interpolate as scint
from datetime import datetime
import warnings
from . import core as nc_core
from . import constants as nc_constants
from . import vdos as nc_vdos
from ._common import print

try:
    # TODO: temporary fix to avoid syntax warning from endf-parserpy
    # https://github.com/IAEA-NDS/endf-parserpy/issues/10
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=SyntaxWarning)
        import endf_parserpy
    from endf_parserpy.interpreter.fortran_utils import read_fort_floats
    from endf_parserpy.interpreter.fortran_utils import write_fort_floats
except ImportError:
    raise SystemExit('Could not import endf_parserpy. Check the package was'+
                     'correctly installed, with a version equal or higher'+
                     'than 0.11.0.\nhttps://endf-parserpy.readthedocs.io/')
version_num = sum([int(x)*10**(3*n)
                   for x,n in zip(endf_parserpy.__version__.split('.'),
                   reversed(range(3)))])
assert version_num >=  1100, 'Too old endf-parserpy found. '+\
                             'Version 0.11.0 or above required'


available_elastic_modes = ('greater', 'scaled', 'mixed')
mass_neutron = (nc_constants.const_neutron_mass_amu*
               nc_constants.constant_dalton2eVc2/
               ((nc_constants.constant_c*1e-12)**2)) # eV*ps^2*Angstrom^-2
hbar = nc_constants.constant_planck/nc_constants.k2Pi*1e12 # eV*ps
T0 = 293.6 # K
kT0 = nc_constants.constant_boltzmann*T0 # eV

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
    return np.array(read_fort_floats(write_fort_floats(x, {'width':11}),
                                     n=len(x),read_opts={'width':11}))

def _wrap_string(inp, lim=66):
    """Return string with lines wrapped to a given width.

    Parameters
    ----------
    inp : str
        String to wrap

    lim : int
        Maximum line width

    Returns
    -------
    s: str
        Output string with lines wrapped

    """
    word_list = []
    for line in inp.split("\n"):
        if line == "":
            word_list.append('')
            continue

        total_width = 0
        current_line_words = []

        if len(line.split()) == 1:
            word = line.split()[0]
            num_segments = len(word) // lim
            for i in range(num_segments):
                word_list.append(word[i*lim:(i+1)*lim])
            word_list.append(word[num_segments*lim:])
        else:
            for word in line.split():
                if total_width + len(word) + 1 <= lim:
                    current_line_words.append(word)
                    total_width += len(word) + 1
                else:
                    word_list.append(" ".join(current_line_words))
                    current_line_words = [word]
                    total_width = len(word)
            if current_line_words:
                word_list.append(" ".join(current_line_words))
    return('\n'.join(word_list))

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
    sab : iterable of numpy array
        symmetric S(alpha, beta) table
    dwi : iterable of float
        Debye-Waller integral
    teff : iterable of float
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
        self._alpha = np.array([])
        self._beta = np.array([])
        self._beta_total = np.array([])
        self._sab_total = []
        self._teff = []
        self._elastic = None
        self._sym = ad.displayLabel()
        Z = '{:3d}'.format(ad.Z())
        A = '   ' if ad.isNaturalElement() else '{:3d}'.format(ad.A())
        sym = self._get_symbol(self._sym).ljust(2)
        self._zsymam = '{:3s}-{:2s}-{:3s}'.format(Z,sym,A)
        self._za = ad.Z()*1000+ad.A()

    def _get_symbol(self, isotope_name):
        """Get element symbol from isotope name. E.g. Be9 -> Be

        Parameters
        ==========
        isotope_name: string

        Returns
        ==========
        element_symbol: string

        """

        symbol = ''
        for c in isotope_name:
            if ((ord(c) >= 97 and ord(c) <= 122) or
                (ord(c) >= 65 and ord(c) <= 90)):
                symbol += c
            else:
                break
        return symbol


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
    ncrystal_comments : string
        Comments in the ncmat file
    temperatures : iterable of flotat
        List of temperatures to process
    ncmat_fn : string
        NCrystal ncmat filename to convert
    composition : iterable of tuples (float, NCrystal AtomData)
        Composition of the material
    elements : iterable of ElementData
        Nuclear data for each of the elements of isotopes in the materal
    edges : iterable of numpy array
        Energies for the Bragg edges for each temperature
    sigmaE : iterable of numpy array
        XS*E for the Bragg edges for each temperature
    """

    def __init__(self, ncmat_fn, temperatures, elastic_mode,
                 vdoslux=3, verbosity=1):
        r"""
        Parameters
        ----------
        ncmat_fn : string
            NCrystal ncmat filename to convert
        temperatures : iterable of float
            List of temperatures to process
        elastic_mode : string
            Elastic approximation used in the material
            (greater, scaled or mixed)
        vdoslux : integer
            Level of luxury to generate data in NCrystal
        verbosity : integer
            Level of verbosity for the output
        """
        self._temperatures = np.sort(np.asarray(temperatures))
        self._ncmat_fn = ncmat_fn
        mat = nc_core.load(ncmat_fn+f';vdoslux={vdoslux}')
        self._composition = mat.info.composition
        self._elems = {}
        self._ncrystal_comments = None
        self._vdoslux = vdoslux
        self._verbosity = verbosity
        # _combine_temperatures:
        # False: use (alpha, beta) grid for lowest temperature
        # True: combine all temperatures
        self._combine_temperatures = False
        if elastic_mode not in available_elastic_modes:
            raise ValueError(f'Elastic mode {elastic_mode}'+
                             f' not in {available_elastic_modes}')
        for frac, ad in self._composition:
            sym = ad.displayLabel()
            self._elems[sym] = ElementData(ad)
        if mat.info.hasAtomInfo():
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
                sym = ad.displayLabel()
                if (frac/(1.0-frac)*ad.incoherentXS() <
                    self._incoherent_fraction or
                    self._incoherent_fraction == -1):
                    self._incoherent_fraction = (frac/
                                                 (1.0-frac)*ad.incoherentXS())
                    self._designated_coherent_atom = sym
            if self._verbosity > 1:
                print(f'Designated incoherent: {sym}')
        self._get_alpha_beta_grid()
        self._get_elastic_data(elastic_mode)
        self._get_inelastic_data()
        self._get_ncrystal_comments()

    @property
    def ncrystal_comments(self):
        return self._ncrystal_comments

    @property
    def temperatures(self):
        return self._temperatures

    @property
    def ncmat_fn(self):
        return self._ncmat_fn

    @property
    def composition(self):
        return self._composition
    @composition.setter
    def composition(self, x):
        self._composition = x

    @property
    def elements(self):
        return self._elems
    @elements.setter
    def elements(self, x):
        self._elems = x

    @property
    def edges(self):
        return self._edges
    @edges.setter
    def edges(self, x):
        self._edges = x

    @property
    def sigmaE(self):
        return self._sigmaE
    @sigmaE.setter
    def sigmaE(self, x):
        self._sigmaE = x

    def _get_alpha_beta_grid(self):
        T = self._temperatures[0]
        cfg = f'{self._ncmat_fn};temp={T}K;vdoslux={self._vdoslux}'
        m = nc_core.load(cfg)
        kT = kT0*T/T0 # eV
        for di in m.info.dyninfos:
            sym = di.atomData.displayLabel()
            if type(di) in (nc_core.Info.DI_VDOS, nc_core.Info.DI_VDOSDebye):
                sctknl = di.loadKernel(vdoslux=self._vdoslux)
                self._elems[sym].alpha = sctknl['alpha']*kT/kT0
                self._elems[sym].beta_total = sctknl['beta']*kT/kT0
            else:
                raise NotImplementedError('Conversion supported only for VDOS'+
                                          ' and VDOSDebye dyninfos')
        if self._combine_temperatures:
            #
            # Combine (alpha, beta) grids from different temperatures
            # into a single grid. This usually results in a huge grid and
            # it is only kept as an option to debug libraries.
            #
            for T in self._temperatures[1:]:
                cfg = f'{self._ncmat_fn};temp={T}K;vdoslux={self._vdoslux}'
                m = nc_core.load(cfg)
                kT = kT0*T/T0 # eV
                for di in m.info.dyninfos:
                    sym = di.atomData.displayLabel()
                    sctknl = di.loadKernel(vdoslux=self._vdoslux)
                    _ =  np.unique(np.concatenate((self._elems[sym].alpha,
                                                   sctknl['alpha']*kT/kT0)))
                    self._elems[sym].alpha = _
                    _ = np.unique(np.concatenate((self._elems[sym].beta_total,
                                                  sctknl['beta']*kT/kT0)))
                    self._elems[sym].beta_total = _
        for frac, ad in self._composition:
            sym = ad.displayLabel()
            #
            # Remove points that cannot be represented in ENDF data
            #
            _ = np.unique(_endf_roundoff(self._elems[sym].alpha))
            self._elems[sym].alpha = _
            _ = np.unique(_endf_roundoff(self._elems[sym].beta_total))
            self._elems[sym].beta_total = _
            x = self._elems[sym].beta_total[
                np.where(self._elems[sym].beta_total<=0)] # get negative beta
            self._elems[sym].beta = -x[::-1] # Invert beta and change sign
            self._elems[sym].beta[0] = 0.0
            if self._verbosity > 2:
                print(f'>>> alpha points: {len(self._elems[sym].alpha)}, '+
                      f'alpha range: ({np.min(self._elems[sym].alpha*kT0/kT)}'+
                      f', {np.max(self._elems[sym].alpha*kT0/kT)})')
                print(f'>>> beta points: {len(self._elems[sym].beta)}, '+
                      f'beta range: ({np.min(self._elems[sym].beta*kT0/kT)},'+
                      f' {np.max(self._elems[sym].beta*kT0/kT)})')

    def _get_coherent_elastic(self, m, T):
        #
        # Load coherent elastic data
        #
        if T == self._temperatures[0]:
            # Find unique Bragg edges, as represented in ENDF-6 floats
            edges = np.array([nc_constants.wl2ekin(2.0*e.dspacing)
                              for e in m.info.hklObjects()])
            self._edges = np.unique(_endf_roundoff(edges))
            # Coherent scattering XS is evaluated between edges
            eps = 1e-3
            _ = np.concatenate((self._edges[:-1]**(1-eps)*
                              self._edges[1:]**eps, [self._edges[-1]]))
            self._evalpoints = _
        sigmaE = _endf_roundoff(m.scatter.xsect(self._evalpoints)
                                *self._evalpoints)
        assert np.all(sigmaE[:-1] <= sigmaE[1:]),\
               'Sigma*E in Bragg edges not cummulative'
        self._sigmaE.append(sigmaE)

    def _get_elastic_data(self, elastic_mode):
        for T in self._temperatures:
            cfg = f'{self._ncmat_fn};temp={T}K;vdoslux={self._vdoslux};'+\
                   'comp=bragg;dcutoff=0.1'
            m = nc_core.load(cfg)
            if m.info.hasAtomInfo():
                self._get_coherent_elastic(m, T)
            #
            for di in m.info.dyninfos:
                sym = di.atomData.displayLabel()
                emin = di.vdosData()[0][0]
                emax = di.vdosData()[0][1]
                rho = di.vdosData()[1]
                res = nc_vdos.analyseVDOS(emin, emax, rho, di.temperature,
                                          di.atomData.averageMassAMU())
                #
                # Load incoherent elastic data
                #
                if m.info.stateOfMatter().name == 'Solid':
                    msd = res['msd']
                    self._elems[sym].dwi.append(msd*2*mass_neutron/hbar**2)
        if self._verbosity > 1:
            print('>> Prepare elastic approximations')
        if elastic_mode == 'scaled' and self._incoherent_fraction < 1e-6:
            elastic_mode = 'greater'
            # TODO: replace this by a warning
            if self._verbosity>1:
                print('>> Scaled elastic mode requested'+
                      'but all elements are coherent.')
        for frac, ad in self._composition:
            sym = ad.displayLabel()
            if elastic_mode == 'mixed': # iel = 100
                if (self._sigmaE is None):
                    # mixed elastic requested but only incoherent available
                    # TODO: replace this by a warning
                    if self._verbosity>1:
                        print(f'>> Mixed elastic mode for {sym} but no '+
                               'Bragg edges found: incoherent approximation')
                    self._elems[sym].sigma_i = (ad.incoherentXS() +
                                                ad.coherentXS())
                    self._elems[sym].elastic = 'incoherent'
                else:
                    if self._verbosity>1:
                        print(f'>> Mixed elastic mode for {sym}')
                    self._elems[sym].elastic = 'mixed'
            if elastic_mode == 'greater': # iel = 98
                if ((ad.incoherentXS() > ad.coherentXS()) or
                    (self._sigmaE is None)):
                    if self._verbosity>1:
                        print( '>> Principal elastic mode '+
                              f'for {sym}: incoherent')
                    self._edges = None
                    self._sigmaE = None
                    self._elems[sym].elastic = 'incoherent'
                else:
                    if self._verbosity>1:
                        print(f'>> Principal elastic mode for {sym}: coherent')
                    self._elems[sym].sigma_i =  None
                    self._elems[sym].dwi =  None
                    self._elems[sym].elastic = 'coherent'
            elif elastic_mode == 'scaled':
                if len(self._composition) == 1: # iel = 99, single atomic case
                    if ((ad.incoherentXS() > ad.coherentXS()) or
                        (self._sigmaE is None)):
                        if self._verbosity>1:
                            print( '>> Scaled elastic mode for '+
                                  f'single atom {sym}: incoherent')
                        self._edges = None
                        self._sigmaE = None
                        self._elems[sym].sigma_i = (ad.incoherentXS() +
                                                    ad.coherentXS())
                        self._elems[sym].elastic = 'incoherent'
                    else:
                        if self._verbosity>1:
                            print( '>> Scaled elastic mode for '+
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
                        print(f'>> Scaled elastic mode for {sym}: '+
                               'incoherent approximation')
                    self._elems[sym].sigma_i = (ad.incoherentXS() +
                                                ad.coherentXS())
                    self._elems[sym].elastic = 'incoherent'
                else:
                    # iel = 99, multi atomic case
                    if sym == self._designated_coherent_atom:
                        if self._verbosity>1:
                            print(f'>> Scaled elastic mode for {sym} in '+
                                   'compound: designated coherent atom, '+
                                  f'dividing by frac^2={frac**2}')
                        self._elems[sym].sigma_i =  None
                        self._elems[sym].dwi =  None
                        self._elems[sym].elastic = 'coherent'
                        self._sigmaE = [x/frac for x in self._sigmaE]
                    else:
                        if self._verbosity>1:
                            print(f'>> Scaled elastic mode for {sym} '+
                                   'in compound: incoherent')
                        self._elems[sym].elastic = 'incoherent'
                        _ = ((1.0+self._incoherent_fraction/ad.incoherentXS())*
                             self._elems[sym].sigma_i)
                        self._elems[sym].sigma_i = _

    def _get_inelastic_data(self):
        for T in self._temperatures:
            m = nc_core.load(f'{self._ncmat_fn};temp={T}K')
            for di in m.info.dyninfos:
                sym = di.atomData.displayLabel()
                #
                # Load incoherent inelastic data
                #
                sctknl = di.loadKernel(vdoslux=self._vdoslux)
                if self._verbosity > 2:
                    print(f'>>> Interpolating T={T}K for {sym}')

                alpha = sctknl['alpha']
                beta = sctknl['beta']
                sab = sctknl['sab']
                sab.shape = (len(beta), len(alpha))
                kT = kT0/T0*T # eV
                _ = np.meshgrid(self._elems[sym]._alpha*kT0/kT,
                                self._elems[sym]._beta_total*kT0/kT)
                alpha_grid, beta_grid = _
                #
                # We need to interpolate S(a,b) because the NCrystal grid might
                # contain numbers that cannot be represented
                # as distinct FORTRAN reals in the ENDF-6 file
                #
                _ = np.column_stack((alpha_grid.ravel(), beta_grid.ravel()))
                sab_int = scint.interpn((alpha, beta), sab.transpose(), _,
                                        bounds_error=False, fill_value=0.0,
                                        method='linear')
                sab_int.shape = np.shape(alpha_grid)
                self._elems[sym]._sab_total.append(sab_int)
                emin = di.vdosData()[0][0]
                emax = di.vdosData()[0][1]
                rho = di.vdosData()[1]
                res = nc_vdos.analyseVDOS(emin, emax, rho, di.temperature,
                                          di.atomData.averageMassAMU())
                self._elems[sym]._teff.append(res['teff'])
    def _get_ncrystal_comments(self):
        _ = [line[1:] for line in
             nc_core.createTextData(self._ncmat_fn).rawData.split('\n')[:]
             if (len(line) > 0 and line[0] == '#')]
        self._ncrystal_comments = _wrap_string("\n".join(_),66)

class EndfFile():
    r"""Creates thermal ENDF file.
       using endf-parserpy

    Methods
    ----------
    write(endf_fn)
        Write ENDF file.
    """
    def __init__(self, element, data, mat, endf_parameters,
                 include_gif=False, isotopic_expansion=False,
                 parameter_description=None, verbosity=1):
        r"""
        Parameters
        ----------
        element : string
            Element to be output

        data : NuclearData
            Nuclear data for the material

        mat: int
            ENDF material number

        endf_parameters : EndfParameters
            Parameters for the ENDF-6 file

        include_gif: boolean
            Include the generalized information in MF=7/MT=451 in isotopes

        isotopic_expansion: boolean
            Expand the information in MF=7/MT=451 in isotopes

        parameter_description: iterable of string
            List of parameters used to generate the file

        verbosity : int
            Level of verbosity of the output (0: quiet)
        """
        self._endf_dict = endf_parserpy.EndfDict()
        self._parser = endf_parserpy.EndfParser(explain_missing_variable=True,
                                                cache_dir=False)
        self._sym = element
        self._mat = mat
        assert ((not isotopic_expansion) or include_gif),\
                'Isotopic expansion requires generalized information file'+\
                ', use --gif'
        self._include_gif = include_gif
        self._isotopic_expansion = isotopic_expansion
        self._parameter_description = parameter_description
        self._verbosity = verbosity
        self._endf_dict['0/0'] = {}
        self._endf_dict['0/0']['MAT'] = self._mat
        self._endf_dict['0/0']['TAPEDESCR'] = 'Created with ncmat2endf'
        self._createMF1(data, endf_parameters)
        self._createMF7(data, endf_parameters)
        endf_parserpy.update_directory(self._endf_dict, self._parser)

    def _createMF7(self, data, endf_parameters):
        """Creates MF=7 file of a thermal ENDF file.
           See ENDF-102, sect. 7.

        Parameters
        ----------
        data : NuclearData
            Nuclear data for the material

        endf_parameters : EndfParameters
            Parameters for the ENDF-6 file

        verbosity : int
            Level of verbosity of the output (0: quiet)
        """
        if self._verbosity > 1:
            print('> Generate MF7')
        awr = data.elements[self._sym].awr
        mat = self._mat
        za = data.elements[self._sym].za
        temperatures = data.temperatures
        #
        # Prepare dictionary for thermal elastic reaction
        #
        elastic = data.elements[self._sym].elastic
        if elastic is not None:
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
                d['Tint'] = temperatures.tolist()
                d['INT'] = [2]
                d['NBT'] = [len(temperatures)]
            lthr_values = {'coherent':1, 'incoherent':2, 'mixed':3}
            d['LTHR'] = lthr_values[data.elements[self._sym].elastic]
        #
        # Prepare dictionary for thermal inelastic reaction
        #
        self._endf_dict['7/4'] = {}
        d = self._endf_dict['7/4']
        d['MAT'] = mat
        d['ZA'] = za
        d['AWR'] = awr
        d['LAT'] =  1 # (alpha, beta) grid written for kT0 = 0.0253 eV
        d['LASYM'] = endf_parameters.lasym # symmetric/asymmetric S(a,b)
        d['LLN'] = 0 # linear S is stored
        d['NI'] = 6
        d['NS'] = 0
        d['B'] = {1:data.elements[self._sym].sigma_free,
                  2:endf_parameters.emax,
                  3:awr,
                  4:endf_parameters.emax,
                  5:0,                                # unused
                  6:1                                 # natom
                 }
        alpha = data.elements[self._sym].alpha
        beta = data.elements[self._sym].beta
        sab_data = []
        for sab_total, T in zip(data.elements[self._sym].sab_total,
                                temperatures):
            kT = T/T0*kT0
            _ = np.meshgrid(alpha*kT0/kT,
                            data.elements[self._sym].beta_total*kT0/kT)
            alpha_grid, beta_grid = _
            if (endf_parameters.lasym == 0) or (endf_parameters.lasym == 1):
                detailed_balance_factor = np.exp(beta_grid/2)
            if endf_parameters.lasym == 3:
                # S(a,b) for all beta
                sab_data.append(sab_total.transpose())
                continue
            if endf_parameters.lasym == 2:
                # S(a,b) for negative beta
                # get negative branch of S(a,b)
                sab_sym2 = sab_total[np.where(beta_grid<=0)]
                sab_sym2.shape = (len(beta), len(alpha))
                sab_sym3 = sab_sym2[::-1,:]  # Invert S(a,b) for negative beta
                sab_data.append(sab_sym3.transpose())
                continue
            if endf_parameters.lasym == 1:
                # S(a,b)*exp(-b/2) for all beta
                sab_sym = sab_total*detailed_balance_factor
                sab_data.append(sab_sym.transpose())
                continue
            if endf_parameters.lasym == 0:
                # S(a,b)*exp(-b/2) for negative beta
                sab_sym = sab_total*detailed_balance_factor
                sab_sym2 = sab_sym[np.where(beta_grid<=0)] # get negative branch of S(a,b)
                sab_sym2.shape = (len(beta), len(alpha))
                sab_sym3 = sab_sym2[::-1,:]  # Invert S(a,b) for negative beta
                sab_data.append(sab_sym3.transpose())
                continue

        if (endf_parameters.lasym == 0) or (endf_parameters.lasym == 2):
            # Save S(a,b) or S(a,b)*exp(-b/2) for negative beta
            d['NB'] = len(beta)
            d['beta_interp/NBT'] = [len(beta)]
            d['beta_interp/INT'] = [4]

            d['T0'] = temperatures[0]
            d['beta'] = {k:v for k, v in enumerate(beta,start=1)}
            d['LT'] = {k:len(temperatures)-1 for k, v in enumerate(beta,start=1)}
            d['T'] = {k:v for k,v in enumerate(temperatures[1:], start=1)}
            d['LI'] = {k:4 for k,v in enumerate(temperatures[1:], start=1)}
            d['NP'] = len(alpha)
            S1 = {}
            sab = sab_data[0]
            sab[sab < endf_parameters.smin] = 0.0
            for j,v in enumerate(beta, start=1):
                S1[j] = {}
                S1[j]['NBT'] = [len(alpha)]
                S1[j]['INT'] = [4]
                S1[j]['alpha'] = (alpha/awr).tolist()
                S1[j]['S'] = sab[:,j-1].tolist()
            d['S_table'] = S1

            S2 = {}
            if len(temperatures) > 1:
                for q, v in enumerate(alpha, start=1):
                    S2[q] = {}
                    for j,v in enumerate(beta, start=1):
                        sab = []
                        for i,v in enumerate(temperatures[1:], start=1):
                                sval = sab_data[i][q-1,j-1]
                                if sval < endf_parameters.smin:
                                    sval = 0.0
                                sab.append(sval)
                        S2[q][j] = {k: v for k,v in enumerate(sab, start =1)}
            d['S'] = S2
        elif (endf_parameters.lasym == 1) or (endf_parameters.lasym == 3):
            # Save S(a,b) or S(a,b)*exp(-b/2) for all beta
            alpha = data.elements[self._sym].alpha
            beta = data.elements[self._sym].beta_total

            d['NB'] = len(beta)
            d['beta_interp/NBT'] = [len(beta)]
            d['beta_interp/INT'] = [4]

            d['T0'] = temperatures[0]
            d['beta'] = {k:v for k, v in enumerate(beta,start=1)}
            d['LT'] = {k:len(temperatures)-1 for k, v in enumerate(beta,start=1)}
            d['T'] = {k:v for k,v in enumerate(temperatures[1:], start=1)}
            d['LI'] = {k:4 for k,v in enumerate(temperatures[1:], start=1)}
            d['NP'] = len(alpha)
            S1 = {}
            sab = sab_data[0]
            sab[sab < endf_parameters.smin] = 0.0
            for j,v in enumerate(beta, start=1):
                S1[j] = {}
                S1[j]['NBT'] = [len(alpha)]
                S1[j]['INT'] = [4]
                S1[j]['alpha'] = (alpha/awr).tolist()
                S1[j]['S'] = sab[:,j-1].tolist()
            d['S_table'] = S1

            S2 = {}
            if len(temperatures) > 1:
                for q, v in enumerate(alpha, start=1):
                    S2[q] = {}
                    for j,v in enumerate(beta, start=1):
                        sab = []
                        for i,v in enumerate(temperatures[1:], start=1):
                                sval = sab_data[i][q-1,j-1]
                                if sval < endf_parameters.smin:
                                    sval = 0.0
                                sab.append(sval)
                        S2[q][j] = {k: v for k,v in enumerate(sab, start =1)}
            d['S'] = S2

        d['teff0_table/Teff0'] = data.elements[self._sym].teff
        d['teff0_table/Tint'] = temperatures.tolist()
        d['teff0_table/NBT'] = [len(temperatures)]
        d['teff0_table/INT'] = [2]
        if self._include_gif:
            if self._isotopic_expansion:
                raise NotImplementedError('Isotopic expansion not yet implemented')
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

    def _createMF1(self, data, endf_parameters):
        """Creates MF=1 file of a thermal ENDF file.
           See ENDF-102, sect. 1.

        Parameters
        ----------
        endf_dict:
            endf-parserpy dictionary

        element : str
            Element to write

        data : dictionary
            Dictionary containing the nuclear data extrated by get_nuclear_data()

        self._verbosity : int
            Level of verbosity of the output (0: quiet)
        """
        awr = data.elements[self._sym].awr
        mat = self._mat
        za = data.elements[self._sym].za
        zsymam = data.elements[self._sym].zsymam
        self._endf_dict['1/451'] = {}
        d = self._endf_dict['1/451']
        d['MAT'] = mat
        d['ZA'] = za
        d['AWR'] = awr
        d['LRP'] = -1
        d['LFI'] = 0
        d['NLIB'] = endf_parameters.nlib
        d['NMOD'] = 0
        d['ELIS'] = 0
        d['LIS'] = 0
        d['LISO'] = 0
        d['STA'] = 0
        d['NFOR'] = 6
        d['AWI'] = 1.0
        d['EMAX'] = endf_parameters.emax
        d['LREL'] = endf_parameters.lrel
        d['NSUB'] = 12
        d['NVER'] = endf_parameters.nver
        d['TEMP'] = 0.0
        d['LDRV'] = 0
        d['HSUB/1'] = f'----{endf_parameters.libname:18s}MATERIAL {mat:4d}'.ljust(66)
        d['HSUB/2'] =  '-----THERMAL NEUTRON SCATTERING DATA'.ljust(66)
        d['HSUB/3'] =  '------ENDF-6 FORMAT'.ljust(66)
        d['NXC'] = 1
        d['ZSYMAM'] = zsymam.ljust(11)
        d['ALAB'] = endf_parameters.alab.ljust(11)
        d['AUTH'] = endf_parameters.auth.ljust(33)
        d['REF'] = endf_parameters.reference.ljust(21)
        now = datetime.now()
        months = ('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC')
        edate= f'EVAL-{months[now.month-1]}{now.strftime("%y")}'
        ddate= f'DIST-{months[now.month-1]}{now.strftime("%y")}'
        rdate= f'REV{endf_parameters.lrel:1d}-{months[now.month-1]}{now.strftime("%y")}'
        d['EDATE'] = edate
        d['DDATE'] = ddate
        d['RDATE'] = rdate
        d['ENDATE'] = endf_parameters.endate.ljust(8)
        description = []
        description.append(66*'*')
        description.append(''.ljust(66))
        description.append(' This file was converted from the following NCMAT [1] file:'.ljust(66))
        description.append(''.ljust(66))
        description.append(data.ncmat_fn.center(66))
        description.append(''.ljust(66))
        description.append(f' using NCrystal {nc_core.get_version()} and endf-parserpy {endf_parserpy.__version__} [2] with the '.ljust(66))
        description.append(' following options:'.ljust(66))
        description.append(''.ljust(66))
        for line in self._parameter_description:
            description.append(line)
        description.append(''.ljust(66))
        description.append(' Temperatures:'.ljust(66))
        for T in data.temperatures:
            description.append(f'       {T:.2f} K'.ljust(66))
        description.append(''.ljust(66))
        description.append('References:'.ljust(66))
        description.append('[1] https://github.com/mctools/ncrystal'.ljust(66))
        description.append('[2] https://endf-parserpy.readthedocs.io/en/latest/'.ljust(66))
        description.append(''.ljust(66))
        description.append(66*'*')
        description.append(''.ljust(66))
        description.append('Comments from NCMAT file:'.ljust(66))
        description.append(''.ljust(66))
        for line in data.ncrystal_comments.split('\n'):
            description.append(line.ljust(66))
        # description.append(''.ljust(66))
        description.append(66*'*')
        d['DESCRIPTION'] = {k:v for k, v in enumerate(description, start=1)}
        d['NWD'] = 5+len(description)
        d['MFx/1'] = 1
        d['MTx/1'] = 451
        d['NCx/1'] = 5
        d['MOD/1'] = d['NMOD']

    def write(self, endf_fn, force_save):
        import pathlib
        from ._common import write_text as nc_write_text
        if self._verbosity > 0:
            print(f'Write ENDF file {endf_fn}...')
        outfile = pathlib.Path(endf_fn)
        if outfile.exists() and not force_save:
            raise SystemExit('Error: output file already exists'
                             ' (run with --force to overwrite)')
        if not outfile.parent.is_dir():
            raise SystemExit('Error: output directory does not exist:'
                             f' { outfile.parent }')
        text = '\n'.join(self._parser.write(self._endf_dict, zero_as_blank=True))
        nc_write_text(outfile,text)

class EndfParameters():
    """Parameters for the ENDF-6 file

    Attributes
    ----------
    alab : string
        Mnemonic for the originating laboratory

    smin : float
        Minimum value of S(alpha, beta) to be stored in the file

    libname : string
        Name of the library

    nlib : int
        Library identifier (e.g. NLIB= 0 for ENDF/B).

    auth : string
        Author(s) name(s).

    reference : string
        Primary reference for the evaluation.

    emax : float
        Upper limit of the energy range for evaluation (eV).

    lrel : int
        Library release number.

    nver : int
        Library version number.

    endate: string
        Master File entry date in the form YYYYMMDD.
    """

    def __init__(self):
        self._alab = 'MyLAB'
        self._auth = 'NCrystal'
        self._reference = 'REFERENCE'
        self._nver = 1
        self._libname = 'MyLib'
        self._endate = 'YYYYMMDD'
        self._nlib = 0
        self._lrel = 0
        self._smin = 1e-100
        self._emax = 5.0
        self._lasym = 0 # Symmetric S(a,b) as default

    @property
    def alab(self):
        return self._alab
    @alab.setter
    def alab(self, x):
        self._alab = x

    @property
    def smin(self):
        return self._smin
    @smin.setter
    def smin(self, x):
        self._smin = x

    @property
    def libname(self):
        return self._libname
    @libname.setter
    def libname(self, x):
        self._libname = x

    @property
    def nlib(self):
        return self._nlib
    @nlib.setter
    def nlib(self, x):
        self._nlib = x

    @property
    def auth(self):
        return self._auth
    @auth.setter
    def auth(self, x):
        self._auth = x

    @property
    def reference(self):
        return self._reference

    @property
    def emax(self):
        return self._emax

    @property
    def lrel(self):
        return self._lrel

    @property
    def nver(self):
        return self._nver

    @property
    def endate(self):
        return self._endate

    @property
    def lasym(self):
        return self._lasym
    @lasym.setter
    def lasym(self, x):
        self._lasym = x

def ncmat2endf( ncmat_fn,
                name,
                endf_parameters,
                temperatures=(293.6,),
                mat_numbers=None,
                elastic_mode='scaled',
                include_gif=False,
                isotopic_expansion=False,
                vdoslux=3,
                force_save=False,
                verbosity=1):
    """Generates a set of ENDF-6 formatted files for a given NCMAT file.

    Parameters
    ----------
    ncmat_fn : str
        Filename of the ncmat file to convert

    temperatures : float or iterable of float
        Temperatures in Kelvin to generate the nuclear data

    mat_numbers : dict of str to int
        Material number for each element

    elastic_mode : str
        Treatment mode for the elastic component
        "greater" = only the greater ellastic component (coherent or incoherent) is saved
        "mixed"   = both the coherent and incoherent inelastic components are saved
        "scaled"  = for monoatomic scatterers, the major component is saved, scaled to the total bound XS
                    for polyatomic scatterers, coherent scattering for the whole system is assigned to the
                    atom with minimum incoherent cross section, and its incoherent contribution is distributed
                    among the other atoms

    include_gif: boolean
        Include the generalized information in MF=7/MT=451 in isotopes

    isotopic_expansion: boolean
        Expand the information in MF=7/MT=451 in isotopes

    vdoslux : integer
        Level of luxury to generate data in NCrystal

    verbosity : int
        Level of verbosity of the output (0: quiet)

    Returns
    -------
    file_names: list of (str, float)
        List of tuples contanining the ENDF-6 files and their fraction in the composition

    """
    if verbosity > 0 and endf_parameters.lasym > 0:
        print(f'Creating non standard S(a,b) with LASYM = {endf_parameters.lasym}')

    if type(temperatures) in [int, float]:
        temperatures = (temperatures,)

    if len(temperatures) > 1:
        warnings.warn('Multiple temperatures requested. Although this is supported, '
        +'it is not recommended because NCrystal generates a custom (alpha,beta) grid for each temperature. '
        +'The (alpha,beta) grid for first temperature will be used, and S(alpha, beta) for other temperatures will be interpolated.', stacklevel=2)
    if verbosity > 0:
        print('Get nuclear data...')

    data = NuclearData(ncmat_fn, temperatures, elastic_mode, vdoslux, verbosity)

    if mat_numbers is not None:
        n = len(mat_numbers)
        for frac, ad in data.composition:
            if ad.displayLabel() in mat_numbers.keys():
                n = n - 1
        assert n==0, 'Incorrect material number assignement'

    parameter_description = []
    parameter_description.append(f'  smin:{endf_parameters.smin}'.ljust(66))
    parameter_description.append(f'  emax:{endf_parameters.emax}'.ljust(66))
    parameter_description.append(f'  lasym:{endf_parameters.lasym}'.ljust(66))
    parameter_description.append(f'  include_gif:{include_gif}'.ljust(66))
    parameter_description.append(f'  vdoslux:{vdoslux}'.ljust(66))
    parameter_description.append(f'  isotopic_expansion:{isotopic_expansion}'.ljust(66))
    parameter_description.append(f'  elastic_mode:{elastic_mode}'.ljust(66))

    file_names = []
    for frac, ad in data.composition:
        sym = ad.displayLabel()
        mat = 999 if mat_numbers is None else mat_numbers[sym]
        endf_fn = f'tsl_{name}.endf' if sym == name else f'tsl_{sym}_in_{name}.endf'
        if data.elements[sym].sab_total is not None:
            endf_file = EndfFile(sym, data, mat, endf_parameters, include_gif, isotopic_expansion, parameter_description, verbosity)
            endf_file.write(endf_fn, force_save)
            file_names.append((endf_fn, frac))
        else:
            if verbosity > 0:
                print(f'Scattering kernel not available for: {endf_fn}')
    return(file_names)
