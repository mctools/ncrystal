
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

Module for creating a set of ENDF-6 thermal scattering files from a .ncmat
file. Most important parameters are supported as arguments, but additional
metadata parameters for the ENDF-6 format can be set by using the Python API and
pasing a custom EndfMetaData object or dictionary.

The module allows to handle multiple temperatures in one ENDF-6 file,
but this is not recommended, because NCrystal computes an optimal (alpha, beta)
 grid for each material and temperature, while the ENDF format imposes the
 same grid on all temperatures.

Ths module uses the endf-parserpy package from IAEA to format and check the
syntax of the ENDF-6 file.

G. Schnabel, D. L. Aldama, R. Capote, "How to explain ENDF-6 to computers:
A formal ENDF format description language", arXiv:2312.08249,
DOI:10.48550/arXiv.2312.08249

"""

available_elastic_modes = ('greater', 'scaled', 'mixed')

class EndfMetaData():
    """Optional MetaData Parameters for the ENDF-6 file describing the origin and
       authorship of the file.  For more information see the ENDF-6 format
       manual: https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf

    Attributes
    ----------
    alab : string
        Mnemonic for the originating laboratory

    libname : string
        Name of the nuclear data library

    nlib : int
        Nuclear data library identifier (e.g. NLIB= 0 for ENDF/B).

    auth : string
        Author(s) name(s).

    reference : string
        Primary reference for the evaluation.

    lrel : int
        Nuclear data library release number.

    nver : int
        Nuclear data library version number.

    endate: string
        Master File entry date in the form YYYYMMDD.

    edate: string Evaluation date in the form MMMYY. The special string "NOW"
        can be used to select the current date.

    ddate: string Distribution date in the form MMMYY. The special string "NOW"
        can be used to select the current date.

    rdate: string
        Revision date in the form MMMYY. The special string "NOW" can be used to
        select the current date.

    FIXME: Add parameter for material number assignment.

    #FIXME: Move these 3 parameters out of the metadata class:

    smin : float (fixme: should not be here??)
        Minimum value of S(alpha, beta) to be stored in the file

    emax : float
        Upper limit of the energy range for evaluation (eV).

    lasym : int
        Flag indicating whether an asymmetric S(a,b) is given.

    """

    def __init__(self):
        #fixme: to _ncmat2endf_impl.py. And encoding all data in dictionary will
        #make the rest of the class easier.
        self.__alab = 'MyLAB'
        self.__auth = 'NCrystal'
        self.__reference = 'REFERENCE'
        self.__nver = 1
        self.__libname = 'MyLib'
        self.__endate = ''
        self.__nlib = 0
        self.__lrel = 0
        self.__edate = 'MMMYY'
        self.__rdate = 'MMMYY'
        self.__ddate = 'MMMYY'

        #Ensure that all "NOW" dates will be evaluated from the same datetime
        #instance by capturing it here in init:
        from datetime import datetime
        now = datetime.now()
        now_MMMYY = now.strftime('%b%y').upper()
        def fmtdate(val, fmt):
            if not isinstance( val, str ):
                from .exceptions import NCBadInput
                raise NCBadInput('ENDF date value must be a string')
            assert fmt == 'MMMYY'
            if val.lower()=='now':
                return now_MMMYY
            if len(fmt) != len(val):
                from .exceptions import NCBadInput
                raise NCBadInput('ENDF date value is not in expected'
                                 f' format ("NOW" or "{fmt}")')
            return val
        self.__fmtdate = fmtdate

        #fixme: to be removed from this class:
        self.__smin = 1e-100
        self.__emax = 5.0
        self.__lasym = 0 # Symmetric S(a,b) as default

    @property
    def alab(self):
        return self.__alab

    def set_alab(self, x):
        self.__alab = x

    @property
    def libname(self):
        return self.__libname

    def set_libname(self, x):
        self.__libname = x

    @property
    def nlib(self):
        return self.__nlib

    def set_nlib(self, x):
        self.__nlib = x

    @property
    def auth(self):
        return self.__auth

    def set_auth(self, x):
        self.__auth = x

    @property
    def reference(self):
        return self.__reference

    def set_reference( self, x ):
        self.__reference = x

    @property
    def lrel(self):
        return self.__lrel

    def set_lrel(self,x):
        self.__lrel = x

    @property
    def nver(self):
        return self.__nver

    def set_nver(self, x):
        self.__nver = x

    def set_all_dates_as_now(self):
        self.set_edate('NOW')
        self.set_ddate('NOW')
        self.set_rdate('NOW')

    @property
    def endate(self):
        return self.__endate

    @property
    def edate(self):
        return self.__edate

    def set_edate(self, x):
        self.__edate = self.__fmtdate(x, 'MMMYY')

    @property
    def ddate(self):
        return self.__ddate

    def set_ddate(self, x):
        self.__ddate = self.__fmtdate(x, 'MMMYY')

    @property
    def rdate(self):
        return self.__rdate

    def set_rdate(self, x):
        self.__ddate = self.__fmtdate(x, 'MMMYY')

    def update_from_dict(self,d):
        """fixme: todo. Mention that the special key-value 'date':'NOW' or
        'date':'YYYYMMDD' can be used as a shorthand for setting all dates to
        that value (todo implement)"""
        if not hasattr(d,'items'):
            from .exceptions import NCBadInput
            raise NCBadInput('Parameter must be a dict')
        for k,v in d.items():
            setfct=getattr(self,f'set_{k}',None)
            if not setfct:
                from .exceptions import NCBadInput
                raise NCBadInput('fixme')
            setfct(v)

    #fixme: to be removed from this class:
    @property
    def smin(self):
        return self.__smin
    def set_smin(self, x):
        self.__smin = x
    @property
    def emax(self):
        return self.__emax
    def set_emax( self, x ):
        self.__emax = x
    @property
    def lasym(self):
        return self.__lasym
    def set_lasym(self, x):
        self.__lasym = x

def ncmat2endf( ncmat_cfg, *,
                material_name='NCrystalMaterial',
                endf_metadata=None,
                temperatures=None,
                mat_numbers=None,#fixme: to EndfMetaData
                elastic_mode='scaled',
                include_gif=False,
                isotopic_expansion=False,
                force_save=False,
                set_date_to_now=False,#fixme: remove
                verbosity=1):
    """Generates a set of ENDF-6 formatted files for a given NCMAT file.

    Parameters
    ----------
    ncmat_cfg : str
        Filename of the ncmat file to convert

    material_name : str
        name of the compound to be processed. ENDF files will be named
        tsl_element_in_name.endf for compounds or tsl_element.endf for
        elements. E.g. tsl_H_in_CH2.endf or tsl_Cu.endf

    endf_metadata : EndfMetaData object or a dictionary
        Metadata parameters for the ENDF file.
        https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf

    temperatures : int, float, tuple or list
        Temperature(s) in Kelvin to generate the nuclear data,
        in addition to the temperature defined in the cfg string.
        (The default temperature in the cfg string is 293.15 K)

    mat_numbers : dict of str to int
        Material number for each element

    elastic_mode : str
        Treatment mode for the elastic component
        "greater" = only the greater elastic component
                    (coherent or incoherent) is saved
        "mixed"   = both the coherent and incoherent inelastic
                    components are saved
        "scaled"  = for monoatomic scatterers, the major component is
                    saved, scaled to the total bound XS
                    for polyatomic scatterers, coherent scattering for
                    the whole system is assigned to the
                    atom with minimum incoherent cross section, and its
                    incoherent contribution is distributed
                    among the other atoms

    include_gif: boolean
        Include the generalized information in MF=7/MT=451 in isotopes

    isotopic_expansion: boolean
        Expand the information in MF=7/MT=451 in isotopes

    force_save: boolean
        Overwrite existing file if it already exists.

    set_date_to_now: boolean
        Set ENDF6 fields EDATE, DDATE and RDATE to current month and year.

    verbosity : int
        Level of verbosity of the output (0: quiet)

    Returns
    -------
    output_composition: list of (str, float)
        List of tuples contanining the ENDF-6 files and
        their fraction in the composition

    """
    from . import exceptions as nc_exceptions
    from . import core as nc_core
    from . import misc as nc_misc
    from . import cfgstr as nc_cfgstr
    from ._common import warn as ncwarn
    from ._common import print as ncprint
    from ._numpy import _ensure_numpy
    _ensure_numpy()

    if not endf_metadata:
        endf_metadata = EndfMetaData()
    elif not isinstance(endf_metadata,EndfMetaData):
        endf_metadata = EndfMetaData()
        endf_metadata.update_from_dict(endf_metadata)

    if set_date_to_now:
        endf_metadata.set_all_dates_as_now()

    if elastic_mode not in available_elastic_modes:
        raise nc_exceptions.NCBadInput(f'Elastic mode {elastic_mode}'
                                       f' not in ({available_elastic_modes})')
    info_obj = nc_core.createInfo(ncmat_cfg)
    if not info_obj.isSinglePhase():
        raise nc_exceptions.NCBadInput('Only single phase materials supported')
    if endf_metadata.lasym > 0:
        ncwarn( 'Creating non standard S(a,b)'
               f' with LASYM = {endf_metadata.lasym}')

    base_temp = info_obj.dyninfos[0].temperature
    if all(['temp' not in _ for _ in nc_cfgstr.decodeCfg(ncmat_cfg)['pars']]):
        ncwarn( 'Temperature not explicitly given in the cfg-string, '
               f' using T = {base_temp:.2f}')

    if temperatures is None:
        temperatures = tuple()
    else:
        if type(temperatures) in (int, float):
            temperatures = (temperatures,)
        elif type(temperatures) in (list, tuple):
            if any(type(T) not in (int, float) for T in temperatures):
                raise nc_exceptions.NCBadInput('Something wrong with the '
                                               'temperatures parameter: '
                                               f'({temperatures})')
            else:
                temperatures = tuple(float( T ) for T in temperatures )
        else:
            raise nc_exceptions.NCBadInput('temperatures parameter: '
                                           'should be a list or tuple '
                                           'of float or int')
    if base_temp in temperatures:
        raise nc_exceptions.NCBadInput('Repeated temperatures: '
                                       'temperatures parameter must not '
                                       'include the temperature defined '
                                       'in the cfg string')
    temperatures = sorted(temperatures + (base_temp,))
    if len(temperatures) > 1:
        ncwarn('Multiple temperatures requested. Although this is supported, '
               'it is not recommended because NCrystal generates '
               'a custom (alpha,beta) grid for each temperature. '
               'The (alpha,beta) grid for first temperature will '
               'be used, and S(alpha, beta) for other temperatures '
               'will be interpolated.')
    if any( T<=0 for T in temperatures ):
        raise nc_exceptions.NCBadInput('Non positive temperatures')

    if nc_core.createScatter(ncmat_cfg).isOriented():
        raise nc_exceptions.NCBadInput('Oriented materials cannot be '
                                       'represented in the ENDF format and '
                                       'are not supported' )
    scattering_components = nc_misc.detect_scattering_components(ncmat_cfg)
    if 'sans' in scattering_components:
        raise nc_exceptions.NCBadInput('SANS cannot be '
                                       'represented in the ENDF format and '
                                       'is not supported' )
    if 'inelas' not in scattering_components:
        raise nc_exceptions.NCBadInput('MF7/MT4 is mandatory in an ENDF file '
                                       'but no inelastic data found' )
    for di in info_obj.dyninfos:
        if type(di) not in (nc_core.Info.DI_VDOS, nc_core.Info.DI_VDOSDebye):
            raise NotImplementedError('Conversion supported only for VDOS'
                                      ' and VDOSDebye dyninfos')
    if (isotopic_expansion and not include_gif):
        raise nc_exceptions.NCBadInput( 'Isotopic expansion requires '
                                        'generalized information file, '
                                        'use --gif' )
    if (isotopic_expansion and include_gif):
        # TODO: remove when implemented
        raise NotImplementedError('Isotopic expansion not yet implemented')

    if type(verbosity) is not int:
        raise nc_exceptions.NCBadInput('Verbosity parameter'
                                       ' must be an integer')
    if verbosity > 0:
        ncprint('Get nuclear data...')
    from ._ncmat2endf_impl import NuclearData
    data = NuclearData(ncmat_cfg=ncmat_cfg, temperatures=temperatures,
                       elastic_mode=elastic_mode,
                       requested_emax=endf_metadata.emax,
                       verbosity=verbosity)

    if mat_numbers is not None:
        n = len(mat_numbers)
        for frac, ad in data.composition:
            if ad.elementName() in mat_numbers.keys():
                n = n - 1
        if n != 0:
            raise nc_exceptions.NCBadInput('Incorrect material number '
                                           'assignement')

    if material_name is None:
        material_name = 'UnknownCompound'

    output_composition = []
    from ._ncmat2endf_impl import EndfFile
    for frac, ad in data.composition:
        sym = ad.elementName()
        mat = 999 if mat_numbers is None else mat_numbers[sym]
        endf_fn = ( f'tsl_{material_name}.endf'
                   if sym == material_name
                   else f'tsl_{sym}_in_{material_name}.endf' )
        if data.elements[sym].sab_total is not None:
            endf_file = EndfFile(sym, data, mat, endf_metadata,
                                 include_gif=include_gif,
                                 isotopic_expansion=isotopic_expansion,
                                 verbosity=verbosity)
            endf_file.write(endf_fn, force_save)
            output_composition.append((endf_fn, frac))
        else:
            if verbosity > 0:
                ncprint(f'Scattering kernel not available for: {endf_fn}')

    if verbosity > 0:
        ncprint('Files created:')
        for fn, frac in output_composition:
            ncprint(f'  {fn} with fraction {frac}')

    return(output_composition)
