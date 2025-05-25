
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
Utilities for converting NCrystal materials to ENDF. See the ncmat2endf
function in this module for more details.
"""

available_elastic_modes = ('greater', 'scaled', 'mixed')
default_smin_value = 1e-100
default_emax_value = 5.0

def ncmat2endf( ncmat_cfg, *,
                material_name='NCrystalMaterial',
                endf_metadata=None,
                temperatures=None,
                elastic_mode='scaled',
                include_gif=False,
                isotopic_expansion=False,
                force_save=False,
                smin=default_smin_value,
                emax=default_emax_value,
                lasym=0,
                verbosity = 1 ):
    """
    Generates a set of ENDF-6 thermal scattering files from an NCMAT
    cfg-string. Most important parameters are supported as arguments, but
    additional metadata parameters for the ENDF-6 format can be set by using
    the Python API and passing a custom EndfMetaData object or dictionary.

    The function allows to handle multiple temperatures in one ENDF-6 file, but
    this is not recommended, because NCrystal computes an optimal (alpha, beta)
    grid for each material and temperature, while the ENDF format imposes the
    same grid on all temperatures.

    This function uses the endf-parserpy package from IAEA to format and check
    the syntax of the ENDF-6 file.

    G. Schnabel, D. L. Aldama, R. Capote, "How to explain ENDF-6 to computers:
    A formal ENDF format description language", arXiv:2312.08249,
    DOI:10.48550/arXiv.2312.08249

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

    verbosity : int
        Level of verbosity of the output (0: quiet, max value: 3)

    Returns
    -------
    output_composition: list of (str, float)
        List of tuples contanining the ENDF-6 files and
        their fraction in the composition

    """

    # NOTICE ^^^^^^^^^^^
    #
    # When updating ncmat2endf there are 2 main doc texts that have to be
    # checked for updates:
    #
    #  1) The ncmat2endf function doc-string in ncmat2endf.py
    #  2) The ncmat2endf CLI --help text in _cli_ncmat2endf.py
    #

    from . import exceptions as nc_exceptions
    from . import core as nc_core
    from . import misc as nc_misc
    from . import cfgstr as nc_cfgstr
    from ._common import warn as ncwarn
    from ._common import print as ncprint
    from ._numpy import _ensure_numpy
    _ensure_numpy()

    if not isinstance(verbosity,int) or not ( 0<=verbosity<=3):
        raise nc_exceptions.NCBadInput('Invalid value of verbosity parameter'
                                       ' (expects value 0, 1, 2, or 3):'
                                       f' {verbosity}')

    if not endf_metadata:
        endf_metadata = EndfMetaData()
    elif not isinstance(endf_metadata,EndfMetaData):
        _ = EndfMetaData()
        _.update_from_dict(endf_metadata)
        endf_metadata = _

    if elastic_mode not in available_elastic_modes:
        raise nc_exceptions.NCBadInput(f'Elastic mode {elastic_mode}'
                                       f' not in ({available_elastic_modes})')
    info_obj = nc_core.createInfo(ncmat_cfg)
    if not info_obj.isSinglePhase():
        raise nc_exceptions.NCBadInput('Only single phase materials supported')
    if lasym > 0:
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
            raise nc_exceptions.NCBadInput('Conversion to ENDF supported only '
                                           'for VDOS and VDOSDebye dyninfos')
    if (isotopic_expansion and not include_gif):
        raise nc_exceptions.NCBadInput( 'Isotopic expansion requires '
                                        'generalized information file, '
                                        'use --gif' )
    if (isotopic_expansion and include_gif):
        raise nc_exceptions.NCBadInput('Isotopic expansion in conversion to'
                                       ' ENDF is not yet supported')

    if verbosity > 0:
        ncprint('Initialise nuclear data...')
    from ._ncmat2endf_impl import NuclearData
    data = NuclearData(ncmat_cfg=ncmat_cfg,
                       temperatures=temperatures,
                       elastic_mode=elastic_mode,
                       requested_emax=emax,
                       verbosity=verbosity)

    if endf_metadata.mat_numbers is not None:
        n = len(endf_metadata.mat_numbers)
        for frac, ad in data.composition:
            if ad.elementName() in endf_metadata.mat_numbers.keys():
                n = n - 1
        if n != 0:
            raise nc_exceptions.NCBadInput('Incorrect material number '
                                           'assignment')

    if material_name is None:
        material_name = 'UnknownCompound'

    output_composition = []
    from ._ncmat2endf_impl import EndfFile
    for frac, ad in data.composition:
        sym = ad.elementName()
        mat = ( 999 if not endf_metadata.mat_numbers
                else endf_metadata.mat_numbers.get(sym))
        if mat is None:
            raise nc_exceptions.NCBadInput('Incorrect material number '
                                           f'assignment for symbol "{sym}"')
        endf_fn = ( f'tsl_{material_name}.endf'
                   if sym == material_name
                   else f'tsl_{sym}_in_{material_name}.endf' )
        if data.elements[sym].sab_total is not None:
            endf_file = EndfFile(sym, data, mat, endf_metadata,
                                 include_gif=include_gif,
                                 isotopic_expansion=isotopic_expansion,
                                 smin=smin, emax=emax, lasym=lasym,
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


_metadata_definitions = dict(
    ALAB = dict( doc = "Mnemonic for the originating laboratory.",
                 defval = 'MyLAB' ),
    AUTH = dict( doc = "Author(s) name(s).",
                 defval = 'NCrystal' ),
    LIBNAME = dict( doc = "Name of the nuclear data library.",
                    defval = 'MyLib' ),
    NLIB = dict( datatype = int,
                 doc = "Nuclear data library identifier (e.g. NLIB=0 for ENDF/B).",
                 defval = 0 ),
    REFERENCE = dict( doc = "Primary reference for the evaluation.",
                      defval = 'REFERENCE' ),
    LREL = dict( datatype = int,
                 doc = "Nuclear data library release number.",
                 defval = 0 ),
    NVER = dict( datatype = int,
                 doc = "Nuclear data library version number.",
                 defval = 1 ),
    MAT_NUMBERS = dict( datatype = 'matnumbers',#dict of (str,int) to be more precise
                        doc = "MISSING? FIXME... should this parameter really be on this class?", #or at least remove the underscore for consistency?
                        defval = {} ),
    ENDATE = dict ( doc = 'Master File entry date in the form YYYYMMDD.',
                    defval = '' ),
    EDATE = dict( datatype = 'datestr',
                  doc = ('Evaluation date in the form MMMYY. The special string'
                         ' "NOW" can be used to select the current date.'),
                  defval = 'MMMYY' ),
    DDATE = dict( datatype = 'datestr',
                  doc = ('Distribution date in the form MMMYY. The special'
                         ' string "NOW" can be used to select the current'
                         ' date.'),
                  defval = 'MMMYY' ),
    RDATE = dict( datatype = 'datestr',
                  doc = ('Revision date in the form MMMYY. The special'
                         ' string "NOW" can be used to select the current'
                         ' date.'),
                  defval = 'MMMYY' ),
)

#def _add_attribs_EMD():
#    for k, v in _metadata_definitions.items():
#        EndfMetaData.
#
#_add_attribs_EMD()
#del _add_attribs_EMD

def _impl_emd_set( now_MMMYY, data, param, value,  ):
    k, v = param.upper(), value
    md = _metadata_definitions.get(k)
    if not md:
        from nc_exceptions import NCBadInput
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
            from .exceptions import NCBadInput
            raise NCBadInput('ENDF date value must be a string')
            if v.lower()=='now':
                v = now_MMMYY
            if len('MMMYY') != len(v):
                from .exceptions import NCBadInput
                raise NCBadInput('ENDF date value is not in expected'
                                 ' format, which is either special value'
                                 ' "NOW" or a date in the format "MMMYY"'
                                 ' (e.g. "Jun25").')
        data[k] = v
        return

    if isinstance(datatype,str) and datatype == 'matnumbers':
        if not hasattr(v,'items'):
            from .exceptions import NCBadInput
            raise NCBadInput('mat_numbers must be a dict')
        for kk,vv in v.items():
            if type(kk) is not str or type(vv) is not int:
                from .exceptions import NCBadInput
                raise NCBadInput('mat_numbers must be a dict from element',
                                 ' labels (str) to material values (int)' )
        data[k] = v
        return

    if not isinstance( v, datatype ):
        from nc_exceptions import NCBadInput
        raise NCBadInput(f'EndfMetaData parameter "{k}" data '
                         f'must be of type {datatype.__name__}')
    data[k] = v

class EndfMetaData():
    """Optional MetaData Parameters for the ENDF-6 file describing the origin
       and authorship of the file. For more information see the ENDF-6 format
       manual: https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf

    """

    def get( self, param ):
        """fixme"""
        v = self.__data.get(param.upper(),None)
        if v is None:
            from nc_exceptions import NCBadInput
            raise NCBadInput(f'Invalid EndfMetaData parameter "{param}"')
        return v

    def set_value( self, param, value ):
        """fixme... also None selects default"""
        _impl_emd_set( self.__now_MMMYY, self.__data, param, value )

    def update_from_dict( self, data ):
        if isinstance(data,EndfMetaData):
            return self.update_from_dict( data.__data )
        for k,v in data.items():
            self.set_value( k,v )

    def __init__(self, data = None):
        #fixme: to _ncmat2endf_impl.py. And encoding all data in dictionary will
        #make the rest of the class easier.

        #fixme: all set_... methods should double-check that types are correct.

        import copy
        self.__data = dict( (k,copy.deepcopy(v['defval']))#fixme: def mat_numbers is a non-immutable type, check that it does not make trouble
                            for k,v in _metadata_definitions.items() )

        #Ensure that all "NOW" dates will be evaluated from the same datetime
        #instance by capturing it here in init:
        from datetime import datetime
        now = datetime.now()
        self.__now_MMMYY = now.strftime('%b%y').upper()
#        def fmtdate(val, fmt):
#            if not isinstance( val, str ):
#                from .exceptions import NCBadInput
#                raise NCBadInput('ENDF date value must be a string')
#            assert fmt == 'MMMYY'
#            if val.lower()=='now':
#                return now_MMMYY
#            if len(fmt) != len(val):
#                from .exceptions import NCBadInput
#                raise NCBadInput('ENDF date value is not in expected'
#                                 f' format ("NOW" or "{fmt}")')
#            return val
#        self.__fmtdate = fmtdate
#
        if data:
            self.update_from_dict(data)

#        self.__alab = 'MyLAB'
#        self.__auth = 'NCrystal'
#        self.__reference = 'REFERENCE'
#        self.__nver = 1
#        self.__libname = 'MyLib'
#        self.__endate = ''
#        self.__nlib = 0
#        self.__lrel = 0
#        self.__mat_numbers = None
#        self.__edate = 'MMMYY'
#        self.__rdate = 'MMMYY'
#        self.__ddate = 'MMMYY'
#

    #fixme: repr should format json dict, __str__ instead?

    def as_json( self ):
        import json
        return json.dumps( self.__data )

    def __repr__(self):
        return '%s(%s)'%(
            self.__class__.__name__,
            self.as_json()
        )
#        return '%s({%s})'%(
#            self.__class__.__name__,
#            ', '.join(f'{repr(k.upper())}:{repr(v)}' for k,v in self.__data.items() )
#        )
#
    def __str__(self):
        return repr(self)


#        s = ('EndfMetaData object: '+
#            f'ALAB:{self.alab}, '+
#            f'AUTH:{self.auth}, '+
#            f'REFERENCE:{self.reference}, '+
#            f'NLIB:{self.nlib}, '+
#            f'NVER:{self.nver}, '+
#            f'LIBNAME:{self.libname}, '+
#            f'LREL:{self.lrel}, '+
#            f'MAT_NUMBERS:{self.mat_numbers}, '+
#            f'ENDATE:{self.endate}, '+
#            f'EDATE:{self.edate}, '+
#            f'RDATE:{self.rdate}, '+
#            f'DDATE:{self.ddate}')
#        return s

    @property
    def alab(self):
        """Mnemonic for the originating laboratory."""
        return self.get('alab')

    @property
    def libname(self):
        """Name of the nuclear data library."""
        return self.get('libname')

    @property
    def nlib(self):
        """Nuclear data library identifier (e.g. NLIB= 0 for ENDF/B)."""
        return self.get('nlib')

    @property
    def auth(self):
        """Author(s) name(s)."""
        return self.get('auth')

    @property
    def reference(self):
        """Primary reference for the evaluation."""
        return self.get('reference')

    @property
    def lrel(self):
        """Nuclear data library release number."""
        return self.get('lrel')

    @property
    def nver(self):
        """Nuclear data library version number."""
        return self.get('nver')

    def set_all_dates_as_now(self):
        """Set edate, ddate and rdate to the current date"""
        self.set_value('edate','NOW')
        self.set_value('ddate','NOW')
        self.set_value('rdate','NOW')

    @property
    def mat_numbers(self):
        return self.get('mat_numbers')

    #def set_mat_numbers(self, x):
    #    if not hasattr(x,'items'):
    #        from .exceptions import NCBadInput
    #        raise NCBadInput('mat_numbers must be a dict')
    #    for k,v in x.items():
    #        if type(k) is not str or type(v) is not int:
    #            from .exceptions import NCBadInput
    #            raise NCBadInput('mat_numbers must be a dict with str keys',
    #                             ' (element labels) and int values',
    #                             ' (material values)' )
    #    self.__mat_numbers = x

    @property
    def endate(self):
        """Master File entry date in the form YYYYMMDD."""
        #fixme: ^^^ mention if this is only for ENDF-B ?
        return self.get('endate')

    @property
    def edate(self):
        """Evaluation date in the form MMMYY. The special string "NOW"
        can be used to select the current date."""
        return self.get('edate')

#    def set_edate(self, x):
#        self.__edate = self.__fmtdate(x, 'MMMYY')

    @property
    def ddate(self):
        """Distribution date in the form MMMYY. The special string "NOW"
        can be used to select the current date."""
        return self.get('ddate')

    #def set_ddate(self, x):
    #    self.__ddate = self.__fmtdate(x, 'MMMYY')

    @property
    def rdate(self):
        """Revision date in the form MMMYY. The special string "NOW" can be used
        to select the current date."""
        return self.get('rdate')

    #def set_rdate(self, x):
    #    self.__ddate = self.__fmtdate(x, 'MMMYY')

    #def update_from_dict(self,d):
    #
    #    """fixme: todo. Mention that the special key-value 'date':'NOW' or
    #    'date':'YYYYMMDD' can be used as a shorthand for setting all dates to
    #    that value (todo implement)"""
    #    if not hasattr(d,'items'):
    #        from .exceptions import NCBadInput
    #        raise NCBadInput('Parameter must be a dict')
    #    for k,v in d.items():
    #        setfct=getattr(self,f'set_{k}',None)
    #        if not setfct:
    #            from .exceptions import NCBadInput
    #            raise NCBadInput(f'Key "{k}" in dict is not a'
    #                             ' supported EndfMetaData key')
    #        setfct(v)

    def get_param_and_docs(self):
        #Fixme: also mention default values?
        d = {}
        for k, v in _metadata_definitions.items():
            doc = getattr(EndfMetaData,k).__doc__
            doc = doc or 'missing'#FIXME remove this after adding doc-strings
            assert doc is not None
            d[k] = doc
        return d
