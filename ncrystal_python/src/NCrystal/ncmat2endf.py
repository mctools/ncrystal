
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

# Note: We must always null out S-values in ENDF files below a minimum, because
#       otherwise unwanted numerical issues might appear during subsequent
#       treatment in e.g. NJOY.
default_smin_value = 1e-100

#ENDF users expect an emax value of 5eV. Note that we will emit a warning if
#this value does not align with the Emax value that would have been normally
#produced from the cfg-string.
default_emax_value = 5.0

def ncmat2endf( ncmat_cfg, *,
                material_name = None,
                endf_metadata = None,
                othertemps = None,
                elastic_mode = 'scaled',
                include_gif = False,
                isotopic_expansion = False,
                force = False,
                smin = default_smin_value,
                emax = default_emax_value,
                lasym = 0,
                outdir = '.',
                verbosity = 1 ):
    """Generates a set of ENDF-6 thermal scattering files from an NCMAT
    cfg-string. Most important parameters are supported as arguments, but
    additional metadata parameters for the ENDF-6 format can be set by using
    the Python API and passing a custom EndfMetaData object or dictionary.

    The function allows to handle multiple temperatures in one ENDF-6 file via
    the othertemps parameter, but this is not recommended, because NCrystal
    computes an optimal (alpha, beta) grid for each material and temperature,
    while the ENDF format imposes the same grid on all temperatures.

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
        Name of the compound to be processed, to be used in the names of the
        produced ENDF files (e.g. the "CH2" in "tsl_H_in_CH2.endf"). For
        monoatomic materials this will default to nothing (e.g. "tsl_Al.endf"),
        and to "UnknownCompound" for polyatomic materials
        ("tsl_Al_in_UnknownCompound.endf").

 . ENDF files will be named
        tsl_element_in_name.endf for compounds or tsl_element.endf for
        elements. E.g. tsl_H_in_CH2.endf or tsl_Cu.endf

    endf_metadata : EndfMetaData object or a dictionary
        Metadata parameters for the ENDF file.
        https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf

    othertemps : int, float, tuple or list
        Temperature(s) in Kelvin to generate the nuclear data,
        in addition to the temperature defined in the cfg string.

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
        Include the generalized information in MF=7/MT=451 in isotopes.

    isotopic_expansion: boolean
        Expand the information in MF=7/MT=451 in isotopes.

    force: boolean
        Overwrite existing file if it already exists.

    outdir : str, path
        Directory for output files.

    verbosity : int
        Level of verbosity of the output (0: quiet, max value: 3).

    Returns
    -------
    summary: dict
        Returns a dictionary with information about the produced files, their
        fractions, temperature(s) [K] and density [g/cm3]. Example:

            {
                "density": 3.9912157584832233,
                "temperature": 293.15,
                "files": [
                    {
                        "file": "tsl_O_in_Al2O3.endf",
                        "fraction": 0.6,
                        "component": "O"
                    },
                    {
                        "file": "tsl_Al_in_Al2O3.endf",
                        "fraction": 0.4,
                        "component": "Al"
                    }
                ]
            }

        If running with more than one temperature, the "temperature" value will
        be a list of temperatures instead, while the density will always be that
        of the base temperature material.

    """

    # NOTICE ^^^^^^^^^^^
    #
    # When updating ncmat2endf there are 2 main doc texts that have to be
    # checked for updates:
    #
    #  1) The ncmat2endf function doc-string in ncmat2endf.py
    #  2) The ncmat2endf CLI --help text in _cli_ncmat2endf.py
    #

    from ._ncmat2endf_impl import _impl_ncmat2endf
    return _impl_ncmat2endf( ncmat_cfg = ncmat_cfg,
                             material_name = material_name,
                             endf_metadata = endf_metadata,
                             othertemps = othertemps,
                             elastic_mode = elastic_mode,
                             include_gif = include_gif,
                             isotopic_expansion = isotopic_expansion,
                             force = force,
                             smin = smin,
                             emax = emax,
                             lasym = lasym,
                             outdir = outdir,
                             verbosity = verbosity )

class EndfMetaData():
    """Optional MetaData Parameters for the ENDF-6 file describing the origin
       and authorship of the file. For more information see the ENDF-6 format
       manual: https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf

    """

    def __init__(self, data = None):
        """Initialise new EndfMetaData object, with default values of all
        parameters. If a data object is given in the form of a dictionary of
        (key,value) pairs, or another EndfMetaData object, the associated values
        will be updated accordingly.
        """
        from ._common import _datetime_now
        import copy
        from ._ncmat2endf_impl import _metadata_definitions
        self.__data = dict( (k,copy.deepcopy(v['defval']))
                            for k,v in _metadata_definitions.items() )
        self.__now_MMMYY = _datetime_now().strftime('%b%y').upper()
        if data:
            self.update_from_dict(data)

    @property
    def alab(self):
        """Mnemonic for the originating laboratory."""
        return self.get_value('alab')

    @property
    def libname(self):
        """Name of the nuclear data library."""
        return self.get_value('libname')

    @property
    def nlib(self):
        """Nuclear data library identifier (e.g. NLIB=0 for ENDF/B)."""
        return self.get_value('nlib')

    @property
    def auth(self):
        """Author(s) name(s)."""
        return self.get_value('auth')

    @property
    def reference(self):
        """Primary reference for the evaluation."""
        return self.get_value('reference')

    @property
    def lrel(self):
        """Nuclear data library release number."""
        return self.get_value('lrel')

    @property
    def nver(self):
        """Nuclear data library version number."""
        return self.get_value('nver')

    @property
    def matnum(self):
        """ENDF-6 MAT number assignment. Dictionary with string keys and int
           values that maps element symbols to MAT numbers, such as
           {"Zn":101,"O":102}. When setting this parameter it is also allowed
           to use a single string formatted like "Zn:101,O:102" See the report
           ENDF-102 Appendix C for numbering recommendations."""
        return self.get_value('matnum')

    @property
    def endate(self):
        """Master File entry date in the form YYYYMMDD (only used for ENDF/B
        releases)."""
        return self.get_value('endate')

    @property
    def edate(self):
        """Evaluation date in the form MMMYY. The special string "NOW"
        can be used to select the current date."""
        return self.get_value('edate')

    @property
    def ddate(self):
        """Distribution date in the form MMMYY. The special string "NOW"
        can be used to select the current date."""
        return self.get_value('ddate')

    @property
    def rdate(self):
        """Revision date in the form MMMYY. The special string "NOW" can be used
        to select the current date."""
        return self.get_value('rdate')

    def set_value( self, name, value ):
        """Set value of named parameter (naming is case insensitive). The
        special value None can be used to revert a parameter to its default
        value."""
        from ._ncmat2endf_impl import _impl_emd_set
        _impl_emd_set( self.__now_MMMYY, self.__data, name, value )

    def set_all_dates_as_now(self):
        """Set edate, ddate and rdate to the current date"""
        for k in 'edate', 'ddate', 'rdate':
            self.set_value(k,'NOW')

    def get_value( self, name ):
        """Get value of parameter (naming is case insensitive)."""
        v = self.__data.get(name.upper(),None)
        if v is None and name.upper() not in self.__data:
            from .exceptions import NCBadInput
            raise NCBadInput('Trying to read invalid EndfMetaData'
                             f' parameter "{name}"')
        return v

    def update_from_dict( self, data ):
        """Update parameters based on (key,value) pairs in dict. This is the
        same as calling .set_value(key,value) for all the entries. This method
        can also be called with data being an EndfMetaData object (essentially
        updating all values)."""
        if isinstance(data,EndfMetaData):
            return self.update_from_dict( data.to_dict() )
        for k,v in data.items():
            self.set_value( k,v )

    def to_dict( self ):
        """Returns dictionary with (key,value) pairs of all parameters and their
        values."""
        import copy
        return copy.deepcopy( self.__data )

    def to_json( self ):
        """Returns string containing a JSON encoded dictionary with (key,value)
        pairs of all parameters and their values.
        """
        import json
        return json.dumps( self.__data )

    def __repr__(self):
        return '%s(%s)'%( self.__class__.__name__, self.to_json() )

    def __str__(self):
        return repr(self)
