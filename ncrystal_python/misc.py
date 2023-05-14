
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2023 NCrystal developers                                   ##
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

"""Various utility functions and classes that do not fit elsewhere.

Of particular note, this module provides several type-erasure classes which can
facilitate the coding of rich and convenient interfaces able to accept many
input formats. For instance, the concept of "a material" in NCrystal could be
either an .ncmat file, a cfg-str, an in-memory string or an NCrystal.TextData
instance with NCMAT data, or even pre-loaded objects like Info, Scatter,
Absorption, or LoadedMaterial objects. Writing a function which works on "a
material" could therefore be quite cumbersome, but by passing the input through
the MaterialSource type erasure class from this module, all of these inputs can
be handled using a single interface.

Similar type-erasure classes exists for other types of data as well (e.g. VDOS
curves, generic text data, ...).

"""

#For client code checking all types of components:
standard_comp_types = ('coh_elas','incoh_elas','inelas','sans')

def detect_scattering_components( material ):
    """Detect which scattering components are present in the material.

    More precisely, the component names listed in the standard_comp_types
    variable will be checked, and a list of those present will be returned. Note
    they are suitable for subsequent usage in the "comp" cfg-string
    parameter.

    Note: This detection only works for material sources which are not
    preloaded. Trying to use this function with a preloaded material sources
    (e.g. an NCrystal.Scatter object) will result in an error.
    """
    from ._miscimpl import detect_scatcomps as _
    return _( standard_comp_types, MaterialSource( material ) )

class MaterialSource:
    """
    Generic type-erasure class for NCrystal materials in a general sence. It can
    be used to wrap cfg-strings, NCMAT data, TextData, LoadedMaterial objects,
    or even individual Info / Scatter / Absorption objects, etc.

    The most important feature of a MaterialSource object is that it has a
    .load(..) method which can be used to load the material into a
    LoadedMaterial object. Some MaterialSource objects will "preloaded", in the
    sense that they wrap already loaded objects (Info, Scatter, Absorption, or
    LoadedMaterial objects). MaterialSource objects that are NOT "preloaded"
    will have to do some work when the .load(..) method is invoked, but can on
    the other hand also offer a bit more flexibility (in particular the
    possiblity for loading with extra cfg-string parameters "appended").
    """

    def __init__( self, data, *, fmt = None, cfg_params = None ):
        """Wrap data object. If the data object is not a preloaded physics
        object, it is also allowed to specify cfg_params like
        cfg_params="temp=200K;vdoslux=2". Finally, the fmt string will be
        detected automatically, unless the fmt parameter is provided.
        """
        from . import _miscimpl
        self.__d = _miscimpl.matsrc_init( MaterialSource, lambda matsrc : matsrc.__d,
                                          data,cfg_params=cfg_params,fmt=fmt)

    @property
    def description( self ):
        """Returns a string with a MaterialSource description."""
        return self.__d['description']

    @property
    def plotlabel( self ):
        """Returns a string with either a specific plotlabel or otherwise just the usual description."""
        return self.__d.get('plotlabel') or self.__d['description']

    def __str__(self):
        return 'MaterialSource(%s)'%self.description

    def __repr__(self):
        return str(self)

    @property
    def is_preloaded( self ):
        """
        Check if material is already loaded into NCrystal physics
        objects. If so, the .load() function will be very cheap to call, but can not
        accept any extra_cfg_params.
        """
        return bool( self.__d.get('preloaded') )

    def load( self, extra_cfg_params = None, doInfo = True, doScatter = True, doAbsorption = True ):
        """
        Loads material (if not preloaded) and in any case returns a
        LoadedMaterial object. If material is preloaded, extra_cfg_params can
        NOT be provided, and the doInfo/doScatter/doAbsorption flags will have
        no effect.
        """
        ecp = (extra_cfg_params or '').strip()
        preloaded = self.__d.get('preloaded')
        if ecp and preloaded:
            from .exceptions import NCBadInput
            raise NCBadInput( 'Can not accept extra_cfg_params on a'
                              ' MaterialSource which is preloaded.' )
        if preloaded:
            return preloaded
        return self.__d['loadfct']( extra_cfg = extra_cfg_params,
                                    doInfo = doInfo,
                                    doScatter = doScatter,
                                    doAbsorption = doAbsorption )

class AnyTextData:
    """
    Immutable type-erasure class for text data (possibly named).
    Be aware, that this class will always keep the entire content in memory.
    """

    def __init__( self, data, is_path = None, name = None ):
        """
        Initialise from data, which can be either in-memory text data (str
        or NCrystal.TextData) or a file path (str or pathlib.Path) from which
        data will be loaded. Data of type bytes will simply be attempted decoded
        as a string.

        A string containing no new-line characters is assumed to be a file-path,
        otherwise it will be interpreted as in-memory data. Set is_path to True
        or False to override this logic.

        The "name" parameter can be used to specify or override the "name" of
        the file source.
        """
        if isinstance(data,AnyTextData):
            _ = ( name or data.__name ), data.__content, data.__keepalive
        else:
            from . import _miscimpl
            _ = _miscimpl.anytextdata_init( data, is_path=is_path, name=name )
        self.__name, self.__content, self.__keepalive = _

    @property
    def name( self ):
        """Access a name of the text data (could for instance be the file name). None if not available."""
        return self.__name

    @property
    def content( self ):
        """Access the contents (always type str)."""
        return self.__content

class AnyVDOS:
    """
    Immutable type-erasure class for phonon DOS (VDOS) curves.
    """

    def __init__( self, data, fmt = None, label = None ):
        """Initialise from data in various formats. Currently, input data can
        either be raw curves in the form of two arrays (egrid,density), DI_VDOS
        objects, or DI_VDOSDebye objects. Raw curves with egrid[0]>0 will be
        extended with a parabola towards (0,0). Furthermore egrid values
        consisting of just two points will be assumed to be a linearly spaced
        grid with those endpoints.

        The label parameter can be used to assign a label to the curves, which
        might for instance be used in utilities dumping or plotting the curves.

        Finally, some curves (notably those from DI_VDOS) will have the notion
        of "original" and processed curves. By default, the processed curves are
        displayed, although several methods have an "orig" parameter which can
        be set to True, in order to access the original curves instead.
        """
        from . import _miscimpl
        self.__d = _miscimpl._anyvdos_init( AnyVDOS, lambda x : x.__d,
                                            data, fmt, label )

    def __choice( self, orig ):
        return '_orig' if ( orig and 'dos_orig' in self.__d ) else ''

    def integral( self, *, orig = False ):
        """Get the integral of the VDOS curve. Note that this includes the area
        under the parabola between (0,0) and the first positive grid point. Any
        non-positive grid points are ignored for the purposes of calculating the
        integral."""
        return self.__d['derived']['dos%s_integral'%self.__choice(orig)]

    def egrid( self, *, orig = False ):
        """Access the energy grid values. Note that this is always the full
        energy grid, even if just a tuple "(emin,emax)" was used to construct
        the object.
        """
        choice = self.__choice(orig)
        eg = self.__d['derived'].get('dos%s_expandedegrid'%choice)
        return eg if eg is not None else self.__d['dos%s'%choice][0]

    def dos( self, *, orig = False, norm = True ):
        """Access the tabulated DOS values. Unless norm=False, this will
        actually be the normalised DOS, i.e. the DOS values divided by the
        integral (c.f. the .integral property)."""
        choice = self.__choice(orig)
        dos = self.__d['dos%s'%choice][1]
        return ( dos / self.__d['derived']['dos%s_integral'%choice] ) if norm else dos

    @property
    def has_orig( self ):
        """Whether or not this particular object has both processed and original curves."""
        return 'dos_orig' in self.__d

    @property
    def label( self ):
        """A string label describing the curve."""
        return self.__d['label']

    @property
    def debye_temperature( self ):
        """
        The Debye temperature in kelvin, if any is associated with the
        curve. Otherwise it will be None. This is usually only set for
        DI_VDOSDebye objects, or when a single floating point value (the Debye
        temperature itself) was used to initialise the AnyVDOS object.
        """
        return self.__d.get('debyetemp')
