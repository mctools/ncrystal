
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

Module with core NCrystal functionality, including the OO classes (Info,
Scatter, Absorption, TextData, AtomData) and related factory methods.

"""

from .exceptions import ( NCrystalUserWarning, # noqa F401
                          NCException,
                          NCFileNotFound, # noqa F401
                          NCDataLoadError, # noqa F401
                          NCMissingInfo, # noqa F401
                          NCCalcError,
                          NCLogicError,
                          NCBadInput,
                          nc_assert )

from ._msg import _setDefaultPyMsgHandlerIfNotSet as _
_()
_=None

from ._chooks import _cstr2str, _get_raw_cfcts, _str2cstr, _get_build_namespace # noqa E402
from . import constants as _nc_constants # noqa E402
from ._numpy import _np,_ensure_numpy,_np_linspace # noqa E402
from . import _coreimpl as _impl # noqa E402
import enum as _enum # noqa E402
import ctypes as _ctypes # noqa E402
import weakref as _weakref # noqa E402
_rawfct = _get_raw_cfcts()

def get_version():
    """Get NCrystal version (same as NCrystal.__version__)"""
    from . import __version__ as _v
    return _v

def get_version_num():
    """Encode version in single integer (same as NCRYSTAL_VERSION C++ macro)
       for easier version comparison. This is also available as a variable named
       version_num in the main NCrystal module.
    """
    return sum(int(i)*j for i,j in zip(get_version().split('.'),(1000000,1000,1)))

def get_version_tuple():
    """Get NCrystal version as a tuple like (3,9,3). This is also available
       as a variable named version_tuple in the main NCrystal module.
    """
    return tuple( int(i) for i in get_version().split('.') )

def get_build_namespace():
    """If compiled with NCRYSTAL_NAMESPACE_PROTECTION, return the namespace here
       (will be an empty string in default installations).
    """
    return _get_build_namespace()

class RCBase:
    """Base class for all NCrystal objects"""
    def __init__(self, rawobj):
        """internal usage only"""
        self._rawobj = rawobj
        #do not ref here, since ncrystal_create_xxx functions in C-interface already did so.
        self._rawunref = _rawfct['ncrystal_unref']#keep fct reference
        self.__rawobj_byref = _ctypes.byref(rawobj)#keep byref(rawobj), since ctypes might
                                                  #disappear before __del__ is called.
    def __del__(self):
        if hasattr(self,'_rawunref') and self._rawunref:
            self._rawunref(self.__rawobj_byref)
    def refCount(self):
        """Access reference count of wrapped C++ object"""
        return _rawfct['ncrystal_refcount'](self._rawobj)

class AtomData(RCBase):
    """Class providing physical constants related to a particular mix of
    isotopes. This can be used to represent elements (i.e. all isotopes having
    same Z) in either natural or enriched form, but can also be used to
    represent atoms in doped crystals. E.g. if a small fraction (0.1%) of
    Cr-ions replace some Al-ions in a Al2O3 lattice, the AtomData could
    represent a mix of 0.1% Cr and 99.9% Al.
    """
    def __init__(self,rawobj):
        """internal usage only"""
        super(AtomData, self).__init__(rawobj)
        f=_rawfct['ncrystal_atomdata_getfields'](rawobj)
        self.__m = f['m']
        self.__incxs = f['incxs']
        self.__cohsl_fm = f['cohsl_fm']
        self.__absxs = f['absxs']
        self.__dl = f['dl']
        self.__descr = f['descr']
        self.__ncomp = f['ncomp']
        self.__z = f['z']
        self.__a = f['a']
        self.__b2f = (self.__m/(self.__m+_nc_constants.const_neutron_mass_amu))**2
        self.__comp = [None]*self.__ncomp
        self.__compalldone = (self.__ncomp==0)

    def averageMassAMU(self):
        """Atomic mass in Daltons (averaged appropriately over constituents)"""
        return self.__m
    def coherentScatLen(self):
        """Coherent scattering length in sqrt(barn)=10fm"""
        return self.__cohsl_fm*0.1#0.1 is fm/sqrt(barn)
    def coherentScatLenFM(self):
        """Coherent scattering length in fm"""
        return self.__cohsl_fm
    def coherentXS(self):
        """Bound coherent cross section in barn. Same as 4*pi*coherentScatLen()**2"""
        return _nc_constants.k4Pidiv100*self.__cohsl_fm**2
    def incoherentXS(self):
        """Bound incoherent cross section in barn"""
        return self.__incxs
    def scatteringXS(self):
        """Bound scattering cross section in barn (same as coherentXS()+incoherentXS())"""
        return self.__incxs+self.coherentXS()
    def captureXS(self):
        """Absorption cross section in barn"""
        return self.__absxs

    def freeScatteringXS(self):
        """Free scattering cross section in barn (same as freeCoherentXS()+freeIncoherentXS())"""
        return self.__b2f * self.scatteringXS()
    def freeCoherentXS(self):
        """Free coherent cross section in barn."""
        return self.__b2f * self.coherentXS()
    def freeIncoherentXS(self):
        """Free incoherent cross section in barn."""
        return self.__b2f * self.incoherentXS()

    def isNaturalElement(self):
        """Natural element with no composition."""
        return self.__z!=0 and self.__ncomp==0 and self.__a==0

    def isSingleIsotope(self):
        """Single isotope with no composition."""
        return self.__a!=0

    def isComposite(self):
        """Composite definition. See nComponents(), getComponent() and components property"""
        return self.__ncomp!=0

    def isElement(self):
        """If number of protons per nuclei is well defined. This is true for natural
           elements, single isotopes, and composites where all components
           have the same number of protons per nuclei."""
        return self.__z!=0

    def Z(self):
        """Number of protons per nuclei (0 if not well defined)."""
        return self.__z

    def elementName(self):
        """If Z()!=0, this returns the corresponding element name ('H', 'He', ...).
           Returns empty string when Z() is 0."""
        if not self.__z:
            return ''
        #NB: We are relying on natural elements to return their element names in
        #description(false). This is promised by a comment in NCAtomData.hh!
        if self.isNaturalElement():
            return self.__descr
        from .atomdata import elementZToName
        return elementZToName(self.__z)

    def A(self):
        """Number of nucleons per nuclei (0 if not well defined or natural element)."""
        return self.__a

    class Component:
        def __init__(self,fr,ad):
            """internal usage only"""
            self.__fr = fr
            self.__ad = ad
            assert not ad.isTopLevel()
        @property
        def fraction(self):
            """Fraction (by count) of component in mixture"""
            return self.__fr
        @property
        def data(self):
            """AtomData of component"""
            return self.__ad
        def __str__(self):
            return '%g*AtomData(%s)'%(self.__fr,self.__ad.description(True))
        def __repr__(self):
            return self.__str__()

    def nComponents(self):
        """Number of sub-components in a mixture"""
        return self.__ncomp
    def getComponent(self,icomponent):
        """Get component in a mixture"""
        c=self.__comp[icomponent]
        if c:
            return c
        rawobj_subc,fraction=_rawfct['ncrystal_atomdata_createsubcomp'](self._rawobj,icomponent)
        ad = AtomData(rawobj_subc)
        c = AtomData.Component(fraction,ad)
        self.__comp[icomponent] = c
        return c
    def getAllComponents(self):
        """Get list of all components"""
        if self.__compalldone:
            return self.__comp
        for i,c in enumerate(self.__comp):
            if not c:
                self.getComponent(i)
        self.__compalldone=True
        return self.__comp
    components = property(getAllComponents)

    def displayLabel(self):
        """Short label which unique identifies an atom role within a particular material."""
        return self.__dl

    def isTopLevel(self):
        """Whether or not AtomData appears directly on an Info object. If not,
         it will most likely either be a component (direct or indirect) of a top
         level AtomData object, or be taken from the composition list of a
         multi-phase object.
        """
        return bool(self.__dl)

    def description(self,includeValues=True):
        """Returns description of material as a string, with or without values."""
        if includeValues:
            zstr=' Z=%i'%self.__z if self.__z else ''
            astr=' A=%i'%self.__a if self.__a else ''
            _=(self.__descr,self.__cohsl_fm,self.coherentXS(),self.__incxs,
               self.__absxs,self.__m,zstr,astr)
            return'%s(cohSL=%gfm cohXS=%gbarn incXS=%gbarn absXS=%gbarn mass=%gamu%s%s)'%_
        return self.__descr

    def __str__(self):
        descr=self.description()
        return '%s=%s'%(self.__dl,descr) if self.__dl else descr

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def fmt_atomdb_str( mass, coh_scat_len, incoh_xs, abs_xs, sep=' ' ):
        """Takes mass value (amu), coherent scattering length (fm), incoherent cross
        section (barn), and absorption cross section @v_n=2200m/s (barn) and return
        data formatted in a string in a format like suitable for the the @ATOMDB
        section of ncmat files, or the atomdb cfg-string parameter. For instance,
        for Si such a data string might look like "28.09u 4.1491fm 0.004b 0.171b" By
        default the four fields are separated by spaces, but the sep parameter can
        be used to change this.
        """
        assert 0.0 < mass < 1e99
        assert -1e99 < coh_scat_len < 1e99
        assert 0.0 <= incoh_xs < 1e99
        assert 0.0 <= abs_xs < 1e99
        return f'{mass}u{sep}{coh_scat_len or 0:g}fm{sep}{incoh_xs or 0:g}b{sep}{abs_xs or 0:g}b'

    def to_atomdb_str( self, sep=' ' ):
        """Return data formatted in a string in a format like suitable for the
        the @ATOMDB section of ncmat files, or the atomdb cfg-string parameter
        (see the static method "fmt_atomdb_str" for details).
        """
        return AtomData.fmt_atomdb_str( mass=self.__m,
                                        coh_scat_len = self.__cohsl_fm,
                                        incoh_xs = self.__incxs,
                                        abs_xs = self.__absxs,
                                        sep = sep )

class StateOfMatter(_enum.Enum):
    """State of matter. Note that Solid's might be either amorphous or crystalline."""
    #NB: List here must be synchronized with list and values in NCInfo.hh:
    Unknown = 0
    Solid = 1
    Gas = 2
    Liquid = 3

class HKLInfoType(_enum.Enum):
    """Describes the kind of information about plane normals and Miller (hkl)
       indices available on each entry in the HKLList."""
    #NB: List here must be synchronized with list and values in NCInfoTypes.hh:
    SymEqvGroup = 0
    ExplicitHKLs = 1
    ExplicitNormals = 2
    Minimal = 3


class Info(RCBase):
    """Class representing information about a given material.

     Objects might represent either multi- or single phase
     materials. Multi-phase objects contain a list of phases (which might
     themselves be either single or multi-phase objects). Most other fields
     (structure, hkl lists, dynamics, etc.) are single-phase specific and will
     be unavailable on multiphase-phase objects. Exceptions are phase-structure
     information (todo) which is only available on multi-phase objects, and
     fields which are available for both multi- and single-phase objects such as
     density, composition, temperature, and state of matter (where such are well
     defined).

    """
    def __init__(self, cfgstr):
        """create Info object based on cfg-string (same as using createInfo(cfgstr))"""
        if isinstance(cfgstr,tuple) and len(cfgstr)==2 and cfgstr[0]=='_rawobj_':
            #Already got an ncrystal_info_t object:
            rawobj = cfgstr[1]
        else:
            rawobj = _rawfct['ncrystal_create_info'](_str2cstr(cfgstr))
        super(Info, self).__init__(rawobj)
        self.__dyninfo=None
        self.__atominfo=None
        self.__custom=None
        self.__atomdatas=[]
        self.__comp=None
        self._som=None
        self._nphases = int(_rawfct['ncrystal_info_nphases'](rawobj))
        assert self._nphases == 0 or self._nphases >= 2
        self._phases = tuple() if self._nphases == 0 else None

    def getUniqueID(self):
        """Unique identifier of object (UID)."""
        return _rawfct['infouid'](self._rawobj)
    uid = property(getUniqueID)

    def _getUnderlyingUniqueID(self):
        """Unique identifier of underlying object, which does not change on simple
        density or cfg-data overrides (expert usage only!)."""
        return _rawfct['infouid_underlying'](self._rawobj)

    def isSinglePhase(self):
        """Single phase object."""
        return self._nphases == 0

    def isMultiPhase(self):
        """Multi phase object."""
        return self._nphases != 0

    def __initPhases(self):
        assert self._phases is None and self._nphases > 1
        ll=[]
        for i in range(self._nphases):
            fraction = _ctypes.c_double()
            ph_info_raw = _rawfct['ncrystal_info_getphase'](self._rawobj,i,fraction)
            ph_info = Info( ('_rawobj_',ph_info_raw) )
            ll.append( ( float(fraction.value), ph_info ) )
        self._phases = tuple(ll)
        return self._phases

    def getPhases(self):
        """Daughter phases in a multi-phase object. Returns a list of fractions and Info
        objects of the daughter phases, in the format [(fraction_i,daughter_i),...]
        """
        return self.__initPhases() if self._phases is None else self._phases
    phases=property(getPhases)

    def _initComp(self):
        assert self.__comp is None
        nc = _rawfct['ncrystal_info_ncomponents'](self._rawobj)
        ll = []
        for icomp in range(nc):
            atomidx,fraction = _rawfct['ncrystal_info_getcomp'](self._rawobj,icomp)
            #NB: atomidx will be invalid in case of multiphase objects!
            if atomidx < 65535:
                #Most likely a single-phase object with valid atomidx, we can
                #use self._provideAtomData and share the AtomData objects also here on the python side:
                ll += [(fraction,self._provideAtomData(atomidx))]
            else:
                #Most likely a multi-phase object with invalid atomidx, we must
                #create new AtomData objects, based on ncrystal_create_component_atomdata:
                raw_ad = _rawfct['ncrystal_create_component_atomdata'](self._rawobj,icomp)
                obj = AtomData(raw_ad)
                assert not obj.isTopLevel()#does not appear directly on Info object
                ll += [(fraction,obj)]
        self.__comp = ll
        return self.__comp

    def stateOfMatter(self):
        """State of matter, i.e. Solid, Liquid, Gas, ... as per the options in the
        StateOfMatter class. Note that the .isCrystalline() method can be used
        to additionally distinguish between amorphous and crystalline
        solids. Return value is an enum object, whose .name() method can be used
        in case a string value is desired.
        """
        if self._som is None:
            self._som = StateOfMatter(_rawfct['ncrystal_info_getstateofmatter'](self._rawobj))
        return self._som

    def isCrystalline(self):
        """Whether or not object is crystalline (i.e. has unit cell structure or list of
        reflection planes)."""
        return self.hasStructureInfo() or self.hasAtomInfo() or self.hasHKLInfo()

    def hasComposition(self):
        """OBSOLETE FUNCTION (always available now)."""
        from ._common import warn
        warn('The .hasComposition method is obsolete'
             ' (it always returns True now).')
        return True

    def getComposition(self):
        """Get basic composition as list of (fraction,AtomData). For a single-phase
        object, the list is always consistent with AtomInfo/DynInfo (if
        present).
        """
        return self._initComp() if self.__comp is None else self.__comp
    composition=property(getComposition)

    def getFlattenedComposition( self,
                                 preferNaturalElements = True,
                                 naturalAbundProvider = None,
                                 asJSONStr=False ):
        """Break down the basic composition of the material into elements and
        isotopes. If an element only occurs as a natural element and has no
        specific isotopes, that element will be returned as a natural isotope
        unless preferNaturalElements=False. Generally, it is best if a
        naturalAbundProvider is given, since if a given Z value has a mix of
        isotopes and the natural elements, the natural element must always be
        broken up.

        If a naturalAbundProvider is given, it must be a function taking a Z
        value and returning the breakdown into isotopes,
        [(A1,frac1),...,(An,fracn)]. It can return None or an empty list to
        indicate missing information.

        Returns a list of (Z,<breakdown>) tuples, where <breakdown> is again a
        list of tuples of A-values and associated abundances (A=0 indicates
        natural elements). If asJSONStr=true, the data structure will be
        returned as a JSON-encoded string, instead of a Python dictionary.
        """
        _js = _rawfct['nc_info_getflatcompos'](self._rawobj,naturalAbundProvider,preferNaturalElements)
        import json
        return _js if asJSONStr else json.loads(_js)

    def dump(self,verbose=0):
        """Dump contained information to standard output. Use verbose argument to set
        verbosity level to 0 (minimal), 1 (middle), 2 (most verbose)."""
        _flush()
        _rawfct['ncrystal_dump_verbose'](self._rawobj,min(999,max(0,int(verbose))))
        _flush()

    def dump_str(self, verbose=0):
        """Return contained information as multi-line string. Use verbose argument to set
        verbosity level to 0 (minimal), 1 (middle), 2 (most verbose)."""
        return _rawfct['nc_dump_tostr'](self._rawobj,min(999,max(0,int(verbose))))

    def hasTemperature(self):
        """Whether or not material has a temperature available"""
        return _rawfct['ncrystal_info_gettemperature'](self._rawobj)>-1

    def getTemperature(self):
        """Material temperature (in kelvin)"""
        t=_rawfct['ncrystal_info_gettemperature'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasGlobalDebyeTemperature(self):
        """OBSOLETE FUNCTION: The concept of global versus per-element Debye
           temperatures has been removed. Please iterate over AtomInfo objects
           instead (see getAtomInfo() function) and get the Debye Temperature
           from those. This function will be removed in a future release.
        """
        from ._common import warn
        warn('The .hasGlobalDebyeTemperature method is obsolete'
             ' (it always returns False now).')
        return False

    def getGlobalDebyeTemperature(self):
        """OBSOLETE FUNCTION: The concept of global versus per-element Debye
           temperatures has been removed. Please iterate over AtomInfo objects
           instead (see getAtomInfo() function) and get the Debye Temperature
           from those. Calling this function will always result in an exception
           thrown for now, and the function will be removed in a future release..
        """
        raise NCLogicError('The concept of global Debye temperatures has been removed. Iterate over'
                           +' AtomInfo objects instead and get the Debye temperature values from those.')
        return None

    def hasAtomDebyeTemp(self):
        """Whether AtomInfo objects are present and have Debye temperatures available
        (they will either all have them available, or none of them will have
        them available).
        """
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[3]

    def hasDebyeTemperature(self):
        """Alias for hasAtomDebyeTemp()."""
        return self.hasAtomDebyeTemp()

    def hasAnyDebyeTemperature(self):
        """OBSOLETE FUNCTION which will be removed in a future release. Please
           call hasDebyeTemperature() instead.
        """
        from ._common import warn
        warn('The .hasAnyDebyeTemperature() method is deprecated.'
             ' Please use .hasAtomDebyeTemp() instead')
        return self.hasAtomDebyeTemp()

    def getDebyeTemperatureByElement(self,atomdata):
        """OBSOLETE FUNCTION which will be removed in a future release. Please access
           the AtomInfo objects instead and query the Debye temperature there.
        """
        from ._common import warn
        warn('The getDebyeTemperatureByElement method is deprecated.'
             ' Please access the AtomInfo objects instead and query'
             ' the Debye temperature there.')
        if atomdata.isTopLevel():
            for ai in self.atominfos:
                if atomdata is ai.atomData:
                    return ai.debyeTemperature
        raise NCBadInput('Invalid atomdata object passed to Info.getDebyeTemperatureByElement'
                         +' (must be top-level AtomData from the same Info object)')

    def hasDensity(self):
        """OBSOLETE FUNCTION (densities are now always available)."""
        from ._common import warn
        warn('The hasDensity() method is deprecated'
             ' (it always returns True now).')
        return True

    def getDensity(self):
        """Get density in g/cm^3. See also getNumberDensity()."""
        t=_rawfct['ncrystal_info_getdensity'](self._rawobj)
        nc_assert(t>0.0)
        return t
    density = property(getDensity)

    def hasNumberDensity(self):
        """OBSOLETE FUNCTION (densities are now always available)."""
        from ._common import warn
        warn('The hasNumberDensity() method is deprecated'
             ' (it always returns True now).')
        return True

    def getNumberDensity(self):
        """Get number density in atoms/angstrom^3. See also getDensity()."""
        t=_rawfct['ncrystal_info_getnumberdensity'](self._rawobj)
        nc_assert(t>0.0)
        return t
    numberdensity = property(getNumberDensity)

    @property
    def factor_macroscopic_xs( self ):
        """Factor needed to convert cross sections from (barns/atom) to inverse
        penetration depth (1/cm). This is actually just the numberdensity value,
        since (barns/atom) * ( atoms/Aa^3 ) = 1e-28m^2/1e-30m^3 = 1/cm.
        """
        return self.numberdensity

    def hasXSectAbsorption(self):
        """OBSOLETE FUNCTION"""
        from ._common import warn
        warn('The hasXSectAbsorption() method is deprecated'
             ' (it always returns True now).')
        return True

    def getXSectAbsorption(self):
        """Absorption cross section in barn (at 2200m/s)"""
        t=_rawfct['ncrystal_info_getxsectabsorption'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasXSectFree(self):
        """OBSOLETE FUNCTION"""
        from ._common import warn
        warn('The hasXSectFree() method is deprecated'
             ' (it always returns True now).')
        return True

    def getXSectFree(self):
        """Saturated (free) scattering cross section in barn in the high-E limit"""
        t=_rawfct['ncrystal_info_getxsectfree'](self._rawobj)
        nc_assert(t>-1)
        return t

    def getSLD(self):
        """Get scattering length density in 1e-6/Aa^2"""
        return _rawfct['ncrystal_info_getsld'](self._rawobj)
    sld = property(getSLD)

    def hasStructureInfo(self):
        """Whether or not material has crystal structure information available."""
        return bool(_rawfct['ncrystal_info_getstructure'](self._rawobj))
    def getStructureInfo(self):
        """Information about crystal structure."""
        d=_rawfct['ncrystal_info_getstructure'](self._rawobj)
        nc_assert(d)
        return d
    structure_info = property(getStructureInfo)

    def _provideAtomData(self,atomindex):
        if atomindex >= len(self.__atomdatas):
            if atomindex >= 65535:
                raise NCLogicError(f'Invalid atomindex ({atomindex}) provided to Info._provideAtomData')
            self.__atomdatas.extend([None,]*(atomindex+1-len(self.__atomdatas)))
        obj = self.__atomdatas[atomindex]
        if obj:
            return obj
        raw_ad = _rawfct['ncrystal_create_atomdata'](self._rawobj,atomindex)
        obj = AtomData(raw_ad)
        assert obj.isTopLevel()
        self.__atomdatas[atomindex] = obj
        return obj

    class AtomInfo:
        """Class with information about a particular atom in a unit cell, including the
        composition of atoms, positions, Debye temperature, and mean-squared-displacements.
        """

        def __init__(self,theinfoobj_wr,atomidx,n,dt,msd,pos):
            """For internal usage only."""
            assert dt is None or ( isinstance(dt,float) and dt > 0.0 )
            assert msd is None or ( isinstance(msd,float) and msd > 0.0 )
            self._info_wr = theinfoobj_wr
            self._atomidx,self.__n,self.__dt,self.__msd,=atomidx,n,dt,msd
            self.__pos = tuple(pos)#tuple, since it is immutable
            self.__atomdata = None
            self.__correspDI_wp = None

        def correspondingDynamicInfo(self):
            """
            Get corresponding DynamicInfo object from the same Info
            object. Returns None if Info object does not have dynamic info
            available
            """
            if self.__correspDI_wp is not None:
                if self.__correspDI_wp is False:
                    return None
                di = self.__correspDI_wp()
                nc_assert(di is not None,"AtomInfo.correspondingDynamicInfo can not"
                          +" be used after associated Info object is deleted")
                return di
            _info = self._info_wr()
            nc_assert(_info is not None,"AtomInfo.correspondingDynamicInfo can not"
                      +" be used after associated Info object is deleted")
            if not _info.hasDynamicInfo():
                self.__correspDI_wp = False
                return None
            for di in _info.dyninfos:
                if di._atomidx == self._atomidx:
                    self.__correspDI_wp = _weakref.ref(di)
                    return di
            nc_assert(False,"AtomInfo.correspondingDynamicInfo: inconsistent internal state (bug?)")
        dyninfo = property(correspondingDynamicInfo)

        @property
        def atomData(self):
            """Return AtomData object with details about composition and relevant physics constants"""
            if self.__atomdata is None:
                _info = self._info_wr()
                nc_assert(_info is not None,"AtomInfo.atomData can not be used after associated Info object is deleted")
                self.__atomdata = _info._provideAtomData(self._atomidx)
                assert self.__atomdata.isTopLevel()
            return self.__atomdata

        @property
        def count(self):
            """Number of atoms of this type per unit cell"""
            return self.__n

        @property
        def debyeTemperature(self):
            """The Debye Temperature of the atom (kelvin). Returns None if not available."""
            return self.__dt

        @property
        def meanSquaredDisplacement(self):
            """The mean-squared-displacement of the atom (angstrom^2), a.k.a. "U_iso". Returns None if not
               available.
            """
            return self.__msd
        msd=meanSquaredDisplacement#alias

        @property
        def positions(self):
            """List (tuple actually) of positions of this atom in the unit cell. Each
            entry is given as a tuple of three values, (x,y,z)"""
            return self.__pos

        @property
        def atomIndex(self):
            """Index of atom on this material"""
            return self._atomidx

        def __str__(self):
            ll=[str(self.atomData.displayLabel()),str(self.__n)]
            if self.__dt>0.0:
                ll.append('DebyeT=%gK'%self.__dt if self.__dt else 'DebyeT=n/a')
            if self.__msd>0.0:
                ll.append('MSD=%gAa^2'%self.__msd if self.__msd else 'MSD=n/a')
            ll.append('hasPositions=%s'%('yes' if self.__pos else 'no'))
            return 'AtomInfo(%s)'%(', '.join(ll))

    def hasAtomInfo(self):
        """Whether or no getAtomInfo()/atominfos are available"""
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[0]

    def hasAtomMSD(self):
        """Whether AtomInfo objects have mean-square-displacements (a.k.a. "U_iso") available"""
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[1]

    def hasAtomPositions(self):
        """OBSOLETE FUNCTION: AtomInfo objects now always have positions
           available. Returns same as hasAtomInfo(). Will be removed in a future
           release.
        """
        from ._common import warn
        warn('The hasAtomPositions() method is deprecated'
             ' (it always returns the same as .hasAtomInfo() now).')
        return self.hasAtomInfo()

    def hasPerElementDebyeTemperature(self):
        """OBSOLETE FUNCTION which will be removed in a future
           release. Please use hasAtomDebyeTemp() instead.
        """
        from ._common import warn
        warn('The hasPerElementDebyeTemperature() method is deprecated.'
             ' Please use the hasAtomDebyeTemp() method instead.')
        return self.hasAtomDebyeTemp()

    def getAtomInfo(self):
        """Get list of AtomInfo objects, one for each atom. Returns empty list if unavailable."""
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[2]
    atominfos = property(getAtomInfo)

    def __initAtomInfo(self):
        assert self.__atominfo is None
        natoms = _rawfct['ncrystal_info_natominfo'](self._rawobj)
        hasmsd = bool(_rawfct['ncrystal_info_hasatommsd'](self._rawobj))
        hasperelemdt=False
        ll=[]
        self_weakref = _weakref.ref(self)
        for iatom in range(natoms):
            atomidx,n,dt,msd = _rawfct['ncrystal_info_getatominfo'](self._rawobj,iatom)
            if dt:
                hasperelemdt=True
            assert hasmsd == (msd>0.0)
            pos=[]
            for ipos in range(n):
                pos.append( _rawfct['ncrystal_info_getatompos'](self._rawobj,iatom,ipos) )
            ll.append( Info.AtomInfo( self_weakref,atomidx, n,
                                      ( dt if ( dt and  dt>0.0) else None),
                                      (msd if (msd and msd>0.0) else None),
                                     pos) )
        self.__atominfo = ( natoms>0, hasmsd, ll, hasperelemdt )

    def hasHKLInfo(self):
        """Whether or not material has lists of HKL-plane info available"""
        return bool(_rawfct['ncrystal_info_nhkl'](self._rawobj)>-1)

    def nHKL(self):
        """Number of HKL planes available (grouped into families with similar
        d-spacing and f-squared)"""
        return int(_rawfct['ncrystal_info_nhkl'](self._rawobj))

    def hklDLower(self):
        """Lower d-spacing cutoff (angstrom)."""
        return float(_rawfct['ncrystal_info_hkl_dlower'](self._rawobj))

    def hklDUpper(self):
        """Upper d-spacing cutoff (angstrom)."""
        return float(_rawfct['ncrystal_info_hkl_dupper'](self._rawobj))

    def hklList(self,all_indices=False):
        """Iterator over HKL info, yielding tuples in the format
        (h,k,l,multiplicity,dspacing,fsquared). Running with all_indices=True to
        get the full list of hkl points in each group - in that case, h, k, and
        l will be numpy arrays of length multiplicity/2 (including just one of
        (h,k,l) and (-h,-k,-l) in the list).
        """
        nc_assert(self.hasHKLInfo())
        return _rawfct['iter_hkllist']( self._rawobj,
                                        all_indices = all_indices )

    def hklObjects( self ):
        """Iterator like .hklList, but with each entry returned as a single
        object, with the information accessible as (hopefully) more userfriendly
        friendly properties.

        Example usage:

        for e in info._hklObjects:
            #help( e );break #<- uncomment for usage info
            print( e )#<- a quick look
            print( e.hkl_label, e.mult, e.d, e.f2 )
            print( e.h, e.k, e.l )#all Miller indices as arrays.

        """
        from ._hklobjects import _iter_hklobjects
        for o in _iter_hklobjects(self):
            yield o

    def getBraggThreshold(self):
        """Get Bragg threshold in Aa (returns None if non-crystalline). This
        method is meant as a fast way to access the Bragg threshold without
        necessarily triggering a full initialisation of all HKL planes.
        """
        bt = float(_rawfct['ncrystal_info_braggthreshold'](self._rawobj))
        return bt if bt > 0.0 else None
    braggthreshold = property(getBraggThreshold)

    def hklInfoType(self):
        """What kind of information about plane normals and Miller indices are
        available in the hklList(). It is guaranteed to be the same for all
        HKLInfo entries, and will return "Minimal" when hklList() is present but
        empty. Like getBraggThreshold(), calling this method will not
        necessarily trigger a full initialisation of the hklList()."""
        return HKLInfoType(int(_rawfct['ncrystal_info_hklinfotype'](self._rawobj)))

    def hklIsSymEqvGroup(self):
        """Returns True if .hklInfoType() equals HKLInfoType.SymEqvGroup."""
        return self.hklInfoType() == HKLInfoType.SymEqvGroup

    def dspacingFromHKL(self, h, k, l): # noqa E741
        """Convenience method, calculating the d-spacing of a given Miller
        index. Calling this incurs the overhead of creating a reciprocal lattice
        matrix from the structure info."""
        return float(_rawfct['ncrystal_info_dspacing_from_hkl'](self._rawobj,h,k,l))

    class DynamicInfo:
        """Class representing dynamic information (related to inelastic scattering)
           about a given atom"""

        def __init__(self,theinfoobj_wr,fr,atomidx,tt):
            """internal usage only"""
            self._info_wr,self.__atomdata = theinfoobj_wr, None
            self.__fraction, self._atomidx, self.__tt = fr,atomidx,tt
            self.__correspAtomInfo_wp = None

        @property
        def _key( self ):
            i = self._info_wr()
            if not i:
                raise NCException('Dynamic info objects can not be used after the associated Info object is'
                                  ' deleted (the solution is normally to keep the Info object around'
                                  ' explicitly while you work on the dynamic info objects)')
            return i, (i._rawobj,self._atomidx)

        def correspondingAtomInfo(self):
            """Get corresponding AtomInfo object from the same Info object. Returns None if Info object does not have AtomInfo available"""
            if self.__correspAtomInfo_wp is not None:
                if self.__correspAtomInfo_wp is False:
                    return None
                ai = self.__correspAtomInfo_wp()
                nc_assert(ai is not None,"DynamicInfo.correspondingAtomInfo can not be used after associated Info object is deleted")
                return ai
            _info = self._info_wr()
            nc_assert(_info is not None,"DynamicInfo.correspondingAtomInfo can not be used after associated Info object is deleted")
            if not _info.hasAtomInfo():
                self.__correspAtomInfo_wp = False
                return None
            for ai in _info.atominfos:
                if ai._atomidx == self._atomidx:
                    self.__correspAtomInfo_wp = _weakref.ref(ai)
                    return ai
            nc_assert(False,"DynamicInfo.correspondingAtomInfo: inconsistent internal state (bug?)")
        atominfo = property(correspondingAtomInfo)

        @property
        def atomIndex(self):
            """Index of atom on this material"""
            return self._atomidx

        @property
        def fraction(self):
            """Atom fraction in material (all fractions must add up to unity)"""
            return self.__fraction

        @property
        def temperature(self):
            """Material temperature (same value as on associated Info object)"""
            return self.__tt

        @property
        def atomData(self):
            """Return AtomData object with details about composition and relevant physics constants"""
            if self.__atomdata is None:
                _info = self._info_wr()
                nc_assert(_info is not None,"DynamicInfo.atomData can not be used after associated Info object is deleted")
                self.__atomdata = _info._provideAtomData(self._atomidx)
                assert self.__atomdata.isTopLevel()
            return self.__atomdata

        def _np(self):
            _ensure_numpy()
            return _np

        def _copy_cptr_2_nparray(self,cptr,n):
            np = self._np()
            return np.copy(np.ctypeslib.as_array(cptr, shape=(n,)))

        def __str__(self):
            n=self.__class__.__name__
            if n.startswith('DI_'):
                n=n[3:]
            s=', %s'%self._extradescr() if hasattr(self,'_extradescr') else ''
            return ('DynamicInfo(%s, fraction=%.4g%%, type=%s%s)'%(self.atomData.displayLabel(),
                                                                   self.__fraction*100.0,
                                                                   n,s))
        def _plotlabel( self ):
            return self.atomData.displayLabel() or self.atomData.description(False)

    class DI_Sterile(DynamicInfo):
        """Class indicating atoms for which inelastic neutron scattering is absent
           or disabled."""
        pass

    class DI_FreeGas(DynamicInfo):
        """Class indicating atoms for which inelastic neutron scattering should be
           modelled as scattering on a free gas."""
        pass

    class DI_ScatKnl(DynamicInfo):
        """Base class indicating atoms for which inelastic neutron scattering will
           be, directly or indirectly, described by a scattering kernel,
           S(alpha,beta). This is an abstract class, and derived classes provide
           actual access to the kernels.
        """

        def __init__(self,theinfoobj_wr,fr,atomidx,tt):
            """internal usage only"""
            super(Info.DI_ScatKnl, self).__init__(theinfoobj_wr,fr,atomidx,tt)
            self.__lastknl,self.__lastvdoslux = None,None

        def _loadKernel( self, vdoslux = 3 ):
            import numbers
            assert isinstance(vdoslux,numbers.Integral) and 0<=vdoslux<=5
            vdoslux=int(vdoslux)
            if self.__lastvdoslux != vdoslux:
                _keepalive, key = self._key
                sugEmax,ne,na,nb,eptr,aptr,bptr,sabptr = _rawfct['ncrystal_dyninfo_extract_scatknl'](key,vdoslux)
                self.__lastvdoslux = vdoslux
                res={}
                assert ne>=0
                res['suggestedEmax'] = float(sugEmax)
                res['egrid'] = self._copy_cptr_2_nparray(eptr,ne) if ne > 0 else self._np().zeros(0)
                assert na>1 and nb>1
                res['alpha'] = self._copy_cptr_2_nparray(aptr,na)
                res['beta']  = self._copy_cptr_2_nparray(bptr,nb)
                res['sab']   = self._copy_cptr_2_nparray(sabptr,na*nb)
                res['temperature'] = self.temperature
                self.__lastknl = res
            assert self.__lastknl is not None
            return self.__lastknl

    class DI_ScatKnlDirect(DI_ScatKnl):
        """Pre-calculated scattering kernel which at most needs a (hidden) conversion to
           S(alpha,beta) format before it is available."""

        def __init__(self,theinfoobj_wr,fr,atomidx,tt):
            """internal usage only"""
            super(Info.DI_ScatKnlDirect, self).__init__(theinfoobj_wr,fr,atomidx,tt)

        def loadKernel( self ):
            """Prepares and returns the scattering kernel in S(alpha,beta) format.

            Note that the sab array is ordered so that
            S(alpha[i],beta[j])=sab[j*len(alpha)+i].
            """
            return self._loadKernel(vdoslux=3)#vdoslux value not actually used

        def plot_knl( self, **kwargs ):
            """Plot the scattering kernel using the NCrystal.plot.plot_knl
               function. Any kwargs are simply passed along.
            """
            from .plot import plot_knl
            plot_knl( self.loadKernel(), **kwargs )

    class DI_VDOS(DI_ScatKnl):
        """Solid state material with a phonon spectrum in the form of a Vibrational
        Density Of State (VDOS) parameterisation. This can be expanded into a
        full scattering kernel. How luxurious this expansion will be is
        controlled by an optional vdoslux parameter in the loadKernel call (must
        be integer from 0 to 5)
        """
        def __init__(self,theinfoobj_wr,fr,atomidx,tt):
            """internal usage only"""
            super(Info.DI_VDOS, self).__init__(theinfoobj_wr,fr,atomidx,tt)
            self.__vdosdata = None
            self.__vdosegrid_expanded = None
            self.__vdosorig = None

        def _extradescr(self):
            return 'npts=%i'%len(self.vdosOrigDensity())

        def vdosData(self):
            """Access the VDOS as ([egrid_min,egrid_max],vdos_density)"""
            if self.__vdosdata is None:
                _keepalive, key = self._key
                emin,emax,nd,dptr = _rawfct['ncrystal_dyninfo_extract_vdos'](key)
                vdos_egrid = (emin,emax)
                vdos_density = self._copy_cptr_2_nparray(dptr,nd)
                self.__vdosdata = (vdos_egrid,vdos_density)
            return self.__vdosdata

        def __loadVDOSOrig(self):
            if self.__vdosorig is None:
                _keepalive, key = self._key
                neg,egptr,nds,dsptr = _rawfct['ncrystal_dyninfo_extract_vdos_input'](key)
                self.__vdosorig = ( self._copy_cptr_2_nparray(egptr,neg),
                                    self._copy_cptr_2_nparray(dsptr,nds) )
            return self.__vdosorig

        def vdosOrigEgrid(self):
            """Access the original un-regularised VDOS energy grid"""
            return self.__loadVDOSOrig()[0]

        def vdosOrigDensity(self):
            """Access the original un-regularised VDOS energy grid"""
            return self.__loadVDOSOrig()[1]

        @property
        def vdos_egrid(self):
            """Access the VDOS energy grid as [egrid_min,egrid_max]"""
            return self.vdosData()[0]

        @property
        def vdos_egrid_expanded(self):
            """Access the egrid expanded into all actual egrid points"""
            if self.__vdosegrid_expanded is None:
                _ = self.vdosData()
                self.__vdosegrid_expanded = _np_linspace(_[0][0],_[0][1],len(_[1]))
            return self.__vdosegrid_expanded

        @property
        def vdos_density(self):
            """Access the VDOS density array"""
            return self.vdosData()[1]

        def loadKernel( self, vdoslux = 3 ):
            """Converts VDOS to S(alpha,beta) kernel with a luxury level given
            by the vdoslux parameter.

            Note that the sab array is ordered so that
            S(alpha[i],beta[j])=sab[j*len(alpha)+i].
            """
            return self._loadKernel(vdoslux=vdoslux)

        def analyseVDOS(self):
            """Same as running the global analyseVDOS function on the contained VDOS."""
            from .vdos import analyseVDOS
            return analyseVDOS(*self.vdos_egrid,self.vdos_density,
                               self.temperature,self.atomData.averageMassAMU())

        plot_knl = _impl.divdos_methods._plot_knl()
        plot_vdos = _impl.divdos_methods._plot_vdos()
        plot_Gn = _impl.divdos_methods._plot_Gn()
        extract_Gn = _impl.divdos_methods._extract_Gn()
        extract_custom_knl = _impl.divdos_methods._extract_custom_knl()

    class DI_VDOSDebye(DI_ScatKnl):
        """Similarly to DI_VDOS, but instead of using a phonon VDOS spectrum provided
           externally, an idealised spectrum is used for lack of better
           options. This spectrum is based on the Debye Model, in which the
           spectrum rises quadratically with phonon energy below a cutoff value,
           kT, where T is the Debye temperature
        """

        def __init__(self,theinfoobj_wr,fr,atomidx,tt):
            """internal usage only"""
            super(Info.DI_VDOSDebye, self).__init__(theinfoobj_wr,fr,atomidx,tt)
            self.__vdosdata = None
            self.__debyetemp = None
            self.__vdosegrid_expanded = None

        def vdosData(self):
            """Access the idealised VDOS as ([egrid_min,egrid_max],vdos_density)"""
            if self.__vdosdata is None:
                from .vdos import createVDOSDebye
                self.__vdosdata = createVDOSDebye(self.debyeTemperature())
            return self.__vdosdata

        def debyeTemperature(self):
            """The Debye temperature of the atom"""
            if self.__debyetemp is None:
                _keepalive, key = self._key
                self.__debyetemp = _rawfct['ncrystal_dyninfo_extract_vdosdebye'](key)
            return self.__debyetemp

        def _extradescr(self):
            return 'TDebye=%gK'%self.debyeTemperature()

        @property
        def vdos_egrid(self):
            """Access the VDOS energy grid as [egrid_min,egrid_max]"""
            return self.vdosData()[0]

        @property
        def vdos_egrid_expanded(self):
            """Access the egrid expanded into all actual egrid points"""
            if self.__vdosegrid_expanded is None:
                _ = self.vdosData()
                self.__vdosegrid_expanded = _np_linspace(_[0][0],_[0][1],len(_[1]))
            return self.__vdosegrid_expanded

        @property
        def vdos_density(self):
            """Access the VDOS density array"""
            return self.vdosData()[1]

        def loadKernel( self, vdoslux = 3 ):
            """Converts VDOS to S(alpha,beta) kernel with a luxury level given by the
               vdoslux parameter, which is similar to the vdoslux parameter used
               in DI_VDOS. Notice that the vdoslux parameter specified here on
               DI_VDOSDebye will be reduced internally by 3 (but not less than
               0), since the Debye model is anyway only a crude approximation
               and it accordingly does not need the same level of precise
               treatment as a full externally specified VDOS.

               Note that the sab array is ordered so that
               S(alpha[i],beta[j])=sab[j*len(alpha)+i].
            """
            return self._loadKernel(vdoslux=vdoslux)

        def analyseVDOS(self):
            """Same as running the global analyseVDOS function on the contained
               VDOS. Note that numbers returned here will be fully based on the
               actual tabulated VDOS curve, and values such as Debye temperature
               can therefore deviate slightly from the original value returned
               by .debyeTemperature().
            """
            from .vdos import analyseVDOS
            return analyseVDOS(*self.vdos_egrid,self.vdos_density,
                               self.temperature,self.atomData.averageMassAMU())


        plot_knl = _impl.divdos_methods._plot_knl()
        plot_vdos = _impl.divdos_methods._plot_vdos()
        plot_Gn = _impl.divdos_methods._plot_Gn()
        extract_Gn = _impl.divdos_methods._extract_Gn()
        extract_custom_knl = _impl.divdos_methods._extract_custom_knl()

    def hasDynamicInfo(self):
        """Whether or not dynamic information for each atom is present"""
        return int(_rawfct['ncrystal_info_ndyninfo'](self._rawobj))>0 if self.__dyninfo is None else bool(self.__dyninfo)

    def getDynamicInfoList(self):
        """Get list of DynamicInfo objects (if available). One for each atom."""
        if self.__dyninfo is None:
            ll = []
            self_weakref = _weakref.ref(self)
            for idx in range(int(_rawfct['ncrystal_info_ndyninfo'](self._rawobj))):
                fr,tt,atomidx,ditype = _rawfct['ncrystal_dyninfo_base']((self._rawobj,idx))
                args=(self_weakref,fr,atomidx,tt)
                if ditype==0:
                    di = Info.DI_Sterile(*args)
                elif ditype==1:
                    di = Info.DI_FreeGas(*args)
                elif ditype==2:
                    di = Info.DI_ScatKnlDirect(*args)
                elif ditype==3:
                    di = Info.DI_VDOS(*args)
                elif ditype==4:
                    di = Info.DI_VDOSDebye(*args)
                else:
                    raise NCLogicError('Unknown DynInfo type id (%i)'%ditype.value)
                ll.append( di )
            self.__dyninfo = ll
        return self.__dyninfo
    dyninfos = property(getDynamicInfoList)

    def findDynInfo( self, display_label ):
        """Look in the dyninfos list for an entry with the given
        display-label. Returns it if found, otherwise returns None."""
        for di in self.dyninfos:
            if di.atomData.displayLabel() == display_label:
                return di

    def findAtomInfo( self, display_label ):
        """Look in the atominfos list for an entry with the given
        display-label. Returns it if found, otherwise returns None."""
        for ai in self.atominfos:
            if ai.atomData.displayLabel() == display_label:
                return ai

    def getAllCustomSections(self):
        """Custom information for which the core NCrystal code does not have any
        specific treatment. This is usually intended for usage by developers adding new
        experimental physics models."""

        if self.__custom is None:
            self.__custom = _rawfct['ncrystal_info_getcustomsections'](self._rawobj)
        return self.__custom
    customsections = property(getAllCustomSections)

class Process(RCBase):
    """Base class for calculations of processes in materials.

    Note that kinetic energies are in electronvolt and direction vectors are
    tuples of 3 numbers.

    """
    def getCalcName(self):
        """Obsolete alias for getName"""
        from ._common import warn
        warn('The .getCalcName() method is deprecated.'
             ' Please use the .getName() method or the .name property instead.')
        return self.getName()

    def getName(self):
        """Process name"""
        return _cstr2str(_rawfct['ncrystal_name'](self._rawobj))
    name = property(getName)

    def getUniqueID(self):
        """UID of underlying ProcImpl::Process object."""
        return _rawfct['procuid'](self._rawobj)
    uid = property(getUniqueID)

    def domain(self):
        """Domain where process has non-vanishing cross section.

        Returns the domain as (ekin_low,ekin_high). Outside this range of
        neutron kinetic energy, the process can be assumed to have vanishing
        cross sections. Thus, processes present at all energies will return
        (0.0,infinity).

        """
        return _rawfct['ncrystal_domain'](self._rawobj)

    def isNull(self):
        """Domain might indicate that this is a null-process, vanishing everywhere."""
        elow,ehigh = self.domain()
        #checking for inf like the following to avoid depending on numpy or math
        #modules just for this:
        return ( elow >= ehigh or ( elow>1e99 and elow==float('inf') ) )

    def isNonOriented(self):
        """opposite of isOriented()"""
        return bool(_rawfct['ncrystal_isnonoriented'](self._rawobj))
    def isOriented(self):
        """Check if process is oriented and results depend on the incident direction of the neutron"""
        return not self.isNonOriented()
    def crossSection( self, ekin, direction ):
        """Access cross sections."""
        return _rawfct['ncrystal_crosssection'](self._rawobj,ekin, direction)
    def crossSectionIsotropic( self, ekin, repeat = None ):
        """Access cross sections (should not be called for oriented processes).

        For efficiency it is possible to provide the ekin parameter as a numpy
        array of numbers and get a corresponding array of cross sections
        back. Likewise, the repeat parameter can be set to a positive number,
        causing the ekin value(s) to be reused that many times and a numpy array
        with results returned.

        """
        return _rawfct['ncrystal_crosssection_nonoriented'](self._rawobj,ekin,repeat)

    def crossSectionNonOriented( self, ekin, repeat = None ):
        """Deprecated method. Please use the crossSectionIsotropic method
        instead."""
        from ._common import warn
        warn('The .crossSectionNonOriented method is deprecated.'
             ' Please use .crossSectionIsotropic or .xsect methods instead')
        return self.crossSectionIsotropic( ekin, repeat )

    def xsect(self,ekin=None,direction=None,wl=None,repeat=None):
        """Convenience function which redirects calls to either crossSectionIsotropic
        or crossSection depending on whether or not a direction is given. It can
        also accept wavelengths instead of kinetic energies via the wl
        parameter. The repeat parameter is currently only supported when
        direction is not provided.
        """
        ekin = Process._parseekin( ekin, wl )
        if direction is None:
            return self.crossSectionIsotropic( ekin, repeat )
        else:
            if repeat is None:
                return self.crossSection( ekin, direction )
            else:
                raise NCBadInput('The repeat parameter is not currently supported when the direction parameter is also provided.')

    @staticmethod
    def _parseekin(ekin,wl):
        if wl is None:
            if ekin is None:
                raise NCBadInput('Please provide either one of the "ekin" or "wl" parameters.')
            return ekin
        else:
            if ekin is not None:
                raise NCBadInput('Do not provide both "ekin" and "wl" parameters')
            return _nc_constants.wl2ekin(wl)

    def getSummary(self,short = False ):
        """By default access a high-level summary of the process in the form of
           a dictionary holding various information which is also available on
           the underlying C++ process object. If instead short==True, what is
           instead returned is simply a short process label along with a
           (recursive) list of sub-components and their scales (if
           appropriate). Finally, if short=='printable', the returned object
           will instead be a list of strings suitable for a quick printout (each
           string is one line of printout).

        """
        #Not caching, method is likely to be called sparringly.

        #short printable:
        if short == 'printable':
            toplbl,comps=self.getSummary(short=True)
            ll=[ toplbl]
            def add_lines( comps, indentlvl = 1 ):
                ncomps = len(comps)
                prefix = '   '*indentlvl
                for i,(scale,(lbl,subcomps)) in enumerate(comps):
                    smb = r'\--' if i+1==ncomps else '|--'
                    scale_str = '' if scale==1.0 else f'{scale:g} * '
                    ll.append(f'{prefix}{smb} {scale_str}{lbl}')
                    if subcomps:
                        add_lines( subcomps, indentlvl + 1 )
            if comps:
                add_lines( comps )
            return ll

        d=_rawfct['nc_dbg_proc'](self._rawobj)
        #full:
        if not short:
            return d
        #short:
        def fmt_lbl(proc):
            name = proc['name']
            summarystr = proc['specific'].get('summarystr','')
            return f'{name}({summarystr})' if summarystr else name
        def extract_subcomponents(proc):
            subprocs = proc.get('specific',{}).get('components',[])
            return (fmt_lbl(proc),list( (scl, extract_subcomponents(sp)) for scl,sp in subprocs ))
        return extract_subcomponents(d)

    def dump(self,prefix=''):
        """
        Prints a quick high level summary of the process. What is printed is in
        fact the lines resulting from a call to
        self.getSummary(short='printable'), with an optional prefix prepended to
        each line.
        """
        from ._common import print
        _flush()
        print(self.dump_str(prefix=prefix),end='')
        _flush()

    def dump_str(self, prefix=''):
        """
        The string-returning sibling of .dump(..). Returns a quick high level
        summary of the process as a multi-line string. What is printed is in
        fact the lines resulting from a call to
        self.getSummary(short='printable'), with an optional prefix prepended to
        each line.
        """
        return prefix+f'\n{prefix}'.join(self.getSummary(short='printable'))+'\n'

    def plot(self, *args, **kwargs ):
        """Convenience method for plotting cross sections. This is the same as
        NCrystal.plot.plot_xsect(material=self,*args,**kwargs), so refer to that
        function for information about allowed arguments."""
        from .plot import plot_xsect
        return plot_xsect( self, *args, **kwargs )

class Absorption(Process):
    """Base class for calculations of absorption in materials"""

    def __init__(self, cfgstr):
        """create Absorption object based on cfg-string (same as using createAbsorption(cfgstr))"""
        if isinstance(cfgstr,tuple) and len(cfgstr)==2 and cfgstr[0]=='_rawobj_':
            #Cloning:
            rawobj_abs = cfgstr[1]
        else:
            rawobj_abs = _rawfct['ncrystal_create_absorption'](_str2cstr(cfgstr))
        self._rawobj_abs = rawobj_abs
        rawobj_proc = _rawfct['ncrystal_cast_abs2proc'](rawobj_abs)
        super(Absorption, self).__init__(rawobj_proc)

    def clone(self):
        """Clone object. The clone will be using the same physics models and sharing any
         read-only data with the original, but will be using its own private copy of any
         mutable caches. All in all, this means that the objects are safe to use
         concurrently in multi-threaded programming, as long as each thread gets
         its own clone. Return value is the new Absorption object.
        """
        newrawobj = _rawfct['ncrystal_clone_absorption'](self._rawobj_abs)
        return Absorption( ('_rawobj_',newrawobj) )

class Scatter(Process):

    """Base class for calculations of scattering in materials.

    Note that kinetic energies are in electronvolt and direction vectors are
    tuples of 3 numbers.

    """

    def __init__(self, cfgstr):
        """create Scatter object based on cfg-string (same as using createScatter(cfgstr))"""
        if isinstance(cfgstr,tuple) and len(cfgstr)==2 and cfgstr[0]=='_rawobj_':
            #Already got an ncrystal_scatter_t object:
            self._rawobj_scat = cfgstr[1]
        else:
            self._rawobj_scat = _rawfct['ncrystal_create_scatter'](_str2cstr(cfgstr))
        rawobj_proc = _rawfct['ncrystal_cast_scat2proc'](self._rawobj_scat)
        super(Scatter, self).__init__(rawobj_proc)


    def clone(self,rng_stream_index=None,for_current_thread=False):
        """Clone object. The clone will be using the same physics models and sharing any
        read-only data with the original, but will be using its own private copy
        of any mutable caches and will get an independent RNG stream. All in
        all, this means that the objects are safe to use concurrently in
        multi-threaded programming, as long as each thread gets its own
        clone. Return value is the new Scatter object.

        If greater control over RNG streams are needed, it is optionally allowed
        to either set rng_stream_index to a non-negative integral value, or set
        for_current_thread=True.

        If rng_stream_index is set, the resulting object will use a specific
        rngstream index. All objects with the same indeed will share the same
        RNG state, so a sensible strategy is to use the same index for all
        scatter objects which are to be used in the same thread:

        If setting for_current_thread=True, the resulting object will use a
        specific rngstream which has been set aside for the current thread. Thus
        this function can be called from a given work-thread, in order to get
        thread-safe scatter handle, with all objects cloned within the same
        thread sharing RNG state.

        """
        if rng_stream_index is not None:
            if for_current_thread:
                raise NCBadInput('Scatter.clone(..): do not set both rng_stream_index and for_current_thread parameters')
            import numbers
            if not isinstance(rng_stream_index, numbers.Integral) or not 0 <= rng_stream_index <= 4294967295:
                raise NCBadInput('Scatter.clone(..): rng_stream_index must be integral and in range [0,4294967295]')
            newrawobj = _rawfct['ncrystal_clone_scatter_rngbyidx'](self._rawobj_scat,int(rng_stream_index))
        elif for_current_thread:
            newrawobj = _rawfct['ncrystal_clone_scatter_rngforcurrentthread'](self._rawobj_scat)
        else:
            newrawobj = _rawfct['ncrystal_clone_scatter'](self._rawobj_scat)
        return Scatter( ('_rawobj_',newrawobj) )

    def sampleScatter( self, ekin, direction, repeat = None ):
        """Randomly generate scatterings.

        Assuming a scattering took place, generate final state of neutron based
        on current kinetic energy and direction. Returns
        tuple(ekin_final,direction_final) where direct_final is itself a tuple
        (ux,uy,uz). The repeat parameter can be set to a positive number,
        causing the scattering to be sampled that many times and numpy arrays
        with results returned.

        """
        return _rawfct['ncrystal_samplesct'](self._rawobj_scat,ekin,direction,repeat)


    def sampleScatterIsotropic( self, ekin, repeat = None ):
        """Randomly generate scatterings (should not be called for oriented processes).

        Assuming a scattering took place, generate final state of
        neutron. Returns tuple(ekin_final,mu) where mu is the cosine of the
        scattering angle. For efficiency it is possible to provide the ekin
        parameter as a numpy array of numbers and get corresponding arrays of
        energies and mu back. Likewise, the repeat parameter can be
        set to a positive number, causing the ekin value(s) to be reused that
        many times and numpy arrays with results returned.

        """
        return _rawfct['ncrystal_samplesct_iso'](self._rawobj_scat,ekin,repeat)

    def generateScattering( self, ekin, direction, repeat = None ):
        """WARNING: Deprecated method. Please use the sampleScatter method instead.

        Randomly generate scatterings.

        Assuming a scattering took place, generate energy transfer (delta_ekin)
        and new neutron direction based on current kinetic energy and direction
        and return tuple(new_direction,delta_ekin). The repeat parameter can be
        set to a positive number, causing the scattering to be sampled that many
        times and numpy arrays with results returned.

        """
        from ._common import warn
        warn('The .generateScattering method is deprecated.'
             ' Please use .sampleScatter or .scatter instead')
        return _rawfct['ncrystal_genscatter'](self._rawobj_scat,ekin,direction,repeat)

    def generateScatteringNonOriented( self, ekin, repeat = None ):
        """WARNING: Deprecated method. Please use the sampleScatterIsotropic method instead.

        Randomly generate scatterings (should not be called for oriented processes).

        Assuming a scattering took place, generate energy transfer (delta_ekin)
        and scatter angle in radians and return tuple(scatter_angle,delta_ekin)
        (this method should not be invoked on oriented processes).  For
        efficiency it is possible to provide the ekin parameter as a numpy array
        of numbers and get corresponding arrays of angles and energy transfers
        back. Likewise, the repeat parameter can be set to a positive number,
        causing the ekin value(s) to be reused that many times and numpy arrays
        with results returned.

        """
        from ._common import warn
        warn('The .generateScatteringNonOriented method is deprecated.'
             ' Please use .sampleScatterIsotropic or .scatter instead')
        return _rawfct['ncrystal_genscatter_nonoriented'](self._rawobj_scat,ekin,repeat)

    def scatter(self,ekin=None,direction=None,wl=None,repeat=None):
        """Convenience function which redirects calls to either
        sampleScatterIsotropic or sampleScatter depending on whether
        or not a direction is given. It can also accept wavelengths instead of
        kinetic energies via the wl parameter.
        """
        ekin = Process._parseekin( ekin, wl )
        return ( self.sampleScatterIsotropic( ekin, repeat )
                 if direction is None
                 else self.sampleScatter( ekin, direction, repeat ) )

    def genscat(self,ekin=None,direction=None,wl=None,repeat=None):
        """WARNING: Deprecated method. Please use the "scatter" method instead.

        Convenience function which redirects calls to either
        generateScatteringNonOriented or generateScattering depending on whether
        or not a direction is given. It can also accept wavelengths instead of
        kinetic energies via the wl parameter.
        """
        from ._common import warn
        warn('The .genscat method is deprecated.'
             ' Please use .scatter instead')
        ekin = Process._parseekin( ekin, wl )
        return ( self.generateScatteringNonOriented( ekin, repeat )
                 if direction is None
                 else self.generateScattering( ekin, direction, repeat ) )

    def rngSupportsStateManipulation(self):
        """Query whether associated RNG stream supports state manipulation"""
        return bool(_rawfct['ncrystal_rngsupportsstatemanip_ofscatter'](self._rawobj_scat))

    def getRNGState(self):
        """Get current RNG state (as printable hex-string with RNG type info
           embedded). This function returns None if RNG stream does not support
           state manipulation
        """
        return _rawfct['nc_getrngstate_scat'](self._rawobj_scat)

    def setRNGState(self,state):
        """Set current RNG state.

           Note that setting the rng state will affect all objects sharing the
           RNG stream with the given scatter object (and those subsequently cloned
           from any of those).

           Note that if the provided state originates in (the current version
           of) NCrystal's builtin RNG algorithm, it can always be used here,
           even if the current RNG uses a different algorithm (it will simply be
           replaced). Otherwise, a mismatch of RNG stream algorithms will result
           in an error.
        """
        _rawfct['ncrystal_setrngstate_ofscatter']( self._rawobj_scat,
                                                   _str2cstr(state) )


def createInfo(cfgstr):
    """Construct Info object based on provided configuration (using available factories)"""
    return Info(cfgstr)

def createScatter(cfgstr):
    """Construct Scatter object based on provided configuration (using available factories)"""
    return Scatter(cfgstr)

def createScatterIndependentRNG(cfgstr,seed = 0):
    """Construct Scatter object based on provided configuration (using available
    factories) and with its own independent RNG stream (using the builtin RNG
    generator and the provided seed)"""
    rawobj = _rawfct['ncrystal_create_scatter_builtinrng'](_str2cstr(cfgstr),seed)
    return Scatter(('_rawobj_',rawobj))

def createAbsorption(cfgstr):
    """Construct Absorption object based on provided configuration (using available factories)"""
    return Absorption(cfgstr)

def createLoadedMaterial(cfgstr):
    """Create and return LoadedMaterial object. A "loaded material" is simply a
    convenience object which combines loaded Info, Scatter, and Absorption
    objects of the same material.

    This function does the same as calling the LoadedMaterial constructor
    directly, and exists purely for consistency.
    """
    return LoadedMaterial(cfgstr)

load = createLoadedMaterial#convenience alias

class LoadedMaterial:
    """This convenience class combines loaded Info, Scatter, and Absorption
    objects of the same material."""

    @staticmethod
    def fromExistingObjects(info=None,scatter=None,absorption=None):
        return LoadedMaterial( ('__fromexistingobjects__',dict(info=info,scatter=scatter,absorption=absorption)) )

    def __init__(self,cfgstr):
        """Instantiate from a cfg-string which will be passed to the createInfo,
        createScatter, and createAbsorption functions.

        As a special case, one can also pass a TextData object to this function,
        and it will be loaded as by calling the directLoad(..)
        function. Likewise, python strings beginning with 'NCMAT' and containing
        at least one newline, will also be assumed to be raw NCMAT data and
        loaded accordingly. For more control, however, the .directLoad function
        is recommended.
        """
        if isinstance(cfgstr,tuple) and len(cfgstr)==2 and cfgstr[0] == '__fromexistingobjects__':
            self.__i,self.__s,self.__a = cfgstr[1]['info'],cfgstr[1]['scatter'],cfgstr[1]['absorption']
            assert self.__i is None or isinstance(self.__i,Info)
            assert self.__s is None or isinstance(self.__s,Scatter)
            assert self.__a is None or isinstance(self.__a,Absorption)
        elif ( isinstance(cfgstr,TextData)
               or ( isinstance(cfgstr,str) and cfgstr.startswith('NCMAT') and '\n' in cfgstr ) ):
            m = directLoad( cfgstr, dtype = cfgstr.dataType if isinstance(cfgstr,TextData) else 'ncmat' )
            self.__i = m.info
            self.__s = m.scatter
            self.__a = m.absorption
        else:
            self.__i = createInfo(cfgstr)
            self.__s = createScatter(cfgstr)
            self.__a = createAbsorption(cfgstr)

    @property
    def info(self):
        """Info object (None if not present)."""
        return self.__i

    @property
    def scatter(self):
        """Scatter object (None if not present)."""
        return self.__s

    @property
    def absorption(self):
        """Absorption object (None if not present)."""
        return self.__a

    def plot(self, *args, **kwargs ):
        """Convenience method for plotting cross sections. This is the same as
        NCrystal.plot.plot_xsect(material=self,*args,**kwargs), so refer to that
        function for information about allowed arguments."""
        from .plot import plot_xsect
        return plot_xsect( self, *args, **kwargs )

    def dump(self, verbose=0 ):
        """
        Convenience method for print information about contained objects. Use
        verbose argument to set verbosity level to 0 (minimal), 1 (middle), 2
        (most verbose).
        """
        from ._common import print
        _flush()
        print(self.dump_str(verbose=verbose),end='')
        _flush()

    def dump_str(self, verbose=0 ):
        """
        The string-returning sibling of .dump(..). Returns information about
        contained objects as a multi-line string. Use verbose argument to set
        verbosity level to 0 (minimal), 1 (middle), 2 (most verbose).
        """
        res=''
        has_any=False
        for name,descr in [ ('info','Material info'),
                            ('scatter','Scattering process (objects tree)'),
                            ('absorption','Absorption process (object tree)') ]:
            o = getattr(self,name)
            if o:
                has_any=True
                res += '\n>>> '+descr+':\n'
                res += '\n'
                res += o.dump_str(**(dict(verbose=verbose) if name=='info' else {}))
        if not has_any:
            res += '<empty>'
            res += '\n'
        return res

    def xsect( self, *args, **kwargs ):
        """Convenience function which adds up the cross sections from any loaded
        absorption and scatter processes. Refer to the Process.xsect method for
        arguments."""
        procs = [p for p in (self.scatter,self.absorption) if p is not None]
        if not procs:
            raise NCCalcError('.xsect(..) can only be called on'
                              ' LoadedMaterial which contains processes.')
        assert len(procs) in (1,2)
        xs = procs[0].xsect(*args,**kwargs)
        if len(procs)==2:
            xs += procs[1].xsect(*args,**kwargs)
        return xs

    def macroscopic_xsect( self, *args, **kwargs ):
        """Convenience function which calculates cross sections from any loaded
        absorption and scatter processes, and converts to macroscopic cross
        sections via the numerical density of the loaded Info object. Returned
        unit is 1/cm.
        """
        if not self.info or not (self.scatter or self.absorption):
            raise NCCalcError('.macroscopic_xsect(..) can only be called on'
                              ' LoadedMaterial which contains both processes'
                              ' and material Info.')
        return self.info.factor_macroscopic_xs * self.xsect( *args,**kwargs )

    def __str__(self):
        def fmt( x ):
            return str(x) if x else 'n/a'
        return 'LoadedMaterial(Info=%s, Scatter=%s, Absorption=%s)'%( fmt(self.__i),
                                                                      fmt(self.__s),
                                                                      fmt(self.__a) )
    def __repr__(self):
        return str(self)

def directLoad( data, cfg_params='', *, dtype='',
                doInfo = True, doScatter = True, doAbsorption = True ):
    """Convenience function which creates Info, Scatter, and Absorption objects
       directly from a text data source (like a data string or file path) rather
       requiring an on-disk or in-memory file. Such usage obviously precludes
       proper caching behind the scenes, and is intended for scenarios where the
       same data should not be used repeatedly.
    """
    if not dtype and hasattr(data,'dataSourceName') and hasattr(data,'rawData') and hasattr(data,'dataType'):
        #TextData might carry the dtype:
        dtype = data.dataType

    from .misc import AnyTextData
    any_data = AnyTextData( data )
    content = any_data.content

    if not dtype:
        if content.startswith('NCMAT'):
            dtype = 'ncmat'
        else:
            #Try to give a nice error for a common user-error:
            if content[0:1024].lstrip().startswith('NCMAT'):
                raise NCBadInput('NCMAT data must have "NCMAT" as the first 5 characters (must not be preceded by whitespace)')

    if not dtype and any_data.name and '.' in any_data.name:
        p = any_data.name.split('.')
        if len(p)>=2 and p[-1] and p[-1].isalpha() and p[-1] not in ('gz','tgz','bz2','zip','tar'):
            dtype = dtype

    rawi,raws,rawa = _rawfct['multicreate_direct'](content,dtype,cfg_params,doInfo,doScatter,doAbsorption)
    info = Info( ('_rawobj_',rawi) ) if rawi else None
    scatter = Scatter( ('_rawobj_',raws) ) if raws else None
    absorption = Absorption( ('_rawobj_',rawa) ) if rawa else None
    return LoadedMaterial.fromExistingObjects(info,scatter,absorption)

MultiCreated = LoadedMaterial#backwards compat alias
directMultiCreate = directLoad#backwards compat alias

class TextData:
    """Text data accessible line by line, with associated meta-data. This always
       include a UID (useful for comparison and downstream caching purposes) and
       the data type (e.g. "ncmat"). Optionally available is the last known
       on-disk path to a file with the same content, which might be useful in
       case the data needs to be passed to 3rd party software which can only
       work with physical files.

       Text data objects are easily line-iterable, easily providing lines
       (without newline characters): for( auto& line : mytextdata ) {...}.  Of
       course, the raw underlying data buffer can also be accessed if needed.

       The raw data must be ASCII or UTF-8 text, with line endings \\n=CR=0x0A
       (Unix) or \\r\\n=LF+CR=0x0D0A (Windows/DOS). Other encodings might work
       only if 0x00, 0x0A, 0x0D bytes do not occur in them outside of line
       endings.

       Notice that ancient pre-OSX Mac line-endings \\r=LF=0x0D are not
       supported, and iterators will actually throw an error upon encountering
       them. This is done on purpose, since files with \\r on unix might hide
       content when inspected in a terminal can be either confusing, a potential
       security issue, or both.
    """

    def __init__(self,name):
        """create TextData object based on string (same as using createTextData(name))"""
        ll=_rawfct['nc_gettextdata'](name)
        assert len(ll)==5
        self.__rd = ll[0]
        self.__uid = int(ll[1])
        self.__dsn = ll[2]
        self.__datatype= ll[3]
        import pathlib
        self.__rp = pathlib.Path(ll[4]) if ll[4] else None

    @property
    def uid(self):
        """Unique identifier. Objects only have identical UID if all contents and
           metadata are identical."""
        return self.__uid

    @property
    def dataType(self):
        """Data type ("ncmat", "lau", ...)."""
        return self.__datatype

    @property
    def dataSourceName(self):
        """Data source name. This might for instance be a filename."""
        return self.__dsn

    @property
    def rawData(self):
        """Raw access to underlying data."""
        return self.__rd

    @property
    def lastKnownOnDiskLocation(self):
        """Last known on-disk location (returns None if unavailable). Note that there
           is no guarantee against the file having been removed or modified since the
           TextData object was created.
        """
        return self.__rp

    def load( self ):
        """Convenience method which loads the TextData with the directLoad(..)
           function and returns the resulting LoadedMaterial object.
        """
        return LoadedMaterial(self)

    def dump( self ):
        """Convenience method which prints the contents.
        """
        from ._common import print
        _flush()
        print( self.rawData )
        _flush()

    def dump_str( self ):
        """Alias which simply returns .rawData (this method exists exclusively
        for API consistency).
        """
        return self.rawData

    def __str__(self):
        return 'TextData(%s, uid=%i, %i chars)'%(self.__dsn,self.__uid,len(self.__rd))

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        """Line-iteration, yielding lines without terminating newline characters"""
        from io import StringIO
        def chomp(x):
            return x[:-2] if x.endswith('\r\n') else (x[:-1] if x.endswith('\n') else x)
        for e in StringIO(self.__rd):
            yield chomp(e)

def createTextData(name):
    """creates TextData objects based on requested name"""
    return TextData(name)

def setDefaultRandomGenerator(rg):
    """Set the default random generator.

    Note that this can only change the random generator for those processes not
    already created.

    To ensure Python does not clean up the passed function object prematurely,
    the NCrystal python code will keep a reference to the passed function
    eternally (or rather, until the Python process shuts down).

    """
    _rawfct['ncrystal_setrandgen'](rg)

def clearCaches():
    """Clear various caches"""
    _rawfct['ncrystal_clear_caches']()

def _flush():
    import sys
    sys.stdout.flush()
    sys.stderr.flush()

def enableFactoryThreads( nthreads = 'auto' ):
    """Enable threading during object initialisation phase. Supply 'auto' or a
    value >= 9999 to simply use a number of threads appropriate for the system.

    The requested value is the TOTAL number of threads utilised INCLUDING the
    main user thread. Thus, a value of 0 or 1 number will disable this thread
    pool, while for instance calling enableFactoryThreads(8) will result in 7
    secondary worker threads being allocated.

    """
    nt = 9999 if nthreads=='auto' else min(9999,max(1,int(nthreads)))
    _rawfct['ncrystal_enable_factory_threadpool'](nt)
