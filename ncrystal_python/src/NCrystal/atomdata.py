
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

Utilities for accessing NCrystal's database of atomic data, with information
about atomic masses, scattering lengths, etc. Also contains a few other related
utilities, like a list of all element names.

"""

__atomdb={}
def atomDB(Z,A=None,throwOnErrors=True):
    """Access internal database with data for isotopes and natural elements.

    If A is provided, both A and Z must be integers, thus defining a specific isotope.

    If Z is an integer and A is 0 or None, the corresponding natural element is provided.

    Finally, the function can be called with a string identifying either natural
    elements or isotopes: atomDB("Al"), atomDB("He3"), ...

    In all cases, in case of errors or missing entries in the database, either
    an NCBadInput exception is thrown (throwOnErrors==True) or None is
    returned (when throwOnErrors==False).
    """
    global __atomdb
    import numbers
    if isinstance(Z,numbers.Integral):
        Z=int(Z)
        key=(Z,int(A or 0))
        strkey=False
    else:
        assert A is None,"Do not supply two arguments unless the first argument is an integer"
        assert isinstance(Z,str),"The first argument to the function must either be of int or str type"
        key=Z
        strkey=True
    obj=__atomdb.get(key,None)
    if obj:
        return obj
    from ._chooks import _get_raw_cfcts,_str2cstr
    _rawfct = _get_raw_cfcts()
    if strkey:
        rawatomdata=_rawfct['ncrystal_create_atomdata_fromdbstr'](_str2cstr(key))
    else:
        rawatomdata=_rawfct['ncrystal_create_atomdata_fromdb'](*key)
    if not _rawfct['ncrystal_valid'](rawatomdata):
        if not throwOnErrors:
            return None
        if strkey:
            s='key="%s"'%key
        else:
            if key[1]==0:
                s='Z=%i'%key[0]
            else:
                s='Z=%i,A=%i'%key
        from .exceptions import NCBadInput
        raise NCBadInput('atomDB: Could not find entry for key (%s)'%s)
    from .core import AtomData
    ad = AtomData(rawatomdata)
    assert ad.isElement()
    Z,A = ad.Z(), (ad.A() if ad.isSingleIsotope() else 0)
    keys=[ (Z,A)]
    if Z==1 and A==2:
        keys+=['H2','D']
    elif Z==1 and A==3:
        keys+=['H3','T']
    else:
        assert ad.isNaturalElement() or ad.isSingleIsotope()
        keys += [ ad.description(False) ]#guaranteed to give just symbol for natelem/singleisotope!
    assert key in keys#Should always be true unless we forgot some keys above
    assert ad.description(False) in keys#Should also be true, given guarantees for AtomData::description(false)
    for k in keys:
        __atomdb[k] = ad
    return ad

def iterateAtomDB(objects=True):
    """Iterate over all entries in the internal database with data for isotopes and
       natural elements. If objects=True, AtomData objects are returned. If
       objects=False, (Z,A) values are returned (A=0 indicates a natural
       element)."""
    from ._chooks import _get_raw_cfcts
    _rawfct = _get_raw_cfcts()
    for z,a in _rawfct['atomdb_getall_za']():
        yield atomDB(z,a) if objects else (int(z),int(a))

__all_element_names=tuple(["H", "He", "Li", "Be", "B", "C",
                           "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                           "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
                           "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
                           "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
                           "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                           "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
                           "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                           "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
                           "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
                           "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",
                           "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts",
                           "Og"])

__elem2z = [None]
def elementNameToZValue( element_name, allow_isotopes = False ):
    """Return Z value for element name (e.g. element_name="C" returns
    6). Returns None if argument is not an element name. If
    allow_isotopes==True, isotope names like "B10", "D", etc. will be decoded as
    well and return the Z value of the corresponding element.
    """
    if allow_isotopes:
        if element_name in ('D','T'):
            return 1
        while element_name and element_name[-1].isdigit():
            element_name = element_name[:-1]
    if __elem2z[0] is None:
        __elem2z[0] = dict( ( e, i+1 ) for i,e in enumerate(__all_element_names) )
    return __elem2z[0].get(element_name, None )

def elementZToName( Z ):
    """
    Return element name corresponding to given Z value, or None in case of invalid or out of range argument.
    """
    import numbers
    i = int(Z)-1 if isinstance(Z,numbers.Integral) else -1
    return __all_element_names[i] if 0 <= i < len(__all_element_names) else None

def allElementNames():
    """Like knownElementNames, but also returns names of elements for which
    NCrystal has no data. See also knownElementNames()."""
    return __all_element_names

def isElementName( label ):
    """Determine whether label is a an element name like "H" or "He". See also
    isKnownElement(..).

    """
    return ( ( elementNameToZValue(label) is not None )
             if __elem2z[0] is None else ( label in __elem2z[0] ) )

__all_known_element_names=[None,None]#list,set
def knownElementNames():
    """Returns tuple of all element names for which NCrystal has data
    values. The elements are ordered by increasing Z value. See also
    allElementNames()."""
    if __all_known_element_names[0] is None:
        ll=[]
        for a in iterateAtomDB():
            if a.isNaturalElement():
                ll.append( (a.Z(),a.elementName()) )
        ll = tuple( name for z,name in sorted(ll) )
        __all_known_element_names[0] = ll
        __all_known_element_names[1] = frozenset(ll)
    return __all_known_element_names[0]

def isKnownElement( label ):
    """Determine whether or not label is an element name for which NCrystal has
    data values. See also isElementName(..).
    """
    if __all_known_element_names[1] is None:
        knownElementNames()#trigger init
    return label in __all_known_element_names[1]
