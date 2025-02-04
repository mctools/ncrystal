
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

Experimental feature for more object-oriented and less efficient access to hkl
lists.

See also https://github.com/mctools/ncrystal/issues/164.

"""

class HKLEntry:

    """A group or family of HKL planes, all sharing the same value of d-spacing
    and structure factor (fsquared). If .is_symequiv evaluates to True, these
    exactly represent a group of symmetry-equivalent planes.
    """
    def __str__( self ):
        hkl=self.hkl_label
        return ( f'HKL( hkl_label=({hkl[0]},{hkl[1]},{hkl[2]}), '
                 f'd={self.d:g}Aa, F2={self.fsquared:g}barn, '
                 f'N={self.mult} )' )

    def __init__(self,hh,kk,ll,mult,dsp,fsq,hklinfotype,issymeqv):
        """For internal usage only, do not create therse objects manually."""
        self.__h = hh
        self.__k = kk
        self.__l = ll
        self.__mult = mult
        self.__dsp = dsp
        self.__fsq = fsq
        self.__hklinfotype = hklinfotype
        self.__issymeqv = issymeqv

    @property
    def hkl_type( self ):
        """The type of HKL group represented by this entry. Returns
        HKLInfoType.SymEqvGroup if this is a group of symmetry-equivalent
        planes.
        """
        return self.__hklinfotype

    @property
    def is_symequiv(self):
        """Returns True if .hkl_type equals HKLInfoType.SymEqvGroup."""
        return self.__issymeqv

    @property
    def hkl_label( self ):
        """Returns a hkl label for the entry, i.e. one of the hkl points in the
        group as a tuple of three integers: (h,k,l)."""
        return (int(self.__h[0]),int(self.__k[0]),int(self.__l[0]))

    @property
    def h( self ):
        """
        An array of h values. Note that this has half the length of
        .multiplicity, since we exclude entries that can be generated from each
        other by a mere sign flip.
        """
        return self.__h

    @property
    def k( self ):
        """
        An array of k values. Note that this has half the length of
        .multiplicity, since we exclude entries that can be generated from each
        other by a mere sign flip.
        """
        return self.__k

    @property
    def l( self ): # noqa E743
        """
        An array of l values. Note that this has half the length of
        .multiplicity, since we exclude entries that can be generated from each
        other by a mere sign flip.
        """
        return self.__l

    @property
    def dspacing( self ):
        """The d-spacing value in units of angstrom."""
        return self.__dsp

    @property
    def d( self ):
        """The d-spacing value in units of angstrom."""
        return self.__dsp

    @property
    def fsquared( self ):
        """The squared structure factor (F^2) in units of barn."""
        return self.__fsq

    @property
    def f2( self ):
        """The squared structure factor (F^2) in units of barn."""
        return self.__fsq

    @property
    def multiplicity( self ):
        """The number of hkl points in the group."""
        return self.__mult

    @property
    def mult( self ):
        """The number of hkl points in the group."""
        return self.__mult

def _iter_hklobjects( info ):
    hklinfotype = info.hklInfoType()
    issymeqv = info.hklIsSymEqvGroup()
    for hh,kk,ll,mult,dsp,fsq in info.hklList(all_indices=True):
        yield HKLEntry(hh,kk,ll,mult,dsp,fsq,
                       hklinfotype,issymeqv)

