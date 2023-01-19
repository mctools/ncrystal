"""Module with internal utilities used by several NCrystal modules"""

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2022 NCrystal developers                                   ##
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

import warnings
import math
try:
    import numpy as _np
except ImportError:
    _np = None
import importlib
__name = 'NCrystal._common' if __name__=='__main__' else __name__
NCrystalUserWarning = importlib.import_module( '..', __name ).NCrystalUserWarning

#Common print function for NCrystal modules (allowing one to capture output of just NCrystal's python layer):
_stdprint = [print]

def print(*args,**kwargs):
    _stdprint[0](*args,**kwargs)

def set_ncrystal_print_fct( fct ):
    _stdprint[0] = fct
    return print


def warn(msg):
    """Emit NCrystalUserWarning via standard warnings.warn function"""
    warnings.warn( NCrystalUserWarning( str(msg) ))

class WarningSpy:
    """Context manager which spies on any warnings emitted via warnings
    module and returns list of messages and categories.

    It does so by temporarily intercepting the warnings.showwarning
    function.

    Usage example:
        with WarningSpy() as warnlist:
            fct()
        print( "Warnings emitted:", warnlist )
    """

    def __init__(self,block=False):
        self.__block = bool(block)

    def __enter__(self):
        self.__l = []
        self.__orig = warnings.showwarning
        warnings.showwarning = self.__spy
        return self.__l

    def __exit__(self,*args,**kwargs):
        warnings.showwarning = self.__orig

    @staticmethod
    def _fmtwarning( message, category ):
        catname = None
        if category is not None:
            catname = str(getattr(category,'__name__',category)).strip()
        elif isinstance(message,Warning):
            catname = str(message.__class__.__name__).strip()
        if not catname:
            catname='Warning'
        return ( ( str(message).strip() or '<no message>'), catname )

    def __spy( self, message, category, *args, **kwargs ):
        #only store string objects, to decouple lifetimes of more complex objects:
        self.__l.append( WarningSpy._fmtwarning( message, category ) )
        if not self.__block:
            return self.__orig( message, category, *args, **kwargs )

_fracdb=[None]
def prettyFmtValue(x):
    """Recognises common fractions in truncated values like '0.33333333' or
    '0.4166667' and formats them as '1/3' and '5/12'. This not only saves space,
    but can also reinject accuracy in atom positions (ncrystal_verifyatompos can
    be used to check this, but it seems to be a clear improvement in all cases
    tested so far).
    """
    global _fracdb
    assert 0.0 <= x < 1.0
    if x==0.0:
        return '0'
    if x==0.5:
        return '1/2'
    def find_nearest_idx(arr,val):
        if _np:
            i = _np.searchsorted(arr, val)
        else:
            if val <= arr[0]:
                i = 0
            elif val >= arr[-1]:
                i = len(arr)-1
            else:
                for i,e in enumerate(arr):
                    if e>=val:
                        break
        if i == 0:
            return 0
        if i == len(arr):
            return len(arr)-1
        return i if abs(val-arr[i])<abs(val-arr[i-1]) else i-1
    def initFractions():
        fractions=set()
        nmax=40
        for a in range(1,nmax):
            for b in range(a+1,nmax+1):
                g=math.gcd(a,b)
                fractions.add( (a//g, b//g ) )
        fractions=list(sorted( (a/b,a,b) for a,b in fractions))
        values = list(e for e,_,_ in fractions)
        if _np:
            values = _np.asarray( values, dtype=float )
        fmt = [ f'{a}/{b}' for _,a,b in fractions ]
        return values,fmt
    xfmt = '%.14g'%x
    if xfmt.startswith('0.'):
        xfmt = xfmt[1:]
    if _fracdb[0] is None:
        _fracdb[0] = initFractions()
    dbvals, dbfmt = _fracdb[0]
    i = find_nearest_idx(dbvals,x)
    v=dbvals[i]
    #check if credible that the numbers are identical (must be close AND the str
    #representations must be consistent with truncation):
    expected_at_prec = ((f'%.{len(xfmt)-1}g')%v)[1:]
    if abs(v-x)<1e-7 and expected_at_prec == xfmt:
        return dbfmt[i]
    if abs(v-x)<1e-14:
        return dbfmt[i]
    return xfmt
