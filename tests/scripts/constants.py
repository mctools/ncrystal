#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

# NEEDS: numpy

import NCTestUtils.enable_fpe
import NCrystal.constants as ncc
import math
import numpy as np

def get_all_vals():
    for k in dir(ncc):
        v = getattr(ncc,k) if not k.startswith('_') else None
        if isinstance(v,float):
            yield k,v
for k,v in sorted(get_all_vals()):
    print( f"NCrystal.constants.{k} = {repr(v)}" )

def reldiff( x, y ):
    if math.isinf(x):
        return float('inf') if ( not math.isinf(y) or ( ( x>0 ) != ( y>0 ) ) ) else 0.0
    return abs(x-y)/(max(1e-300,abs(x)+abs(y)))

def require_flteq( x, y ):
    okfct = lambda a,b : bool( reldiff( a, b ) < tol )
    tol = 1e-13
    if hasattr( x, '__len__' ):
        if not len(x) == len(y) or any( ( not okfct(a,b) ) for a,b in zip(x,y) ):
            raise RuntimeError(f'numpy flteq failed for arrays x={x} and y={y}!')
    elif not okfct(x,y):
        raise RuntimeError(f'require_flteq( x={x}, y={y} ) failed!')

def require_flteq_inv( x, y ):
    if hasattr( x, '__len__' ):
        for a,b in zip(x,y):
            require_flteq_inv( a,b)
        return
    if not x:
        require_flteq( x, 1/y )
    else:
        require_flteq( 1/x, y )

testwls = [ 0.0, 1e-10, 1.8, 10, 1e4, float('inf') ]
testwls.append( np.asarray(testwls,dtype=float))
for wl in testwls:
    is_array_test = hasattr(wl,'__len__')
    ekin = ncc.wl2ekin( wl )
    if is_array_test:
        print( f" ---- test array of wavelengths (printouts disabled)")
    else:
        print( f" ---- test wl = {wl} Aa")

    if not is_array_test:
        print( f"      -> wl2ekin({wl:g} Aa) = {ekin:g} eV")

    require_flteq( wl, ncc.ekin2wl(ekin) )
    if not is_array_test:
        print(  '      -> checked ekin2wl consistency')

    wlsq = ncc.ekin2wlsq(ekin)
    if not is_array_test:
        print( f"      -> ekin2wlsq({ekin:g} eV) = {wlsq:g} Aa^2")
    require_flteq( wl, np.sqrt( wlsq ) )
    if not is_array_test:
        print(  '      -> checked ekin2wlsq consistency')

    wlsqinv = ncc.ekin2wlsqinv(ekin)
    if not is_array_test:
        print( f"      -> ekin2wlsqinv({ekin:g} eV) = {wlsqinv:g} Aa^-2")
    require_flteq_inv( wl, np.sqrt(wlsqinv) )
    if not is_array_test:
        print(  '      -> checked ekin2wlsqinv consistency')

    _ekinalt = ncc.wlsq2ekin(wl*wl)
    if not is_array_test:
        print( f"      -> wlsq2ekin({wl*wl:g} Aa^2) = {_ekinalt:g} eV")
    require_flteq( ekin, _ekinalt )
    if not is_array_test:
        print(  '      -> wlsq2ekin consistency')

    k = ncc.wl2k( wl )
    if not is_array_test:
        print( f"      -> wl2k({wl:g} Aa) = {k:g} Aa^-1")
    require_flteq_inv( (1.0/(2*math.pi))*wl, k )

    _ksq = ncc.wl2ksq( wl )
    if not is_array_test:
        print( f"      -> wl2ksq({wl:g} Aa) = {_ksq:g} Aa^-1")
    require_flteq( k*k,  _ksq )

    _ekin = ncc.ksq2ekin( k*k )
    if not is_array_test:
        print( f"      -> ksq2ekin({k*k:g} Aa) = {_ekin:g} eV")
    require_flteq( ekin, _ekin )

    _k = ncc.ekin2k( ekin )
    if not is_array_test:
        print( f"      -> ekin2k({ekin:g} eV) = {_k:g} Aa^-1")
    require_flteq( k, _k )

    _ekin2 = ncc.k2ekin( k )
    if not is_array_test:
        print( f"      -> k2ekin({k:g} Aa^-1) = {_ekin2:g} eV")
    require_flteq( _ekin2, ekin )

    _wl2 = ncc.k2wl( k )
    if not is_array_test:
        print( f"      -> k2wl({k:g} Aa^-1) = {_wl2:g} Aa")
    require_flteq( _wl2, wl )

    _ksq = ncc.ekin2ksq( ekin )
    if not is_array_test:
        print( f"      -> ekin2ksq({ekin:g} eV) = {_ksq:g} Aa^-2")
    require_flteq( k*k, _ksq )

print("Done")
