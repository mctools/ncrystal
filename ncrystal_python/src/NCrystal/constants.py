
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

Various constants and neutron energy / wavelength / wavenumber conversions
methods like for instance wl2ekin and ekin2wl. Note that all conversion
functions can be used with numpy arrays in addition to scalar numbers.

"""

#some constants (NB: Copied here from NCMath.hh - must keep synchronized!!):

from ._numpy import _ensure_numpy, _np

constant_c  = 299792458e10#  speed of light in Aa/s
constant_dalton2kg =  1.660539040e-27#  amu to kg
constant_dalton2eVc2 =  931494095.17#  amu to eV/c^2
constant_avogadro = 6.022140857e23#  mol^-1
constant_boltzmann = 8.6173303e-5#  eV/K
const_neutron_mass_amu = 1.00866491588#  [amu]
constant_planck = 4.135667662e-15 # [eV*s]

kPi        = 3.1415926535897932384626433832795028841971694
k2Pi       = 6.2831853071795864769252867665590057683943388
k4Pidiv100 = 0.125663706143591729538505735331180115367886776
k4PiSq     = 39.4784176043574344753379639995046045412547976
kInf       = float('inf')

_as_np_array = ( lambda x : _np.asarray(x,dtype=float) ) if _np else ( lambda *a,**kwargs : _ensure_numpy() )
def wlsq2ekin( wlsq ):
    """Neutron wavelength squared (angstrom^2) to energy (eV)"""
    if hasattr(wlsq,'__len__'):
        x = _as_np_array( wlsq )
        wlsqnonzero = ( x != 0.0 )
        wlsqinv = 1.0 / _np.where( wlsqnonzero, x, 1.0)#fallback 1.0 wont be used
        return _const_wlsqekin * _np.where( wlsqnonzero, wlsqinv, kInf )
    else:
        return _const_wlsqekin / wlsq if wlsq else kInf

def wl2ekin( wl ):
    """Neutron wavelength (angstrom) to energy (eV)"""
    if hasattr(wl,'__len__'):
        x = _as_np_array( wl )
        return wlsq2ekin( x*x )
    else:
        return wlsq2ekin( wl*wl )

def ekin2wl( ekin ):
    """Neutron energy (eV) to wavelength (angstrom)"""
    if hasattr(ekin,'__len__'):
        x = _as_np_array( ekin )
        ekinnonzero = x != 0.0
        ekininv = 1.0 / _np.where( ekinnonzero, x, 1.0)#fallback 1.0 wont be used
        return _np.sqrt( _const_wlsqekin * _np.where( ekinnonzero, ekininv, kInf ) )
    else:
        from math import sqrt as _math_sqrt
        return _math_sqrt( _const_wlsqekin / ekin ) if ekin else kInf

def ekin2wlsq( ekin ):
    """Neutron energy (eV) to wavelength squared (angstrom^2)"""
    if hasattr(ekin,'__len__'):
        x = _as_np_array( ekin )
        ekinnonzero = x != 0.0
        ekininv = 1.0 / _np.where( ekinnonzero, x, 1.0)#fallback 1.0 wont be used
        return _const_wlsqekin * _np.where( ekinnonzero, ekininv, kInf )
    else:
        return ( _const_wlsqekin / ekin ) if ekin else kInf

def ekin2wlsqinv( ekin ):
    """Neutron energy (eV) to inverse wavelength squared (1/angstrom^2)"""
    return ekin * _const_inv_wlsqekin#constant is 1/_const_wlsqekin

_const_wlsqekin     = 0.081804209605330899    # ekin = _const_wlsqekin /wl^2
_const_inv_wlsqekin = 12.22430978582345950656 # 1 / _const_wlsqekin
_const_ekin2ksq_factor = k4PiSq * _const_inv_wlsqekin
_const_ksq2ekin_factor = 1.0 / _const_ekin2ksq_factor

def ekin2ksq( ekin ):
    """Neutron energy (eV) to wavenumber squared, (k^2, in units of 1/angstrom^2)"""
    return _const_ekin2ksq_factor * ekin

def ekin2k( ekin ):
    """Neutron energy (eV) to wavenumber, (k, in units of 1/angstrom)"""
    if hasattr(ekin,'__len__'):
        _ = _const_ekin2ksq_factor * _as_np_array( ekin )
        return _np.sqrt( _ )
    else:
        from math import sqrt
        return sqrt( _const_ekin2ksq_factor * ekin )

def ksq2ekin( ksq ):
    """Wavenumber squared, (k^2, in units of 1/angstrom^2) to neutron energy (eV)"""
    return _const_ksq2ekin_factor * ksq

def wl2k( wl ):
    """Neutron wavelength (angstrom) to wavenumber, (k, in units of 1/angstrom)"""
    if hasattr(wl,'__len__'):
        x = _as_np_array( wl )
        wlnonzero = x != 0.0
        wlinv = 1.0 / _np.where( wlnonzero, x, 1.0)#fallback 1.0 wont be used
        return k2Pi * _np.where( wlnonzero, wlinv, kInf )
    else:
        return k2Pi / wl if wl else kInf

def wl2ksq( wl ):
    """Neutron wavelength (angstrom) to wavenumber squared, (k^2, in units of 1/angstrom^2)"""
    return ( wl2k(wl) )**2

def k2wl( k ):
    """Wavenumber, (k, in units of 1/angstrom) to neutron wavelength (angstrom)"""
    return wl2k( k )#using that k2wl = wl2k (both are f(x)=2pi/x)

def k2ekin( k ):
    """Wavenumber, (k, in units of 1/angstrom) to neutron energy (eV)"""
    return ksq2ekin( k * k )
