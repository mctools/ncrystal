
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

Utilities for filters and windows, where a neutron beam is simply attenuated
by a material.

"""

__all__ = [ 'NCrystalFilter' ]

class NCrystalFilter:

    """Macroscopic total cross section (scattering plus absorption) of an
    isotropic material, for attenuating a neutron beam passing through it
    (e.g. a filter or a window). A neutron travelling a distance L (in cm)
    through the material is transmitted with probability exp(-xsect*L).

    The cross sections are evaluated from a piecewise linear table vs. the
    neutron wavelength, which reproduces the cross section of the material
    within a relative tolerance of 1e-3 (relative to max(xs,1e-12 barn per
    atom), so cross sections vanishing at wavelength 0 can be tabulated).
    Beyond the last point of the table, the last segment is extrapolated
    linearly. An exception is raised if the cross section can not be tabulated
    reliably.

    The options are reserved for future use, and must be None or empty.
    """

    def __init__( self, cfgstr, options = None ):
        from ._chooks import _get_raw_cfcts
        self.__cfgstr = str(cfgstr)
        self.__options = options or None
        wl, macroxs = _get_raw_cfcts()['filtertable']( self.__cfgstr,
                                                       self.__options )
        wl.setflags( write = False )
        macroxs.setflags( write = False )
        self.__wl, self.__macroxs = wl, macroxs

    @property
    def cfgstr( self ):
        """The cfg-string of the material."""
        return self.__cfgstr

    @property
    def options( self ):
        """The options (None if not provided)."""
        return self.__options

    @property
    def table( self ):
        """The table, as two read-only numpy arrays: the wavelengths (in Aa)
        and the macroscopic cross sections (in 1/cm). The first point is the
        limit for wavelength -> 0. Discontinuities (e.g. Bragg edges) are
        represented by two points with the same wavelength, where the first has
        the value just below the discontinuity, and the second the value just
        above it."""
        return self.__wl, self.__macroxs

    def xsect( self, ekin = None, wl = None ):
        """Macroscopic cross section (in 1/cm) at the given neutron kinetic
        energy (in eV) or wavelength (in Aa). These can be numbers (returning a
        float), or arrays (returning a numpy array)."""
        from ._numpy import _ensure_numpy, _np
        from .constants import ekin2wl
        from .exceptions import NCBadInput
        if ( ekin is None ) == ( wl is None ):
            raise NCBadInput('Please provide exactly one of the "ekin" or'
                             ' "wl" parameters.')
        _ensure_numpy()
        w = _np.asarray( ekin2wl( ekin ) if wl is None else wl, dtype = float )
        x, y = self.__wl, self.__macroxs
        #Linear interpolation in the segment with x[i] <= w < x[i+1] (at a
        #discontinuity, this is the segment above it):
        wi = _np.clip( w, x[0], x[-1] )
        i = _np.clip( _np.searchsorted( x, wi, side = 'right' ) - 1,
                      0, len(x) - 2 )
        x0, x1, y0, y1 = x[i], x[i+1], y[i], y[i+1]
        dx = _np.where( x1 > x0, x1 - x0, 1.0 )
        res = _np.where( x1 > x0, y0 + ( ( wi - x0 ) / dx ) * ( y1 - y0 ), y1 )
        #Linear extrapolation of the last segment beyond the table (clamped at
        #0):
        dxl = x[-1] - x[-2]
        slope = ( y[-1] - y[-2] ) / dxl if dxl > 0.0 else 0.0
        if slope != 0.0:
            extrap = _np.maximum( y[-1] + ( w - x[-1] ) * slope, 0.0 )
        else:
            extrap = _np.full_like( w, y[-1] )
        res = _np.where( w >= x[-1], extrap, res )
        res = _np.where( w <= x[0], y[0], res )
        return float( res ) if res.ndim == 0 else res
