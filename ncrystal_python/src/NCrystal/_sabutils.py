
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


# Various internal utilities related to SAB arrays.

def check_sab_shape( *, nalpha, nbeta, sab ):
    if not ( ( sab.ndim == 1 and len(sab) == nalpha*nbeta )
             or ( sab.ndim == 2 and  sab.shape == ( nbeta, nalpha ) ) ):
        from .exceptions import NCBadInput
        raise NCBadInput("sab must be 1D array of length nalpha*nbeta"
                         " or a 2D array of shape (nbeta,nalpha).")

def reshape_sab_2D( *, nalpha, nbeta, sab ):
    check_sab_shape( nalpha = nalpha, nbeta = nbeta, sab = sab )
    return sab.reshape( ( nbeta, nalpha) ) if sab.ndim == 1 else sab

def reshape_sab_1D( *, nalpha, nbeta, sab ):
    check_sab_shape( nalpha = nalpha, nbeta = nbeta, sab = sab )
    return sab.reshape( ( nbeta*nalpha, ) ) if sab.ndim == 2 else sab

def trim_knl_edges( *, alphagrid, betagrid, sab, view = False ):
    """Trim edges of a S(alpha,beta) grid where S=0 for a whole grid row
       or column. Specifically, both upper and lower beta edges are
       trimmer, as well as the upper alpha edge. The lower alpha edge is
       always left alone.

       Whether or not such trimming is sensible depends on the context and is
       the responsibility of the calling code to decide.

       Return values are the trimmed (alphagrid, betagrid, sab). If view = True,
       the returned arrays are just a view into the originals if possible,
       otherwise they will be created as a new copy.

       The sab array must be shaped as a 2D (shape=(nbeta,nalpha)) or 1D(
       shape=(nbeta*nalpha,)) array, and the returned sab array will be encoded
       in the same scheme as the input.
    """
    sab_ndim = sab.ndim
    sab = reshape_sab_2D( nalpha = len(alphagrid),
                          nbeta = len(betagrid),
                          sab = sab )
    b, a, s = _trim_edges_2d( betagrid, alphagrid, sab, keep_low_y = True )
    assert s.ndim == 2
    if sab_ndim == 1:
        s = reshape_sab_1D( nalpha = len(a),
                            nbeta = len(b),
                            sab = sab )
    return ( a, b, s ) if view else ( a.copy(), b.copy(), s.copy() )

def _trim_edges_2d( x, y, s, keep_low_y = False ):
    # Assume a function s(x,y) is described by grid values x, y and function
    # values at the corresponding grid points s (so s[i,j] is the value on the
    # grid point (x[j],y[i])). This function then trims the grid to discards any
    # edge-columns or edge-rows where s[i,j] is 0 everywhere.
    #
    #x columns: axis=0
    #y rows: axis=1
    from ._numpy import _np, _ensure_numpy
    _ensure_numpy()

    assert x.ndim==1
    assert y.ndim==1
    assert s.ndim==2
    nx, ny, ns = len(x), len(y), ( s.shape[0] * s.shape[1] )
    assert nx == s.shape[0]
    assert ny == s.shape[1]
    assert ns == nx * ny
    keep_x = _np.any( s != 0, axis=1 )
    keep_y = _np.any( s != 0, axis=0 )
    assert len(keep_x) == nx
    assert len(keep_y) == ny
    if keep_low_y:
        keep_y[0] = True

    def firstlast( v ):
        return ( _np.argmax(v), len(v)-1-_np.argmax(v[::-1]) )
    first_y, last_y = firstlast(keep_y)
    first_x, last_x = firstlast(keep_x)
    return ( x[first_x:last_x + 1],
             y[first_y:last_y + 1],
             s[first_x:last_x + 1,first_y:last_y + 1] )
