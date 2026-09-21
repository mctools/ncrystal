#!/usr/bin/env python3

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

# NEEDS: numpy

# Verify that arrays passed to the C library (energies, wavelengths, VDOS
# spectra) are handled correctly, no matter their memory layout, dtype or
# container type.

import NCrystalDev as NC
import numpy as np


def main():

    mat = NC.load('Al_sg225.ncmat;dcutoff=0.5;temp=300')
    scat = mat.scatter
    absn = mat.absorption

    ekin = np.geomspace(1e-4, 1.0, 40)
    assert ekin.flags['C_CONTIGUOUS'] and ekin.dtype == np.float64

    #Reference results for the plain contiguous float64 array. Also checked against
    #scalar calls, to make sure the reference itself is sound:
    ref_scat = scat.xsect(ekin)
    ref_abs = absn.xsect(ekin)
    assert len(ref_scat) == len(ekin)
    for i, e in enumerate(ekin):
        assert ref_scat[i] == scat.xsect(float(e))
        assert ref_abs[i] == absn.xsect(float(e))
    assert len(set(ref_scat)) > 20#values must actually vary for a test to be useful

    def layout_variants( values ):
        """Returns a dict of arrays or other containers, all with the same values
        but in different memory layouts."""
        n = len(values)
        v = {}
        #column of a 2D array (strided):
        a = np.full( (n,3), -1.0 )
        a[:,1] = values
        v['column of 2D array'] = a[:,1]
        #every second element (strided):
        a = np.full( 2*n, -1.0 )
        a[::2] = values
        v['every second element'] = a[::2]
        #every third, from an offset:
        a = np.full( 3*n+2, -1.0 )
        a[1::3][:n] = values
        v['every third element with offset'] = a[1::3][:n]
        #negative stride (reversed view of reversed data):
        v['negative stride'] = values[::-1].copy()[::-1]
        #contiguous, but read-only:
        a = values.copy()
        a.flags.writeable = False
        v['read-only'] = a
        #Fortran ordered 2D array with a single row is contiguous in memory but
        #two-dimensional, so it is tested for the error below and not here.
        #other dtypes/containers (values are exactly representable in each):
        v['list'] = [float(e) for e in values]
        v['tuple'] = tuple(float(e) for e in values)
        v['float64 non-native byte order'] = values.astype('>f8')
        return v

    def check_variants( name, fct, ref, values, exact_dtype_conv = None ):
        for lbl, arr in layout_variants( values ).items():
            res = fct( arr )
            assert isinstance( res, np.ndarray )
            assert res.shape == ref.shape
            if not np.array_equal( res, ref ):
                print(f'FAILURE: {name} gave wrong result for input: {lbl}')
                raise SystemExit(1)
        print(f'{name}: all memory layouts and containers ok')

    check_variants( 'scatter xsect(ekin)', lambda a : scat.xsect(a), ref_scat, ekin )
    check_variants( 'absorption xsect(ekin)', lambda a : absn.xsect(a), ref_abs, ekin )
    check_variants( 'scatter xsect(ekin=)', lambda a : scat.xsect(ekin=a), ref_scat, ekin )

    #Wavelength interface (the wavelengths are converted to energies internally):
    wl = np.linspace(0.5, 12.0, 40)
    ref_wl = scat.xsect(wl=wl)
    assert len(set(ref_wl)) > 20
    check_variants( 'scatter xsect(wl=)', lambda a : scat.xsect(wl=a), ref_wl, wl )

    #Types which are not float64 must be converted rather than reinterpreted:
    e32 = ekin.astype(np.float32)
    assert np.array_equal( scat.xsect(e32), scat.xsect(e32.astype(np.float64)) )
    eint = np.array([1, 2, 3, 5, 8], dtype=np.int64)
    assert np.array_equal( scat.xsect(eint), scat.xsect(eint.astype(np.float64)) )
    eint32 = eint.astype(np.int32)
    assert np.array_equal( scat.xsect(eint32), scat.xsect(eint.astype(np.float64)) )
    print('float32 and integer arrays are converted correctly')

    #The repeat parameter with a strided array:
    a = np.full( (len(ekin),2), -1.0 )
    a[:,0] = ekin
    res = scat.xsect( a[:,0], repeat=3 )
    assert np.array_equal( res, scat.xsect( ekin, repeat=3 ) )
    assert len(res) == 3*len(ekin)
    print('repeat with strided array ok')

    #The input array must never be modified:
    a = np.full( (len(ekin),2), -1.0 )
    a[:,0] = ekin
    a_orig = a.copy()
    scat.xsect( a[:,0] )
    assert np.array_equal( a, a_orig )
    print('input arrays are not modified')

    #Multi-dimensional arrays are not supported, and should give a clear error
    #instead of being silently misinterpreted:
    for bad in ( np.ones((4,3)), np.ones((1,4)), np.ones((4,1)) ):
        try:
            scat.xsect( bad )
        except NC.NCBadInput as e:
            assert 'one-dimensional' in str(e), str(e)
        else:
            raise SystemExit('FAILURE: Multi-dimensional array was accepted')
    print('multi-dimensional arrays are rejected with a clear error')

    #Sampling functions use the same conversions of their input arrays:
    for arr in layout_variants(np.full(50, 0.0253)).values():
        ef, mu = scat.sampleScatterIsotropic(arr)
        assert len(ef) == 50 and len(mu) == 50
        assert np.all( np.abs(mu) <= 1.0 ) and np.all( ef > 0.0 )
    print('sampling functions ok')

if __name__=='__main__':
    main()
