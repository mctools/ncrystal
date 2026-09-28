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

# Test expansion of VDOS curves into scattering kernels. Apart from basic
# usage, this includes tests of non-legacy expansions (vdoslux 2000..2006) at
# low temperatures, where structure on the scale of kT must be resolved:
#
#  * G1 must be sampled finely enough, even when the input VDOS grid is
#    coarse (as for the H VDOS in the acrylic glass stdlib file, where the
#    grid spacing corresponds to ~5.8kT at 5K).
#  * Gn spectra must only be thinned when smooth. This is tested with a VDOS
#    with features a few meV wide (a soft mode), but extending to 0.55eV (an
#    oscillator), whose Gn spectra at 14K extend over thousands of kT while
#    having structure on the scale of kT. Thinning such spectra based on
#    their number of points only, gave cross sections wrong by 50%.
#  * Expansions needing excessive resources must fail with a CalcError.
#
# The resolution is tested via the detailed balance relation,
# Gn(+E)=exp(-E/kT)*Gn(-E), evaluated in the middle of the Gn bins (i.e.
# testing the linear interpolation).

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
import NCrystalDev as NC
import NCrystalDev.exceptions as nc_exceptions
import NCrystalDev.vdos as nc_vdos
import numpy as np
from NCrystalDev.constants import constant_boltzmann
from NCTestUtils.common import ensure_error

def test( vdos, m, T, do_plot, vdoslux, target_emax = None ):
    if do_plot:
        import NCrystalDev.plot as nc_plot
        nc_plot.plot_vdos(vdos)
    print(f"Calling nc_vdos.extractKnl( m={m:g}, T={T:g}, "
          f"vdoslux={vdoslux}, target_emax={target_emax or 0.0:g} )")
    knl = nc_vdos.extractKnl( vdos,
                              mass_amu = m,
                              temperature = T,
                              target_emax = target_emax,
                              vdoslux = vdoslux,
                              plot = do_plot )
    for k,v in knl.items():
        if hasattr(v,'shape'):
            v = 'NumpyArray( shape={} )'.format(*v.shape)
        print(f" Got {k} : {v}")
    print()

def main( do_plot ):
    def t( *a, **kw ):
        test( *a, **kw, do_plot = do_plot )
    vdos = nc_vdos.createVDOSDebye(400.0)
    t( vdos, m = 27.0, T = 0.5, vdoslux = 1 )
    expected_error = ( "VDOS expansion too slow - can not reach E=5000eV"
                       " after 10000 phonon convolutions (likely causes:"
                       " either the target energy value is too high, vdoslux"
                       " too low, the temperature too high, or the VDOS is"
                       " very unusual).")
    t( vdos, m = 27.0, T = 300.0, vdoslux = 1 )
    t( vdos, m = 27.0, T = 300.0, vdoslux = 1, target_emax = 50 )
    t( vdos, m = 27.0, T = 300.0, vdoslux = 4 )
    t( vdos, m = 27.0, T = 300.0, vdoslux = 4, target_emax = 0.5 )
    with ensure_error(nc_exceptions.NCCalcError,expected_error):
        t( vdos, m = 27.0, T = 0.5, vdoslux = 1, target_emax = 5000 )
    test_g1_resolution()
    test_thinning()
    test_resource_limits()

def gn_dbcheck( egrid, gn, temp, relfloor = 0.0 ):
    #Returns binwidth/kT and largest deviation from detailed balance, in the
    #middle of bins with energies in 0.25kT..10kT:
    kt = constant_boltzmann * temp
    emid = 0.5 * ( egrid[1:] + egrid[:-1] )
    emid = emid[ ( emid >= 0.25*kt ) & ( emid <= 10.0*kt ) ]
    gn_up = np.interp( emid, egrid, gn )
    gn_down = np.interp( -emid, egrid, gn )
    ok = gn_down > relfloor * gn.max()
    assert ok.any()
    ratio = gn_up[ok] * np.exp( emid[ok] / kt ) / gn_down[ok]
    return ( egrid[1] - egrid[0] ) / kt, np.abs( ratio - 1.0 ).max()

def test_g1_resolution():
    #For a pure exp(-E/kT) fall-off and binwidths of 0.25kT, deviations are
    #<0.8%, but slope discontinuities of the piecewise linear VDOS (e.g. at
    #its lowest energy) add up to ~2% more. Before G1 binwidths were limited,
    #the deviations at 5K and 14K were ~800% and ~60%.
    info = NC.createInfo('stdlib::AcrylicGlass_C5O2H8.ncmat')
    di = next( d for d in info.dyninfos
               if d.atomData.displayLabel() == 'H' )
    vdos = ( di.vdos_egrid, di.vdos_density )
    mass = di.atomData.averageMassAMU()
    for temp in ( 0.1, 1.0, 5.0, 14.0, 50.0, 293.15 ):
        egrid, g1 = nc_vdos.extractGn( vdos, n = 1, mass_amu = mass,
                                       temperature = temp )
        bw, maxdev = gn_dbcheck( egrid, g1, temp )
        print(f'T={temp:g}K: G1 binwidth/kT <= 0.25: {bw <= 0.25 + 1e-9}'
              f'  |DB deviation| < 2.5%: {maxdev < 0.025}')
        assert bw <= 0.25 + 1e-9
        assert maxdev < 0.025

    #Temperatures below 0.1K are not supported:
    with ensure_error(NC.NCBadInput,'VDOS expansion not supported for'
                      ' temperatures below 0.1K (requested T=0.09K)'):
        nc_vdos.extractGn( vdos, n = 1, mass_amu = mass, temperature = 0.09 )

    #VDOS extending to 10eV at 0.1K would need ~4.6e6 G1 bins:
    wide_vdos = ( ( 0.01, 10.0 ), [ 1.0 ] * 20 )
    with ensure_error(NC.NCCalcError):
        nc_vdos.extractGn( wide_vdos, n = 1, mass_amu = mass,
                           temperature = 0.1 )

def test_thinning():
    #Test both detailed balance of Gn functions, and that cross sections
    #agree with reference values obtained with all Gn functions kept at a
    #binwidth of kT/8.
    cfgstr = 'plain_osc.ncmat;temp=14K;comp=inelas;vdoslux=2005'
    xs_ref = [ ( 1e-05, 77.934 ), ( 3.798e-05, 44.807 ),
               ( 0.0001442, 34.099 ), ( 0.0005479, 41.112 ),
               ( 0.002081, 40.013 ), ( 0.007903, 31.025 ),
               ( 0.03001, 28.109 ), ( 0.114, 26.906 ), ( 0.433, 24.41 ),
               ( 1.644, 21.742 ) ]
    info = NC.createInfo( cfgstr )
    di = info.dyninfos[0]
    for n in ( 5, 20, 50, 100 ):
        egrid, gn = di.extract_Gn( n, without_xsect = True )
        _, maxdev = gn_dbcheck( egrid, gn, di.temperature, relfloor = 1e-6 )
        print(f'G{n}: |detailed balance deviation| < 12%: {maxdev < 0.12}')
        assert maxdev < 0.12
    sc = NC.createScatter( cfgstr )
    for e, xsref in xs_ref:
        xs = sc.xsect( e )
        ok = abs( xs / xsref - 1.0 ) < 0.015
        print(f'Ekin={e:g}eV: xs within 1.5% of {xsref:g}b: {ok}')
        assert ok

def test_resource_limits():
    #A VDOS extending to 8eV needs more than 1e6 points in G1 at a
    #temperature of 0.4K. The expansion is thus stopped before the second
    #order, leaving it unable to cover neutron energies up to the required
    #0.1eV. The same VDOS works at a higher temperature.
    c = NC.NCMATComposer()
    c.set_dyninfo_vdos( 'H', vdos_egrid = ( 0.01, 8.0 ), vdos = [ 1.0 ] * 20 )
    c.set_density( 1.0 )
    c.set_state_of_matter( 'solid' )
    NC.registerInMemoryFileData( 'widevdos.ncmat', c.create_ncmat() )
    for temp, lux, expect_error in [ ( 0.4, 2000, True ),
                                     ( 5.0, 2000, False ) ]:
        cfgstr = f'widevdos.ncmat;temp={temp}K;comp=inelas;vdoslux={lux}'
        try:
            NC.createScatter( cfgstr )
            print(f'{cfgstr}: OK')
            assert not expect_error
        except NC.NCCalcError as e:
            msg = str(e)
            print(f'{cfgstr}: CalcError: {msg.split("(")[0].strip()}')
            assert expect_error
            assert 'excessive resources' in msg


if __name__ == '__main__':
    import sys
    main( do_plot = '--plot' in sys.argv[1:] )
