#!/usr/bin/env python3

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

# NEEDS: numpy

import NCrystalDev._mmc as ncmmc #fix import
from NCrystalDev.core import NCBadInput
from NCTestUtils.common import ensure_error

def testbad(cfgstr,scenario, expecterr):
    test(cfgstr,scenario, expecterr)

def test(cfgstr,scenario, expecterr = None):
    print()
    print(f'Testing "{scenario}"')
    print(f'  with matcfg "{cfgstr}"')
    def f():
        return ncmmc.minimc_decode_scenario( cfgstr, scenario )

    if expecterr is not None:
        with ensure_error(NCBadInput,expecterr):
            s = f()
        return
    s = f()
    assert isinstance(s,dict)
    assert all(isinstance(k,str) for k,v in s.items())
    assert all(isinstance(v,str) for k,v in s.items())
    assert set(s.keys())==set(['cfgstr','geomcfg','srccfg','enginecfg'])
    assert s['cfgstr'] == cfgstr
    print('  -> geomcfg   = "%s"'%s['geomcfg'])
    print('  -> srccfg    = "%s"'%s['srccfg'])
    print('  -> enginecfg = "%s"'%s['enginecfg'])



def main():
    c = 'stdlib::Al_sg225.ncmat'
    c_mp = ( 'phases<0.001*stdlib::C_sg194_pyrolytic_graphite.ncmat'
             '&0.999*solid::Al/1gcm3>' )
    test(c,'1.0Aa on 2.0mm')
    test(c,'1e-3Aa on 1e2mm')
    testbad(c,'1e-3 Aa on 1e2mm',
            'Invalid MiniMC scenario string: "1e-3 Aa on 1e2mm".'
            ' Expected keyword "on" but got "Aa".')

    testbad(c,'1e-3Aa on 1e2 mm',
            'Invalid MiniMC scenario string: "1e-3Aa on 1e2 mm". '
            'Unexpected parameter: "mm"')

    testbad(c,'1.0 on 2.0mm',
            'Missing energy unit on "1.0". Possible units include Aa, meV, '
            'eV, or BT (=Bragg threshold in Aa, or 4Aa if not available).'
            )

    testbad(c,'1.0ergs on 2.0mm',
            'Invalid energy unit on "1.0ergs". Possible units include Aa, meV, '
            'eV, or BT (=Bragg threshold in Aa, or 4Aa if not available).'
            )
    test(c,'1.1BT on 2.0mm')
    test(c,'1.8Aa on 0.5mfp')
    test(c,'')
    test(c,'1.8meV')
    testbad(c,'[',
            'Forbidden character "[" in MiniMC scenario string: "["')

    weird_char = '\u2300'
    testbad(c,f'4.0{weird_char} on 1mm',
            'Forbidden character "\\x0-30" in MiniMC scenario'
            ' string: "4.0⌀ on 1mm"')

    test(c,'1.0Aa on 2.0mm')
    test(c,'1.0Aa pencil on 2.0mm')
    test(c,'1Aa on slab')
    test(c,'25meV on 2.0mm slab')
    test(c,'25meV pencil on 2.0mm slab')

    test(c,'  25meV_pencil_on:2.0mm:::slab:::_')

    testbad(c,'1.8Aa on 0.5miles',
            'invalid length: "0.5miles". Must be a value followed by'
            ' either the special unit "mfp" (mean-free-path between'
            ' scatterings) or one of the standard length units: Aa, '
            'nm, mu, mm, cm, m');

    testbad(c,'1.8Aa on 0.5e-3.2cm',
            'Invalid thickness specification in "0.5e-3.2cm".')

    testbad(c,'1eV on 1e-250Aa',
            'Could not determine suitable thickness from "1e-250Aa"'
            ' (thickness is invalid or out of range)')

    testbad(c,'1eV2 on 1mm','Invalid energy specification in "1eV2".')

    testbad(c,'1eV on 1mm2','Invalid thickness specification in "1mm2".')

    testbad(c,'1eV on 17',
            'Missing length unit on: "17". Accepts either the special unit'
            ' "mfp" (mean-free-path between scatterings) or one of the '
            'standard length units: Aa, nm, mu, mm, cm, m')

    testbad(c,'1Aa on','Invalid MiniMC scenario string: "1Aa on". '
            'Missing parameters after keyword "on".')

    testbad(c,'1Aa on slab 1mm: ','Invalid MiniMC scenario string:'
            ' "1Aa on slab 1mm: ". Unexpected parameter: "1mm"')

    test(c,'1BT on 1mm')
    test(c_mp,'1BT on 1mm')


if __name__ == '__main__':
    main()
