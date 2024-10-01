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

import NCTestUtils.enable_fpe
import NCrystal as NC

def test(cfgstr,expectBad = False):
    print("="*100)
    print(f'\nTrying to load "{cfgstr}"\n')
    try:
        td = NC.createTextData(cfgstr)
    except NC.NCBadInput as err:
        if not expectBad:
            raise
        else:
            print(f'Got expected error: >>>{err}<<<')
            return
    if expectBad:
        raise SystemExit("Did not end with exception as expected!")
    print(f"Got TextData(uid={td.uid}) with contents:")
    for l in td:
        print("---->",l)
    print()
testbad = lambda cfgstr : test(cfgstr,True)

test('solid::Pb/+2.02e-30gcm3')
test('solid::Pb/+2.02e-30gcm3')#same uid
test(' solid ::  Pb/+2.02e-30gcm3 ')#same uid (spaces here are discarded)
test('solid::Pb/ +2.02e-30gcm3')#too bad, gives different uid since keys differ in internal whitespace!
                                #This would need in-factory cache to avoid and is likely rare
testbad('solid::Pb/+2.02e-30g_per_cm3')
testbad('solid::Pb/2.02kgcm3')

test('solid::Al2O3/0.06perAa3/TDebye100K')
testbad('solid::Al2O3/0.06peraa3/TDebye100K')
test('solid::Al2O3/2.0kgm3/TDebye100K')
testbad('solid::Al2O3/2.0kgm3/TDebye100')
testbad('solid::Al2O3/2.0kgm3/TDebye')
testbad('freegas::Al2O3/0.06perAa3/TDebye100K')
test('freegas::Ar/0.06perAa3')

test('solid::B4C/2.52gcm3')
test('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05_B11')
test('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05_B11')
testbad('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05_____B11')
testbad('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05 B11')
testbad('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05:B11')
test('solid::B4C/2.383gcm3/   B_is_0.95_B10_0.05_B11  ')
test('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05_B11/Al_26.98u_0.3449fm_0.0082b_0.231b')

testbad('solid::B4C/2.383gcm3/B is 0.95 B10 0.05 B11')
testbad('solid::B4C/2.383gcm3/B:is:0.95:B10:0.05:   :B11')
testbad('solid::B4C/2.383gcm3/B is 0.95 B10 0.05 B11/Al 26.98u 0.3449fm:0.0082b:0.231b')

testbad('solid::B4C/2.52gcm3/TDebye123K_U')
testbad('solid::B4C/2.383gcm3/B:is:0.95:B10:0.05:2B11')

test('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05_B11/TDebye100K_C/TDebye123.44K_B')

testbad('solid::B4C/2.52gcm3'+(b'\xc3\xa6\xc3\xb8\xc3\xa5'.decode('utf-8')))#three danish letters
testbad('solid::B4C/2.52gcm3"')
testbad('solid::Al2O3/2.0kgm3/TDebye-100K')
testbad('solid::Al2O3/2.0kgm3/TDebye0.0K')
testbad('solid::Al2O3/2.0kgm3/TDebyeinfK')
testbad('solid::Al2O3/2.0kgm3/TDebyenanK')
testbad('solid::Al2O3/2.0kgm3/TDebye1e-500K')
test('solid::Al2O3/2.0kgm3/TDebye1e-300K')
testbad('solid::Al2O3/2.0kgm3/TDebye300K/TDebye300K')

testbad('solid::Al2O3/2.0kgm3/TDebye200K_AL')
testbad('solid::Gt2O3/2.0kgm3/TDebye200K_Gt')
testbad('solid::Al2O3/2.0kgm3/1.0kgm3')
testbad('solid::Al2O3/1.0kgm3/1.0kgm3')
testbad('solid::Al2O3/0.0kgm3')
testbad('solid::Al2O3/infkgm3')
testbad('solid::Al2O3/Al/1.0kgm3')

testbad('solid::B4C/2.383gcm3/B_is_0.95_B10_0.05_B11_0.1')
testbad('solid::B4C/2.383gcm3/B_is_0.95_B10_1.05_B11')
testbad('solid::2.383gcm3')
testbad('solid::B4C')

test('solid::Al2O3/2.0kgm3/TDebye200K_Al')
test('solid::Al2O3/2.0kgm3/TDebye200K_Al/TDebye1000K_O')
testbad('solid::Al2O3/2.0kgm3/TDebye200K_Al/TDebye1000K_B')
testbad('solid::Al2O3/2.0kgm3/TDebye200K_Al/TDebye1000K')

for fn in ('solid','freegas','gasmix'):
    for f in NC.browseFiles(factory=fn):
        test(f.fullKey)
