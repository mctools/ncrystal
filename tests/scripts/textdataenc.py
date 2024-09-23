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

import NCrystal as NC

# A very simple test, based on something that gave problems on Windows (due to
# encodings or line endings):

def query_print( key ):
    lines=[]
    print()
    print( "===============================================================" )
    print( "===============================================================" )
    print( "===============================================================" )
    print( f"=== NC.createTextData( {repr(key)} ===" )
    td = NC.createTextData(key)
    print( "=== Printing line by line ===" )
    for line in td:
        lines.append(line)
        print(line)
    print( "=============================" )
    assert not '\r' in td.rawData
    print( "=== Printing .rawData ===" )
    print( td.rawData, end='' )
    print( "=============================" )
    return lines,td.rawData

query_print('stdlib::Ti_sg194.ncmat')

dummy_data = [
    b'NCMAT v7',
    b'@DENSITY',
    b'  1 g_per_cm3',
    b'@DYNINFO',
    b'  element  C',
    b'  fraction 1',
    b'  type     vdosdebye',
    b'  debye_temp 300',
]

dummy_data = b'<ENDOFLINE>'.join( dummy_data+[b''] )

dummy_data_decoded = dummy_data.decode('ascii').replace('<ENDOFLINE>','\n')
NC.registerInMemoryFileData( 'unix.ncmat',
                             dummy_data.replace(b'<ENDOFLINE>', b'\x0a' ) )
NC.registerInMemoryFileData( 'dos.ncmat',
                             dummy_data.replace(b'<ENDOFLINE>', b'\x0d\x0a' ) )
NC.registerInMemoryFileData( 'mixed.ncmat',
                             dummy_data
                             .replace( b'@DENSITY<ENDOFLINE>',
                                       b'@DENSITY\x0d\x0a')
                             .replace( b'<ENDOFLINE>',
                                       b'\x0a' )
                            )
NC.registerInMemoryFileData( 'weird.ncmat',
                             dummy_data
                             .replace( b'@DENSITY<ENDOFLINE>',
                                       b'@DENSITY\x0d')
                             .replace( b'<ENDOFLINE>',
                                       b'\x0a' )
                            )

res_unix =  query_print('virtual::unix.ncmat')
res_dos = query_print('virtual::dos.ncmat')
res_mixed = query_print('virtual::mixed.ncmat')
res_weird = query_print('virtual::weird.ncmat')

for enc in ('unix','dos','mixed','weird'):
    lines,rawdata = query_print(f'virtual::{enc}.ncmat')
    assert lines == dummy_data_decoded.splitlines(), f'Unexpected {enc}.ncmat lines'
    assert rawdata == dummy_data_decoded, f'Unexpected {enc}.ncmat .rawData'
