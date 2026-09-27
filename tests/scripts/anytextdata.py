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

# Tests AnyTextData initialisation, in particular error handling for invalid
# data and non-existing files.

import pathlib

import NCrystalDev as NC
import NCTestUtils.enable_fpe  # noqa: F401
from NCrystalDev.misc import AnyTextData
from NCTestUtils.common import work_in_tmpdir


def expect_badinput( fct, descr ):
    try:
        fct()
    except NC.NCBadInput as e:
        print(f'{descr}: NCBadInput: {e}')
    else:
        raise RuntimeError(f'{descr}: expected NCBadInput')

al = NC.createTextData( 'stdlib::Al_sg225.ncmat' ).rawData
with work_in_tmpdir():
    pathlib.Path( 'al.ncmat' ).write_text( al )
    for data, descr in ( ( 'al.ncmat', 'path as str' ),
                         ( pathlib.Path( 'al.ncmat' ), 'pathlib.Path' ),
                         ( al, 'in-memory str' ),
                         ( al.encode(), 'in-memory bytes' ) ):
        t = AnyTextData( data )
        assert t.content == al
        print(f'{descr}: OK (name={t.name!r})')
    expect_badinput( lambda : AnyTextData( 'no_such_file.ncmat' ),
                     'non-existing path' )
    expect_badinput( lambda : AnyTextData( 12345 ), 'int' )
    #A cfg-string is not text data (it is interpreted as a file path):
    expect_badinput( lambda : NC.directLoad( 'stdlib::Al_sg225.ncmat;'
                                             'dcutoff=1' ),
                     'directLoad with cfg-string' )
