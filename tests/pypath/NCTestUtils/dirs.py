
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

   Module providing directories needed during testing.

"""

from pathlib import Path as _Path

_pymoddir = _Path(__file__).resolve().absolute().parent

def _find_data_dir():
    from .modeinfo import is_simplebuild_mode, is_ncrystalverify_mode
    if is_simplebuild_mode():
        import os
        ddir = _Path(os.environ['SBLD_DATA_DIR'])/'NCTestUtils'
    elif is_ncrystalverify_mode():
        ddir = _pymoddir.parent.parent.joinpath('data')
    else:
        ddir = _pymoddir.parent.parent.joinpath('data')
    assert ddir.is_dir()
    return ddir

test_data_dir = _find_data_dir()

def get_named_test_data_dir(name):
    #name of subdir of /tests/data
    from .modeinfo import is_simplebuild_mode, is_ncrystalverify_mode
    if is_simplebuild_mode():
        import os
        ddir = _Path(os.environ['SBLD_DATA_DIR'])/f'NCTestData_{name}'
    elif is_ncrystalverify_mode():
        ddir = _pymoddir.parent.parent.joinpath('data',name)
    else:
        ddir = _pymoddir.parent.parent.joinpath('data',name)
    assert ddir.is_dir()
    return ddir
