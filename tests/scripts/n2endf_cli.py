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

# NEEDS: numpy endf-parserpy

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401

from NCTestUtils.ncmat2endf_utils import test_cli
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
from pathlib import Path


#
# CLI tests
#

ncmat2endf_impl.unit_test_chop_svals[0] = True
if Path('tsl_Al.endf').is_file():
    Path('tsl_Al.endf').unlink()
test_cli('"stdlib::Al_sg225.ncmat;temp=350K;vdoslux=1"'
         ' -vvv -n Al -f -e greater'
        r""" --mdata='{"ALAB": "MyLab"}'""")
assert Path('tsl_Al_in_Al.endf').is_file()
Path('tsl_Al_in_Al.endf').unlink()
test_cli('-h')
test_cli('--mdata=help')
test_cli('--mdata','help')
assert not Path('tsl_Al.endf').is_file()
test_cli('"stdlib::Al_sg225.ncmat;vdoslux=1"'
         ' -vvv -f -e scaled --totsab --asymsab')
assert Path('tsl_Al.endf').is_file()
assert not Path('some/out/dir/tsl_Al.endf').is_file()
test_cli('"stdlib::Al_sg225.ncmat;vdoslux=1" -f --now --dir=some/out/dir')
assert Path('some/out/dir/tsl_Al.endf').is_file()
