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
from NCTestUtils.ncmat2endf_utils import test_cfg
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl

ncmat2endf_impl.unit_test_chop_vals[0] = True

test_cfg('stdlib::Al_sg225.ncmat;vdoslux=1', material_name='Al',
         check_teff=True,
         ref_parsed={'tsl_Al_in_Al.endf':'0 0 1 451 7 2 7 4'},
         check_edge_positions=True,
         othertemps=350, elastic_mode='scaled', compare_xsec=False,
         dump_file=True)
