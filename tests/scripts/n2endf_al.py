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
ncmat2endf_impl.unit_test_chop_svals[0] = True

# fixme: Replace with NCrystal values
ref_Eint = (0.003741252, 0.004988336, 0.009976672, 0.01371792, 0.01496501,
            0.01995334, 0.0236946, 0.02494168, 0.02993002, 0.03367127,
            0.03990669, 0.04364794, 0.04489502, 0.04988336, 0.05362461,
            0.05487169, 0.05986003, 0.06360128, 0.06484837, 0.0698367)
ref_S0 = ((0.005134741, 0.00839204, 0.01258346, 0.01924357, 0.02131949,
           0.02254634, 0.02674178, 0.03073558, 0.03405288, 0.03793785,
           0.03912707, 0.04336451, 0.04591494, 0.04767616, 0.04925869,
           0.05078662, 0.05123031, 0.05363637, 0.05479998, 0.05684046),
           ((0.005071159,), (0.00827479,), (0.01232923,), (0.0186919,),
            (0.0206669,), (0.02181487,), (0.02569196,), (0.02936742,),
            (0.03236999,), (0.03584283,), (0.03688405,), (0.04054818,),
            (0.04274442,), (0.04423605,), (0.04555975,), (0.04683249,),
            (0.04719598,), (0.04914274,), (0.05008033,), (0.05169737,)))

test_cfg('Al_sg225.ncmat;vdoslux=1', material_name='Al',
         check_teff=True,
         ref_parsed={'tsl_Al.endf':'0 0 1 451 7 2 7 4'},
         ref_bragg_edges={'tsl_Al.endf':(ref_Eint, ref_S0)},
         temperatures=[350], elastic_mode='scaled', compare_xsec=False)

