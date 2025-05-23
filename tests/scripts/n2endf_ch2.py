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

# fixme: add tests that dump the data before calling endf-parserpy
# NEEDS: numpy endf-parserpy

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
from NCTestUtils.ncmat2endf_utils import test_cfg
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
from NCrystalDev.ncmat2endf import EndfMetaData
ncmat2endf_impl.unit_test_chop_svals[0] = True

d = {'mat_numbers':{"C":37, "H": 38},
     'edate':'JUL01',
     'ddate':'JUL01',
     'rdate':'JUL01',
     'alab':'TestLab',
     'libname':'TestLib',
     'auth': 'Jane Doe',
     'reference': 'Aaaaaa, et al.',
     'nlib': 0,
     'nver': 4
}
m = EndfMetaData()
m.update_from_dict(d)

metadata = EndfMetaData()
metadata.set_mat_numbers( {"C":37, "H": 38} )

test_cfg('Polyethylene_CH2.ncmat;vdoslux=1', material_name='CH2',
         check_teff=True,
         ref_parsed={'tsl_H_in_CH2.endf':'0 0 1 451 7 2 7 4',
                     'tsl_C_in_CH2.endf':'0 0 1 451 7 2 7 4'},
         endf_metadata=metadata)

