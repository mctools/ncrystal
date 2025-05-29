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
from NCrystalDev.ncmat2endf import EndfMetaData
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
ncmat2endf_impl.unit_test_not_write_version[0] = True
# We need precision to compute the XS, but it would produce
# differences in the least significant digits
ncmat2endf_impl.unit_test_dump[0] = False

d = {'matnum':{"Si":37},
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

test_cfg('Si_sg227.ncmat;vdoslux=1;incoh_elas=false', elastic_mode='scaled',
          check_edge_positions=True, endf_metadata=d,
          ref_parsed={'tsl_Si.endf':'0 0 1 451 7 2 7 4 7 451'},
          check_teff=True,
          compare_xsec=True, dump_file=False, include_gif=True)

ncmat2endf_impl.unit_test_chop_svals[0] = True
ncmat2endf_impl.is_unit_test[0] = True
ncmat2endf_impl.unit_test_dump[0] = True

test_cfg('SiO2-alpha_sg154_AlphaQuartz.ncmat', elastic_mode='scaled',
          check_edge_positions=False,
          compare_xsec=False, dump_file=False, include_gif=True, verbosity=3)
test_cfg('SiO2-alpha_sg154_AlphaQuartz.ncmat', elastic_mode='greater',
          check_edge_positions=False, verbosity=3,
          compare_xsec=False, dump_file=False, include_gif=True)
test_cfg('SiO2-alpha_sg154_AlphaQuartz.ncmat', elastic_mode='greater',
          check_edge_positions=False, verbosity=3, othertemps=400,
          compare_xsec=False, dump_file=False, include_gif=True, lasym=3)

