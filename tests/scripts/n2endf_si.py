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
          compare_xsec=True, dump_file=False, include_gif=True)

