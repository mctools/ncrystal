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

from NCrystalDev.exceptions import NCBadInput
from NCrystalDev.ncmat2endf import EndfMetaData
from NCTestUtils.ncmat2endf_utils import ( test_cfg_fail,
                                           test_cli )

import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
ncmat2endf_impl.unit_test_chop_svals[0] = True


#
# Error handling tests
#

# Oriented materials not supported
test_cfg_fail( NCBadInput, 'Ge_sg227.ncmat;dcutoff=0.5;mos=40arcsec;'
               'dir1=@crys_hkl:5,1,1@lab:0,0,1;'
               'dir2=@crys_hkl:0,-1,1@lab:0,1,0')
# Wrong material number assignment
metadata = EndfMetaData()
metadata.set_mat_numbers( {"Ge":99} )
test_cfg_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1', endf_metadata=metadata)
# Negative temperatures
test_cfg_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1', temperatures=[-100])
# Repeated temperatures
test_cfg_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1', temperatures=[293.15])
# No inelastic
test_cfg_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1;comp=coh_elas')
# Wrong elastic mode
test_cfg_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1',
               elastic_mode='something wrong')

#
# CLI tests
#

test_cli('"stdlib::Al_sg225.ncmat;temp=350K;vdoslux=1"'
         ' -vvv -m Al -f -e greater')

#fixme: more CLI tests
#fixme: add tests that dump the data before calling endf-parserpy
