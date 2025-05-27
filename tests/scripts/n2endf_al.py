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
import NCrystalDev.ncmat as nc_ncmat

ncmat2endf_impl.unit_test_chop_svals[0] = True

test_cfg('Al_sg225.ncmat;vdoslux=1', material_name='Al',
         check_teff=True,
         ref_parsed={'tsl_Al.endf':'0 0 1 451 7 2 7 4'},
         check_edge_positions=True,
         temperatures=350, elastic_mode='scaled', compare_xsec=False,
         dump_file=True)

ncmat2endf_impl.unit_test_chop_svals[0] = True
ncmat2endf_impl.is_unit_test[0] = True
ncmat2endf_impl.unit_test_dump[0] = True

test_cfg('AlN_sg186_AluminumNitride.ncmat;vdoslux=1', material_name='AlN',
         temperatures=350, elastic_mode='scaled', compare_xsec=False,
         dump_file=False, verbosity=3)
test_cfg('Al2O3_sg167_Corundum.ncmat;vdoslux=1', material_name='Al2O3',
         elastic_mode='mixed', compare_xsec=False,
         dump_file=False, verbosity=3)
test_cfg('Al4C3_sg166_AluminiumCarbide.ncmat;vdoslux=1', material_name='Al4C3',
         temperatures=350, elastic_mode='greater', compare_xsec=False,
         dump_file=False, verbosity=3)
# Long config strings
test_cfg('Al_sg225.ncmat;temp=300K;temp=300K;temp=300K;'
         'temp=300K;temp=300K;temp=300K;temp=300K;temp=300K;'
         'temp=300K;temp=300K;vdoslux=1', material_name='Al',
         elastic_mode='scaled', compare_xsec=False,
         dump_file=False)
c=nc_ncmat.NCMATComposer('Al_sg225.ncmat')
longname = ('Al_with_a_Very_Very_Ver_Very_Very_'
            'Very_Very_Ver_Very_Very_Long.ncmat')
verylongname = ('Al_with_a_Very_Very_Ver_Very_Very_Very_Very_'
                'Ver_Very_Very_Very_Very_Ver_Very_Very_Long_Name.ncmat')
c.write(longname)
c.write(verylongname)
test_cfg(longname, material_name='Al',
         elastic_mode='scaled', compare_xsec=False,
         dump_file=False)
test_cfg(verylongname+';temp=300K;vdoslux=1', material_name='Al',
         elastic_mode='scaled', compare_xsec=False,
         dump_file=False)

