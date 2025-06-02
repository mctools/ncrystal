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

from argparse import ArgumentError
from NCrystalDev.ncmat2endf import EndfMetaData
from NCTestUtils.ncmat2endf_utils import ( test_cfg_fail,
                                           test_cli_fail )
from NCTestUtils.dirs import test_data_dir

from NCrystalDev.exceptions import NCBadInput
import NCrystalDev.ncmat as nc_ncmat
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
from pathlib import Path
from NCrystalDev._numpy import _np

# Check setting the date to now without logging the current date
metadata = EndfMetaData()
metadata.set_all_dates_as_now()
MONTHS = ('jan', 'feb', 'mar', 'apr', 'may', 'jun',
          'jul', 'ago', 'sep', 'oct', 'nov', 'dec')
assert metadata.edate[:3].lower() in MONTHS
assert int(metadata.edate[3:].lower()) >= 00
assert int(metadata.rdate[3:].lower()) <= 99
assert metadata.edate.lower() == metadata.ddate.lower()
assert metadata.edate.lower() == metadata.rdate.lower()

#
# Error handling tests
#

# Oriented materials not supported
test_cfg_fail( 'Ge_sg227.ncmat;dcutoff=0.5;mos=40arcsec;'
               'dir1=@crys_hkl:5,1,1@lab:0,0,1;'
               'dir2=@crys_hkl:0,-1,1@lab:0,1,0')
# Wrong material number assignment
metadata = EndfMetaData()
metadata.set_value('MATNUM', {"Ge":99} )
metadata.set_value('ALAB', None )

#test a few conversion functions, including that repr(metadata) can be evaluated
#as an EndfMetaData object:
assert EndfMetaData(metadata).to_json() == metadata.to_json()
assert EndfMetaData(metadata.to_dict()).to_json() == metadata.to_json()
assert eval(repr(metadata)).to_json() == metadata.to_json() # round trip
print(metadata)

test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1', endf_metadata=metadata)

# Negative temperatures
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1', othertemps=[-100])
# Repeated temperatures
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1', othertemps=[293.15])
# No inelastic
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1;comp=coh_elas')
# Wrong verbosity parameter
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1', verbosity='WRONGVALUE')
# Wrong elastic mode
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               elastic_mode='something wrong')
# Conversion implemented only for natural elements
test_cfg_fail( 'MgD2_sg136_MagnesiumDeuteride.ncmat;vdoslux=1')
# Wrong temperature parameter
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               othertemps='WRONGVALUE')
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               othertemps=['WRONGVALUE'])
# SANS not supported (currently being caught as multiphase)
c=nc_ncmat.NCMATComposer('Al_sg225.ncmat')
c.add_secondary_phase(0.01,'void.ncmat')
c.add_hard_sphere_sans_model(50)
c.write('Al_with_SANS.ncmat')
test_cli_fail( 'Al_with_SANS.ncmat')
# Conversion only supported for VDOS and VDOSDebye dyninfos
test_cfg_fail( 'freegas::Ar/2.5e-5perAa3;vdoslux=1')
# Isotopic expansion requires generalized information file
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1', include_gif=False,
               isotopic_expansion=True)
# Isotopic expansion not implemented
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1', include_gif=True,
               isotopic_expansion=True)

# Try to set wrong parameter in metadata
d = {'WRONGPARAM':0}
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               endf_metadata=d)
# Invalid symbol in metadata
d = {'ALAB':r'WRONGSYMBOL"'}
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               endf_metadata=d)
# Invalid material number assignment
d = {'MATNUM': 0}
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               endf_metadata=d)
d = {'MATNUM': {'Al':'0'}}
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               endf_metadata=d)
# Try to set wrong parameter for dates
d = {'EDATE':'WRONGDATE'}
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               endf_metadata=d)
d = {'EDATE':0}
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               endf_metadata=d)
# Invalid data type for metadata
d = {'ALAB':0}
test_cfg_fail( 'Al_sg225.ncmat;vdoslux=1',
               endf_metadata=d)
# Try to read wrong parameter from metadata
try:
    metadata.get_value('WRONGPARAM')
except NCBadInput as e:
    print("FAILED (as expected): %s"%e)
else:
    raise SystemExit('Did not fail as expected')
# Incompatible arguments in CLI
test_cli_fail( '"stdlib::Al_sg225.ncmat"', '-v', '-q',
              exception_type=ArgumentError)
# Wrong metadata in CLI
test_cli_fail( '"stdlib::Al_sg225.ncmat"',
              '--mdata', 'WRONGINPUT', exception_type=ArgumentError )
# Multiphase material
test_cli_fail( '"phases<0.1*Al_sg225.ncmat&0.9*Si_sg227.ncmat>"')
# File already exists
Path('tsl_Al_in_Bla.endf').touch()
test_cli_fail('"stdlib::Al_sg225.ncmat;vdoslux=1"'
              ' -q -n Bla -e mixed --totsab', exception_type=RuntimeError)
test_cli_fail('"stdlib::Al_sg225.ncmat;vdoslux=1"'
              ' -q -n Bla -e scaled --asymsab', exception_type=RuntimeError)

Path('Al_multirole.ncmat').write_text(
    test_data_dir.joinpath('Al_multirole.ncmat').read_text()
)

test_cli_fail('Al_multirole.ncmat')

#
# Unit tests
#
assert (ncmat2endf_impl._interp2d(0.5,1,[0,1],[0,1],
        _np.array([[0,0],[1,1]]))[0,0] == 0.5)
assert (ncmat2endf_impl._interp2d([1],[0.5],[0,1],[0,1],
        _np.array([[0,1],[0,1]]))[0,0] == 0.5)
cfg = 'ThO2_sg225_ThoriumDioxide.ncmat;vdoslux=1'
data = ncmat2endf_impl.NuclearData(ncmat_cfg=cfg,
                                   temperatures=(293.15,),
                                   elastic_mode='scaled',
                                   requested_emax=1.0,
                                   verbosity=3)
cfg = 'SiO2-alpha_sg154_AlphaQuartz.ncmat;comp=inelas'
data = ncmat2endf_impl.NuclearData(ncmat_cfg=cfg,
                                   temperatures=(293.15,),
                                   elastic_mode='scaled',
                                   requested_emax=1.0,
                                   verbosity=3)
cfg = 'Al_sg225.ncmat;coh_elas=false'
data = ncmat2endf_impl.NuclearData(ncmat_cfg=cfg,
                                   temperatures=(293.15,),
                                   elastic_mode='mixed',
                                   requested_emax=1.0,
                                   verbosity=3)
cfg = 'Al_sg225.ncmat;coh_elas=false'
data = ncmat2endf_impl.NuclearData(ncmat_cfg=cfg,
                                   temperatures=(293.15,),
                                   elastic_mode='greater',
                                   requested_emax=1.0,
                                   verbosity=3)
cfg = 'Al_sg225.ncmat;incoh_elas=false'
data = ncmat2endf_impl.NuclearData(ncmat_cfg=cfg,
                                   temperatures=(293.15, 300),
                                   elastic_mode='mixed',
                                   requested_emax=1.0,
                                   verbosity=3)
data._combine_temperatures = True
data._get_alpha_beta_grid()
cfg = 'V_sg229.ncmat;vdoslux=1'

ncmat2endf_impl.unit_test_chop_vals[0] = True

data = ncmat2endf_impl.NuclearData(ncmat_cfg=cfg,
                                   temperatures=(293.15,),
                                   elastic_mode='scaled',
                                   requested_emax=1.0,
                                   verbosity=3)

ncmat2endf_impl.unit_test_chop_vals[0] = False

cfg = 'SiO2-alpha_sg154_AlphaQuartz.ncmat;coh_elas=false'
data = ncmat2endf_impl.NuclearData(ncmat_cfg=cfg,
                                   temperatures=(293.15,),
                                   elastic_mode='scaled',
                                   requested_emax=1.0,
                                   verbosity=3)
