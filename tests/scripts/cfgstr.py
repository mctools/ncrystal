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

# NEEDS: numpy

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev.cfgstr as nc_cfgstr
from NCrystalDev.exceptions import NCBadInput, NCFileNotFound
from NCTestUtils.common import ensure_error

def main():
    import pprint
    for c in [ 'stdlib::Al_sg225.ncmat',
               'Al_sg225.ncmat;temp=293.15K',
               'Al_sg225.ncmat;temp=300K',
               'Al_sg225.ncmat;temp=20C',
               'Al_sg225.ncmat;temp=293.15K; vdoslux =3;dcutoff=   0.1 nm',
               'Al_sg225.ncmat;temp=300K',
               'C_sg194_pyrolytic_graphite.ncmat',#embedded cfg
               ('phases<0.7*C_sg227_Diamond.ncmat&0.28*solid::V/6.12gcm3'
                '&0.02*Polyethylene_CH2.ncmat;comp=inelas>;temp=380K'),
               ('phases<0.5*solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11'
                '&0.5*gasmix::BF3/2bar>'),
               'phases<1.0*C_sg227_Diamond.ncmat>;temp=380K',
               ('phases<0.4*C_sg227_Diamond.ncmat'
                '&0.6*C_sg227_Diamond.ncmat>;temp=380K'),
              ]:
        print(f"Decoding cfg {repr(c)}:")
        pprint.pprint(nc_cfgstr.decodeCfg(c),
                      indent=4,
                      sort_dicts=False)

    with ensure_error(NCFileNotFound,
                      'Could not find data: "foobarnonexistent.ncmat"'):
        nc_cfgstr.decodeCfg('foobarnonexistent.ncmat;temp=300K')

    with ensure_error(NCBadInput,
                      'Syntax error - invalid value'
                      ' "warm" provided for parameter "temp"'):
        nc_cfgstr.decodeCfg('stdlib::Al_sg225.ncmat;temp=warm')

if __name__ == '__main__':
    main()
