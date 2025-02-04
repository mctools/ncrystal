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

import NCrystalDev as NC
import NCrystalDev.mcstasutils as ncm

import NCTestUtils.enable_testdatapath # noqa F401

ncversion = NC.__version__
for cfg in ( ('phases<0.9*stdlib::Al_sg225.ncma'
              't&0.1*stdlib::void.ncmat;temp=50K>;vdoslux=5;temp = 200K'),
             'stdlib::Al_sg225.ncmat',
             'solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11',
             'stdlib::Polyethylene_CH2.ncmat',
             'customdirs::Al_sg225_with_voids.ncmat' ):
    for split in (True,False):
        args=dict( cfgstr = cfg,
                   name = 'mymaterial',
                   split_by_physics = split )
        hr = '='*80
        print()
        print(hr)
        print('===>',cfg)
        print('===> split =',('yes' if split else 'no'))
        print(hr)
        print()
        a,ll=ncm.cfgstr_2_unioncfg(**dict((k,v) for k,v in args.items() if k!='name'))
        print( '%.13g'%a,ll )
        print()
        print(hr)
        print()
        s = ncm.cfgstr_2_union_instrument_code( **args )
        print( s.replace(f'NCrystal v{ncversion} via',
                         'NCrystal v<current> via',) )

