#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

import NCTestUtils.enable_fpe
import NCrystalDev as NC

def test_dis(info):
    for i,di in enumerate(info.dyninfos):
        if i%2:
            di.atominfo#pre-access, to check caching
        print(di.atomIndex,di,di.atominfo)

def test_ais(info):
    for i,ai in enumerate(info.atominfos):
        if i%2:
            ai.dyninfo#pre-access, to check caching
        print (ai.atomIndex, ai, ai.dyninfo )

for d in ( 'Na4Si3Al3O12Cl_sg218_Sodalite.ncmat',
           'LiquidHeavyWaterD2O_T293.6K.ncmat',
           'He_Gas_STP.ncmat',
           'Ti_sg194.ncmat' ):
    print('-------------',d)
    for di_first in (True,False):
        if not di_first:
            print('----------')
        NC.clearCaches()
        info = NC.createInfo(f'{d};dcutoff=0.8')
        if di_first:
            test_dis(info)
            test_ais(info)
        else:
            test_ais(info)
            test_dis(info)
