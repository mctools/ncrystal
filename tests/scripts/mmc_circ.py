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

def main():
    from NCTestUtils.minimc_ref import main_minimc_unittest_stdsphere as m
    m( cfgstr='LiquidWaterH2O_T293.6K.ncmat',
       neutron_energy='1.0Aa',
       key='<auto>_h20',
       illuminate_uniformly=True )
    m( cfgstr='LiquidHeavyWaterD2O_T293.6K.ncmat',
       neutron_energy='3.0Aa',
       key='<auto>_d20',
       illuminate_uniformly=True )

if __name__ == '__main__':
    main()
