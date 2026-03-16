#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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
    roulopts = [
        ( None,'def' ),
        ( '0.5,1.0,0',#Psurv=0.5 after all steps
          'roul50p' ),
        ( '0.001,0.001,0',#~always kill low weight
          'roulw' ),
        ( '1e-12,1.0,5',#kill ~all after 5 steps
          'roul5' ),
    ]

    for tally in ('theta','nscat','nscat_uw','w'):
        for roulopt, key in roulopts:
            print(f"Testing tally={tally} roulette={roulopt or 'default'}")
            m( extra_enginecfg=( '' if roulopt is None
                                 else f'roulette={roulopt}' ),
               key=f'<auto>_{tally}_{key}',
               cfgstr='Al_sg225.ncmat',
               neutron_energy='2.0Aa',
               tally=tally )

if __name__ == '__main__':
    main()
