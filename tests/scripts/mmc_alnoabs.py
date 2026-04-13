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

# NEEDS: numpy spglib

import NCTestUtils.enable_fpe # noqa F401

def main():
    import NCrystalDev as NC
    #Prepare Aluminium without absorption (NB: we could also use the MiniMC
    #enginecfg!)
    c = NC.NCMATComposer('stdlib::Al_sg225.ncmat')
    a = NC.atomDB('Al')
    c.update_atomdb( 'Al',
                     mass = a.averageMassAMU(),
                     coh_scat_len = a.coherentScatLenFM(),
                     incoh_xs = a.incoherentXS(),
                    abs_xs = 0.0 )
    c.register_as('Al_noabsn.ncmat')

    #Launch the test:
    from NCTestUtils.minimc_ref import main_minimc_unittest_stdsphere as m
    m( cfgstr='virtual::Al_noabsn.ncmat',
       neutron_energy='2.5Aa',
       sphere_diam_meter=0.2 )

if __name__ == '__main__':
    main()
