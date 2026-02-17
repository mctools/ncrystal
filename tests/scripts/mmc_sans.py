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

def main():
    import NCrystalDev as NC
    from NCTestUtils.minimc_ref import main_minimc_unittest_stdsphere as m

    #First a simple test that our q-values are in the right position. Bragg
    #peaks should be at q=2pi/dsp, i.e. at 2.687... for the longest dsp peak in
    #aluminium (neutron wavelength should not matter as long as it is less that
    #2dsp):
    m( cfgstr='stdlib::Al_sg225.ncmat;dcutoff=2.1;comp=bragg',
       neutron_energy='3Aa',
       sphere_diam_meter=0.001,
       tally='q',
       tallybins='q:100:0:5',
       key='<auto>_braggq')

    #Now for the SANS sphere model. First prepare Aluminium with spherical
    #voids:
    c = NC.NCMATComposer('stdlib::Al_sg225.ncmat')
    c.add_secondary_phase(0.1,'void.ncmat')
    c.add_hard_sphere_sans_model( 10 )
    c.register_as('Al_sans.ncmat')

    #Launch the test:
    m( cfgstr='virtual::Al_sans.ncmat',
       neutron_energy='5Aa',
       sphere_diam_meter=0.002,
       tally='q',
       tallybins='q:200:0.0:3.5',
       key='<auto>_sans',
       n=1e6)

if __name__ == '__main__':
    main()
