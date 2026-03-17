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

    #Test basic aluminium at 4Aa and 0.8Aa, as well as temperature
    #dependence. Also test the tallybreakdown parameter:

    #4Aa + 700K:
    from NCTestUtils.minimc_ref import main_minimc_unittest_stdsphere as m
    m( cfgstr='Al_sg225.ncmat;temp=700K',
       neutron_energy='4.0Aa', key='<auto>_4Aa700K' )

    #4Aa at default temp:
    res = m( cfgstr='Al_sg225.ncmat',
             neutron_energy='4.0Aa', key='<auto>_4Aa' )
    assert res.tally_names == ['theta']
    assert ( list(res.tally('theta').hist_breakdown.keys())
             == ['NOSCAT', 'SINGLESCAT_ELAS', 'SINGLESCAT_INELAS',
                 'MULTISCAT_PUREELAS', 'MULTISCAT_OTHER'] )

    #Same at 0.8Aa +tallybreakdown=0.
    res = m( cfgstr='Al_sg225.ncmat',neutron_energy='1.0Aa', key='<auto>_0.8Aa',
             extra_enginecfg = 'tallybreakdown=0'  )
    assert res.tally_names == ['theta']
    assert res.tally('theta').hist_breakdown is None

if __name__ == '__main__':
    main()
