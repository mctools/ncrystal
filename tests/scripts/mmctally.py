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
    #Try to test all tallies:
    from NCrystalDev.core import enableFactoryThreads
    from NCrystalDev.minimc import tally_info
    from NCTestUtils.minimc_ref import main_minimc_unittest_stdsphere as m
    from NCrystalDev.constants import wl2ekin

    enableFactoryThreads(3)
    tallies = ('de', 'e', 'l', 'mu', 'nscat', 'nscat_uw', 'q', 'theta', 'w')
    import pprint
    pprint.pp(tally_info())
    assert tallies == tally_info()['flags']['ALLHISTS']

    #Test all tallies:
    for tally in tallies:
        print(f"Testing tally: {tally}")
        m( cfgstr='Al_sg225.ncmat;temp=200K',
           neutron_energy='2.5Aa',
           key=f'<auto>_{tally}',
           tally=tally )

    #Test tallyref on the tallies possibly affected:

    #NB: We have no 'divergence' in src direction currently, apart from
    #    the isotropic source, so just testing effect of energy spread.

    for tally in ('de','mu','theta','q'):
        kw = dict(cfgstr='Al_sg225.ncmat;temp=200K',
                  neutron_energy=(wl2ekin(2.5),'wl=2.5+-0.1'),
                  tally=tally)
        for tr in ('src','truth'):
            print(f"Testing tally: {tally} [tallyref={tr}]")
            m( key=f'<auto>_{tally}_tallyref{tr}',
               extra_enginecfg=f'tallyref={tr}',
               **kw)

if __name__ == '__main__':
    main()
