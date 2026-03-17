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
    from NCTestUtils.minimc_ref import main_minimc_unittest as m
    def test( key, srccfg_energypart, hist_elow, hist_emax ):
        kwargs = dict( cfgstr = 'void.ncmat',
                       srccfg = f'constant;z=-0.001;{srccfg_energypart};n=1e6',
                       geomcfg = 'slab;dz=0.001',
                       key=f'<auto>_{key}',
                       tally='e',
                       tallybins=f'e:100:{hist_elow}:{hist_emax}',
                       extra_enginecfg='nscatlimit=0;absorption=0' )
        m(**kwargs)
        #Again, with auto-binning
        kwargs['srccfg'] += ';n=1e4'
        kwargs['tally'] = 'q,e,de,l'
        kwargs['tallybins'] = None
        kwargs['key'] += '_autobin'
        m(**kwargs)
    test('fixe','ekin=0.025', 0.0, 0.04 )
    test('fixwl','wl=1.8', 0.0, 0.04 )
    test('rangee','ekin=0.02-0.03', 0.0, 0.04 )
    test('rangewl','wl=1.0-5e0', 0.0, 0.1 )
    test('lognormale','ekin=0.025+-0.0025', 0.0, 0.06 )
    test('lognormalwl','wl=1.8+-0.18', 0.0, 0.06 )


if __name__ == '__main__':
    main()
