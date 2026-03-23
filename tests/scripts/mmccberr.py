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

from NCrystalDev.minimc import run as mmcrun
from NCTestUtils.common import ensure_error

def main():

    #Test errors during simulation, in particular in the callback function
    #itself.

    n_tot = [0]
    n_cb = [0]
    cb_opt = ['std']
    def cb( data ):
        n_cb[0] += 1
        w = data['w']
        n_tot[0] += len(w)
        if n_tot[0] < 5e4 or cb_opt[0] == 'std':
            print("UserCB: Emitting STD after %i neutrons"%n_tot[0])
            return False if n_tot[0]%2 else None
        elif cb_opt[0] == 'error':
            print("UserCB: Emitting HALTERR after %i neutrons"%n_tot[0])
            ekin0 = data['ekin0']#should not be there with basket_type=basic and
                                 #tally=e, so triggers an exception
            n_tot[0] += len(ekin0)#this line should never happen
        else:
            assert cb_opt[0] == 'haltsrc'
            if n_tot[0]%2:
                #It should not make a difference if we continue to emit haltsrc
                #or not:
                cb_opt[0]='std'
            print("UserCB: Emitting HALTSRC after %i neutrons"%n_tot[0])
            return True

    cb_opt = ['error']
    n_tot[0] = 0
    with ensure_error(KeyError,'ekin0'):
        res = mmcrun( "Al_sg225.ncmat", scenario='2Aa on 2cm 2e5 times',
                      callback = cb,
                      enginecfg = 'nthreads=2;tally=e',
                      callback_options='basket_type=basic;cachelen=1e4' )
    print("ran callback %i times (%i neutrons)"%(n_cb[0],n_tot[0]))
    assert 5 <= n_cb[0] <= 10 # very rough, saw 7 in tests
    assert 50e3 <= n_tot[0] <= 65e3 # very rough, saw 57k in tests

    #Try again, this time halting source:
    cb_opt = ['haltsrc']
    n_tot[0] = 0
    res = mmcrun( "Al_sg225.ncmat", scenario='2Aa on 2cm 2e5 times',
                  callback = cb,
                  enginecfg = 'nthreads=2;tally=e',
                  callback_options='basket_type=basic;cachelen=1e4' )
    print("ran callback %i times (%i neutrons)"%(n_cb[0],n_tot[0]))
    assert 14 <= n_cb[0] <= 22 # very rough, saw 17-18 in tests
    assert 58e3 <= n_tot[0] <= 92e3 # very rough, saw 73k in tests
    print(res.output_metadata['tallied']['count'],n_tot[0])
    assert res.output_metadata['tallied']['count']==n_tot[0]

    #Try again, this time with no error:
    cb_opt = ['std']
    n_tot[0] = 0
    res = mmcrun( "Al_sg225.ncmat", scenario='2Aa on 2cm 2e5 times',
                  callback = cb,
                  enginecfg = 'nthreads=2;tally=e',
                  callback_options='basket_type=basic;cachelen=1e4' )
    print("ran callback %i times (%i neutrons)"%(n_cb[0],n_tot[0]))
    assert 80 <= n_cb[0] <= 120 # very rough, saw 100 in tests
    assert 550e3 <= n_tot[0] <= 900e3 # very rough, saw 656k in tests
    print(res.output_metadata['tallied']['count'],n_tot[0])
    assert res.output_metadata['tallied']['count']==n_tot[0]

    #Test a few other things:
    cb_opt = ['std']
    n_tot[0] = 0
    from NCrystalDev.exceptions import NCBadInput
    with ensure_error(NCBadInput,
                      'Invalid cachelen "1000000001"'
                      ' (must be less than 1000000000)'):
        res = mmcrun( "Al_sg225.ncmat", scenario='2Aa on 2cm 2e5 times',
                      callback = cb,
                      enginecfg = 'nthreads=2;tally=e',
                      callback_options='cachelen=1000000001' )


if __name__ == '__main__':
    main()
