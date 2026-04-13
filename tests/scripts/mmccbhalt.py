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

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.minimc import run as mmcrun
from NCrystalDev.hist import HistFiller1D
from NCrystalDev.constants import ekin2wl

def main(do_plot):

    #Test our ability to halt the simulation once we feel we have collected
    #enough statistics. We will shoot a thermal source of neutrons on 5mm cold
    #PE, and stop only once we have collected a sufficient amount of
    #long-wavelength (>15Aa) neutrons headed roughly forward (~20deg cone).

    #hist + running counter of how many "useful" neutrons we collected:
    hist_lambda = HistFiller1D(100,0.0,50.0,'lambda (fwd)')
    n_useful_tot = [0]

    #Define data processing callback function which fills the histograms:
    def cb( data ):
        w, wl = data['w'], ekin2wl(data['ekin'])
        mask_fwd = data['uz'] > 0.94 # ~20 degree cone
        n_useful = ( mask_fwd & (wl>15.0) ).sum()
        print(f'Callback processing {len(w)} neutrons (useful: {n_useful})')
        hist_lambda.fill( wl[mask_fwd], w[mask_fwd] )
        #If we return True, source will stop producing more neutrons:
        n_useful_tot[0] += n_useful
        return n_useful_tot[0] > 500.0

    #Run simulation (note we set a high value of n since we halt the src
    #manually):
    res = mmcrun( cfgstr = "stdlib::Polyethylene_CH2.ncmat;temp=20",
                  geomcfg = "slab;dz=0.001",
                  srccfg = "constant;n=1e18;z=-0.05;ekin=thermal:300K",
                  callback = cb,
                  #Just to test these options somewhere:
                  enginecfg = 'nthreads=2;roulette=0.1,1.0,0',
                  callback_options='basket_type=basic;ncaches=1;cachelen=1e5'
                 )

    #Summarise and plot:
    def stats(key):
        return ( res.output_metadata[key]['count'],
                 res.output_metadata[key]['weight'] )
    print('Neutron counts:')
    print("   Produced by source:  %i (flux=%.8g)"%stats('provided'))
    print("   Missing geometry:    %i (flux=%.8g)"%stats('miss'))
    print("   Available for tally: %i (flux=%.8g)"%stats('tallied'))
    print('   Considered "useful": %i'%n_useful_tot[0])
    if do_plot:
        hist_lambda.to_hist1d().plot(logy=True)

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:])
