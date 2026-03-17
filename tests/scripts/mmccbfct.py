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

import NCrystalDev.minimc as ncmmc
import NCrystalDev.core as nccore

def main(do_plot):
    #FIXME: This is not testing much. We should monitor the resulting histograms
    #for compatibility with a reference.

    from NCrystalDev.hist import HistFiller1D as Hist

    hist_e0_fwd = Hist( 500, 0.0, 0.1, title='e0 (fwd)' )
    hist_e0_back = Hist( 500, 0.0, 0.1, title='e0 (back)' )
    hist_e_fwd = Hist( 400, 0.0, 0.1, title='E (fwd)' )
    hist_e_back = Hist( 400, 0.0, 0.1, title='E (back)' )
    hist_nscat_fwd = Hist( 22, -1.5, 20.5, title='nscat (fwd)' )
    hist_nscat_back = Hist( 22, -1.5, 20.5, title='nscat (back)' )
    hist_nscat = Hist( 22, -1.5, 20.5, title='nscat' )


    ntot = [0]
    def cb( data ):
        nscat,e,w,uz,e0 = (data['nscat'],data['ekin'], data['w'], data['uz'],
                           data['ekin0'])
        print('Callback processing %i neutrons'%len(w))
        print('    Available fields:',' '.join(data.keys()))
        if False:
            print('   -> nscat range: %g to %g'%( nscat.min(), nscat.max() ))
            print('   -> w range: %g to %g'%( w.min(), w.max() ))
            print('   -> uz range: %g to %g'%( uz.min(), uz.max() ))
            print('     -> counts: %i (nscat=-1), %i (nscat=0), %i (nscat=1), %i (nscat>1)'%(len(e[nscat==-1]),
                                                                                             len(e[nscat==0]),
                                                                                             len(e[nscat==1]),
                                                                                             len(e[nscat>1]),))
            print('     -> uz @ nscat=0:',uz[nscat==0])
            print('     -> w @ nscat=0:',w[nscat==0])
            print('     -> uz @ nscat=1:',uz[nscat==1])
            print('     -> w @ nscat=1:',w[nscat==1])
        ntot[0] += len( w )
        mask_fwd = uz > 0.5
        mask_back = uz < -0.5
        #hist_x0_fwd.fill( x0[mask_fwd], w[mask_fwd] )
        #hist_x0_back.fill( x0[mask_back], w[mask_back] )
        hist_e0_fwd.fill( e0[mask_fwd], w[mask_fwd] )
        hist_e0_back.fill( e0[mask_back], w[mask_back] )
        hist_e_fwd.fill( e[mask_fwd], w[mask_fwd] )
        hist_e_back.fill( e[mask_back], w[mask_back] )
        hist_nscat.fill( nscat )
        hist_nscat_fwd.fill( nscat[mask_fwd], w[mask_fwd] )
        hist_nscat_back.fill( nscat[mask_back], w[mask_back] )

    nccore.enableFactoryThreads()
    res = ncmmc.run(
        #'void.ncmat',
        'Al_sg225.ncmat;temp=300',
        scenario='6Aa on 10cm slab 1e6 times',
        enginecfg='nthreads=auto',#;absorption=0',#;nscatlimit=1',#fixme: allow nthreads=0?
        callback = cb
    )
    print(res.setup['geom']['decoded'])
    print(res.setup['src']['decoded'])
    tallied_stats = res.output_metadata['tallied']

    #FIXME put back in: assert tallied_stats['count'] == ntot[0]

    print('Total neutrons tallied (count): %i'%tallied_stats['count'])
    print('Total neutrons tallied (sumw): %i'%tallied_stats['weight'])
    print('Total neutrons tallied in E(fwd) hist (sumw) %i'
          % hist_e_fwd.to_hist1d().contents.sum())
    print('Total neutrons tallied in E(back) hist (sumw) %i'
          % hist_e_back.to_hist1d().contents.sum())
    #res.dump()
    #res.tally('theta').dump()

    if not do_plot:
        return
    #res.tally('theta').plot()

    #hist_nscat.to_hist1d().plot()
    #raise SystemExit

    common = dict( error_bands=2.0, alpha = 0.5, do_show=False)

    def pl( hfwd, hback ):
        hfwd.to_hist1d().plot(label='fwd',**common,color='red')
        plt = hback.to_hist1d().plot(label='back',**common,color='green')
        plt.gca().legend()
        plt.gca().semilogy()
        plt.grid()
        plt.show()
    pl(hist_e0_fwd,hist_e0_back)
    #pl(hist_x0_fwd,hist_x0_back)
    pl(hist_e_fwd,hist_e_back)
    pl(hist_nscat_fwd,hist_nscat_back)

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:])
