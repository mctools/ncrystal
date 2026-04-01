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
from NCrystalDev.constants import ekin2wl
from NCTestUtils.dirs import get_named_test_data_dir
from NCrystalDev.hist import Hist1D
from NCrystalDev.plot import PlotContext

def main(do_plot, do_update):
    from NCrystalDev.hist import HistFiller1D

    hists = {}
    def add_hist( key, nbins, xmin, xmax ):
        assert key not in hists
        hists[key] = HistFiller1D(nbins,xmin,xmax,title=key)
    add_hist( 'wl0_fwd', 100, 2.5, 5.5 )
    add_hist( 'wl0_back', 100, 2.5, 5.5 )
    add_hist( 'e_fwd', 100, 0.0, 0.1 )
    add_hist( 'e_back', 100, 0.0, 0.1 )
    add_hist( 'nscat_fwd', 22, -1.5, 20.5 )
    add_hist( 'nscat_back', 22, -1.5, 20.5 )
    #add_hist( 'nscat', 22, -1.5, 20.5 )

    ntot = [0]
    def cb( data ):
        nscat,e,w,uz,wl0 = (data['nscat'],data['ekin'], data['w'], data['uz'],
                            ekin2wl(data['ekin0']))
        print('Callback processing %i neutrons'%len(w))
        print('    Available fields:',' '.join(data.keys()))
        assert len(data)==18, "callback should trigger extended baskets"
        ntot[0] += len( w )
        mask_fwd = uz > 0.5
        mask_back = uz < -0.5
        hists['wl0_fwd'].fill( wl0[mask_fwd], w[mask_fwd] )
        hists['wl0_back'].fill( wl0[mask_back], w[mask_back] )
        hists['e_fwd'].fill( e[mask_fwd], w[mask_fwd] )
        hists['e_back'].fill( e[mask_back], w[mask_back] )
        hists['nscat_fwd'].fill( nscat[mask_fwd], w[mask_fwd] )
        hists['nscat_back'].fill( nscat[mask_back], w[mask_back] )

    nccore.enableFactoryThreads()
    res = ncmmc.run(
        'Al_sg225.ncmat;temp=300',
        geomcfg="slab;dz=0.05",
        srccfg="constant;n=6e5;z=-0.05;wl=4+-0.2",
        enginecfg='nthreads=2',
        callback = cb
    )
    assert ( res.setup['geom']['decoded']['short_description']
             == 'slab with thickness 100mm' )
    tallied_stats = res.output_metadata['tallied']

    assert tallied_stats['count'] == ntot[0]

    #convert to Hist1D:
    for k in hists.keys():
        hists[k] = hists[k].to_hist1d()

    print()
    print('Total neutrons tallied (count): %i'%tallied_stats['count'])
    print('Total neutrons tallied (sumw): %i'%tallied_stats['weight'])
    print('Total neutrons tallied in E(fwd) hist (sumw) %i'
          % hists['e_fwd'].contents.sum())
    print('Total neutrons tallied in E(back) hist (sumw) %i'
          % hists['e_back'].contents.sum())


    #Find refs:
    td = get_named_test_data_dir('mmcref',for_updates = do_update)
    def key2path(key):
        return td.joinpath(f'href_cbfct_{key}.json')

    #update refs:
    if do_update:
        for key, h in sorted(hists.items()):
            f_href = key2path(key)
            print(f"Updating {f_href}")
            f_href.write_text(h.to_json())
        raise SystemExit('Updates done. Aborting.')


    hists_ref = {}
    for key in sorted(hists.keys()):
        f_href = key2path(key)
        if not f_href.is_file():
            raise SystemExit(f'Reference file {f_href} not found.'
                             ' Run with --update to generate')
        hists_ref[key] = Hist1D(f_href.read_text())

    worst_pval = 1.0
    print()
    for key in sorted(hists.keys()):
        h, h_ref = hists[key], hists_ref[key]
        pval = h.check_compat( h_ref, return_pval = True )
        print( f"{key}: Pvalue for comp. with ref"
               f" (higher is more compatible): {pval:g}" )
        worst_pval = min(worst_pval,pval)

    def result():
        if worst_pval < 0.01:
            raise SystemExit('Incompatibility detected')

    if not do_plot:
        return result()

    common = dict( error_bands=2.0, alpha = 0.5)

    def pl( hfwd_key, hback_key ):
        hfwd = hists[hfwd_key]
        hback = hists[hback_key]
        pctx=PlotContext()
        hfwd.plot(label=hfwd.title,**common,color='red',
                  **pctx.kwargs_subcontext())
        hback.plot(label=hback.title,**common,color='green',
                   **pctx.kwargs_subcontext())
        pctx.axis.semilogy()
        pctx.finalise(do_legend=True,do_grid=True)
    pl('wl0_fwd','wl0_back')
    pl('e_fwd','e_back')
    pl('nscat_fwd','nscat_back')

    for key in sorted(hists.keys()):
        h, h_ref = hists[key], hists_ref[key]
        pctx=PlotContext()
        h_ref.plot(error_bands=1.0,alpha=0.3,color='blue',label='ref',
                   **pctx.kwargs_subcontext())
        h.plot(color='none',logy=True,label='new',
               **pctx.kwargs_subcontext())
        pctx.axis.set_title(key)
        pctx.finalise(do_legend=True,do_grid=True)

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:],
         do_update = '--update' in sys.argv[1:])
