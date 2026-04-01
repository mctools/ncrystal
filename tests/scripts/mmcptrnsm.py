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
from NCTestUtils.dirs import get_named_test_data_dir
import NCrystalDev.core as nccore
from NCrystalDev.hist import Hist1D
from NCrystalDev.plot import PlotContext
import math
import numpy as np

def get_Al_noabs():
    from NCrystalDev.ncmat import NCMATComposer
    from NCrystalDev.atomdata import atomDB
    c, a = NCMATComposer('stdlib::Al_sg225.ncmat'), atomDB('Al')
    c.update_atomdb( 'Al',
                     mass = a.averageMassAMU(),
                     coh_scat_len = a.coherentScatLenFM(),
                     incoh_xs = a.incoherentXS(),
                     abs_xs = 0.0 )
    n = 'Al_noabsn.ncmat'
    c.register_as(n)
    return f'virtual::{n}'

def main(do_plot, do_update):
    from NCrystalDev.hist import HistFiller1D as Hist
    nccore.enableFactoryThreads(2)

    #Define parameters and histograms:
    cfgstr = 'Al_sg225.ncmat'
    #cfgstr = get_Al_noabs() #uncomment to investigate Al with no absorption
    transm_def_degree = 20
    thickness_cm = 10
    elow, ehigh = 0.003, 0.0055
    nbins = 100
    n = nbins*4000

    hist_e0_transm = Hist( nbins, elow, ehigh, title='ekin_initial (transm)' )

    #Define data processing callback function which fills the histograms:
    transm_uz_threshold = math.cos( transm_def_degree*math.pi/180)

    def cb( data ):
        print('Callback processing %i neutrons'%len(data['w']))
        print('   got keys:',data.keys())
        w,uz,e0 = data['w'], data['uz'], data['ekin0']
        mask_transm = uz > transm_uz_threshold
        hist_e0_transm.fill( e0[mask_transm], w[mask_transm] )

    #Run simulation:
    dt = thickness_cm*0.01/2#half-thickness in meter
    res = mmcrun( cfgstr,
                  geomcfg = f"slab;dz={dt}",
                  srccfg = f"constant;z=-{dt};n={n};ekin={elow}-{ehigh}",
                  enginecfg="nthreads=3;tally=nscat,nscat_uw",
                  callback = cb )

    tallied_stats = res.output_metadata['tallied']
    print('Total neutrons tallied (count): %i'%tallied_stats['count'])
    print('Total neutrons tallied (sumw): %i'%tallied_stats['weight'])

    h = hist_e0_transm.to_hist1d()
    #normalise to expected sum(weight)/bin in absence of interactions:
    h.scale(100.0*nbins/n)

    #Find reference hist:
    td = get_named_test_data_dir('mmcref',for_updates = do_update)
    f_href = td.joinpath('href_ptrnsm.json')
    if do_update:
        f_href.write_text(h.to_json())
        print(f'Updated {f_href}. Aborting.')
        return
    if not f_href.is_file():
        raise SystemExit(f'Reference file {f_href} not found.'
                         ' Run with --update to generate')
    h_ref = Hist1D(f_href.read_text())
    pval = h.check_compat( h_ref, return_pval = True )
    print( "Pvalue for comp. with ref"
           f" (higher is more compatible): {pval:g}" )

    def result():
        if pval < 0.01:
            raise SystemExit('Incompatibility detected')


    if not do_plot:
        return result()

    pctx=PlotContext()
    h_ref.plot(error_bands=1.0,alpha=0.3,color='blue',label='ref',
               **pctx.kwargs_subcontext())
    h.plot(color='none',logy=True,label='new',**pctx.kwargs_subcontext())
    pctx.finalise(do_legend=True,do_grid=True)
    pctx=PlotContext()
    h.plot(error_bands=1.0,label='MiniMC (includes rescattering)',
           color='blue',alpha=0.5,**pctx.kwargs_subcontext())
    h_ref.plot(error_bands=1.0,label='Reference',color='green',
               alpha=0.5,**pctx.kwargs_subcontext())
    e = np.linspace(elow,ehigh,2000)
    m = nccore.load(cfgstr)
    def ptransm(xs):
        return 100.0*np.exp(-thickness_cm*m.info.factor_macroscopic_xs*xs)
    xs_abs =  m.absorption.xsect(e)
    xs_scat = m.scatter.xsect(e)
    pctx.finalise(do_legend=True,do_grid=True)


    pctx.axis.plot(e,ptransm(xs_abs),label='Theory (no rescattering, absorption only)',color='orange')
    pctx.axis.plot(e,ptransm(xs_abs+xs_scat),label='Theory (no rescattering)',color='red')
    pctx.axis.plot(e,ptransm(xs_scat),label='Theory (no rescattering, scattering only)',color='brown')
    pctx.axis.set_title(f'Ptransmission (theta_scat<{transm_def_degree}degree)')
    pctx.axis.set_ylabel('(%)')
    pctx.axis.set_xlabel('E (eV)')
    pctx.axis.set_ylim(0,100.0)
    pctx.axis.get_figure().tight_layout()
    pctx.finalise(do_legend=True,do_grid=True)

    res.tally('nscat').plot(logy=False)
    res.tally('nscat_uw').plot(logy=False)
    return result()

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:],
         do_update = '--update' in sys.argv[1:])
