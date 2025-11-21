#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
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

import NCrystalDev._mmc as ncmmc
from NCrystalDev._mmc_results import MMCDiffractionResults as Results
from NCTestUtils.env import ncsetenv
import NCTestUtils.dirs as dirs

def main(do_plot, do_update):
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')

    if do_update:
        datadir = dirs._test_data_dir_for_updates()
        assert datadir, "can only update in ncrystal dev repo"
    else:
        datadir = dirs.test_data_dir

    reffileA = datadir.joinpath('mmcfmwkA.json')
    reffileB = datadir.joinpath('mmcfmwkB.json')

    resA = ncmmc.runsim_diffraction_pattern(
        cfgstr='Al_sg225.ncmat;comp=bragg;dcutoff=0.8',
        srccfg = 'constant;wl=1.8;z=-0.009999999;n=1e4',
        geomcfg = 'sphere;r=0.1',
        nthreads = 1
    ).rebin(50)
    resB = ncmmc.quick_diffraction_pattern( 'Al_sg225.ncmat;comp=inelas',
                                            neutron_energy = "2.0Aa",
                                            material_thickness = "0.1mfp",
                                            nthreads = 2,
                                            nstat = 1e4).rebin(50).clone()
    #Basic serialisation/deserialisation check:
    j = resA.to_json()
    d = resA.to_dict()
    res_j = Results.from_json(j)
    res_d = Results.from_dict(d)
    n = len(resA.histograms)
    assert len(res_j.histograms) == n
    assert len(res_d.histograms) == n
    for i in range(n):
        for r in [res_j,res_d]:
            resA.histograms[i].check_compat(r.histograms[i],
                                            threshold = 1.0,
                                            check=True)

    #Rely on reference files, to check stability and to have stable output for
    #further testing:
    if do_update:
        reffileA.write_text(resA.to_json())
        print(f"Updated {reffileA}")
        reffileB.write_text(resB.to_json())
        print(f"Updated {reffileB}")
        print("Aborting due to update")
        return
    #load refs and check compatibility:
    resA_ref = Results.from_json(reffileA.read_text())
    resB_ref = Results.from_json(reffileB.read_text())

    def check_compat(h1,h2):
        h1.check_compat( h2, check=True, threshold = 0.001 )

    print("Verifying resA stability.")
    check_compat(resA.histogram_main,resA_ref.histogram_main)
    print("Verifying resB stability.")
    check_compat(resB.histogram_main,resB_ref.histogram_main)
    #Proceed with ref data, for stable values in printouts:
    print("Done. Proceed with ref data.")
    resA = resA_ref
    resB = resB_ref

    def plotcmp( h, href, title ):
        href.plot(
            do_show=False,error_bands=1.0,
            alpha=0.3,color='blue',label='ref'
        )
        plt=h.plot(
            do_show=False,color='none',logy=True,label='new'
        )
        plt.title(title)
        plt.legend()
        plt.grid()
        plt.show()
    plotcmp( resA.histogram_main, resA_ref.histogram_main, 'A' )
    plotcmp( resB.histogram_main, resB_ref.histogram_main, 'B' )

    resA.plot_breakdown(rebin_factor=2,logy=True)
    resB.plot_breakdown(rebin_factor=2,logy=True)

    resA.plot(rebin_factor=2,logy=False)
    resA.histogram_main.clone(rebin_factor=2).dump()
    print("resA.histogram_main.title=%s"%repr(resA.histogram_main.title))
    print(resA.histogram_titles)
    h = resA.histogram_sum(select=['NOSCAT','MULTISCAT_PUREELAS'],
                           exclude='SINGLESCAT_ELAS')
    h.clone(rebin_factor=2).plot(logy=True)
    h = resA.histogram_sum(select='NOSCAT')
    assert h is resA.histogram('NOSCAT')
    assert resA.histogram_sum() is resA.histogram('MAIN')
    assert resA.histogram_sum() is resA.histogram_main
    print( [k for k in sorted(resA.histogram_breakdown.keys())] )


    print(ncmmc.quick_diffraction_pattern_autoparams('Be_sg194.ncmat'))
    for x in [0.0,1e-12, 1e-10, 2.123456789e-10, 1e-6, 2.123456789e-6, 1e-3,
              2.123456789e-3, 1e-2, 2.123456789e-2, 1e-1, 2.123456789e-2,
              1.0, 2.1234567893, 1004343.0 ]:
        for r2d in [True,False]:
            s=ncmmc._encode_length_to_str( x, round2digits=r2d )
            print( 'format length=%.14g (m) round2=%i -> %s'%(x,r2d,repr(s)))

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:],
         do_update = '--update' in sys.argv[1:])
