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
from NCTestUtils.env import ncsetenv
import NCTestUtils.dirs as dirs

def main(do_plot, do_update):
    mmc_run = ncmmc.run

    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')

    if do_update:
        datadir = dirs._test_data_dir_for_updates()
        assert datadir, "can only update in ncrystal dev repo"
    else:
        datadir = dirs.test_data_dir

    reffileA = datadir.joinpath('mmcfmwkA.json')
    reffileB = datadir.joinpath('mmcfmwkB.json')


    print("Run two MiniMC scenarios")

    resA = mmc_run(
        cfgstr='Al_sg225.ncmat;comp=bragg;dcutoff=0.8',
        srccfg = 'constant;wl=1.8;z=-0.009999999;n=1e4',
        geomcfg = 'sphere;r=0.1',
        enginecfg='nthreads=1;tally=theta;tallybins=theta:36:0:180'
    )
    resB = mmc_run( 'Al_sg225.ncmat;comp=inelas',
                    scenario = '2.0Aa on 0.1mfp 1e4 times',
                    enginecfg = (';nthreads=2;tally=theta'
                                 ';tallybins=theta:36:0:180') )

    print("Basic serialisation/deserialisation check")
    j = resA.to_json()
    d = resA.to_dict()
    res_j = ncmmc.MMCResults(j)
    res_d = ncmmc.MMCResults(d)

    def assert_dict_same(d1,d2,_dir=''):
        assert isinstance(d1,dict)
        assert isinstance(d2,dict)
        if not d1.keys()==d2.keys():
            raise RuntimeError('Incompatible keys at lvl'
                               '="%s": %s vs %s'%(_dir,
                                                  d1.keys(),
                                                  d2.keys()))
        for k in d1.keys():
            if d1[k] != d2[k]:
                keyprint = k
                if _dir:
                    keyprint=_dir+'/'+keyprint
                print('Issues with  key=%s'%keyprint)
                if isinstance(d1[k],dict):
                    assert_dict_same(d1[k],d2[k],_dir=_dir+'/%s'%k)
                else:
                    print(type(d1[k]))
                    print(type(d2[k]))
                    print('equal:',d1[k]==d2[k])
                raise RuntimeError('Issues with key=%s'%keyprint)
        assert d1 == d2

    assert_dict_same(res_j.to_dict(), d)
    assert_dict_same(res_d.to_dict(), d)
    assert res_j.to_dict() == d
    assert res_d.to_dict() == d
    nt = len(resA.tallies)
    assert len(res_j.tallies) == nt
    assert len(res_d.tallies) == nt
    tns = resA.tally_names
    assert len(tns)==nt
    for tn in tns:
        for r in [res_j,res_d]:
            t1 = resA.tally(tn)
            t2 = r.tally(tn)
            h1 = t1.histograms
            h2 = t2.histograms
            assert h1.keys()==h2.keys()
            assert 'SINGLESCAT_ELAS' in h1
            assert 'total' in h1
            for hk in h1.keys():
                assert h1[hk] == h2[hk]

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

    print("Loading reference files")
    resA_ref = ncmmc.MMCResults(reffileA.read_text())
    resB_ref = ncmmc.MMCResults(reffileB)

    print("Test plotting code")
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
    plotcmp( resA.tally('theta').hist_total,
             resA_ref.tally('theta').hist_total, 'A' )
    if do_plot:
        #resB fluctuates too much since nthreads=2, so should not leave plot
        #curves in reflogs.
        plotcmp( resB.tally('theta').hist_total,
                 resB_ref.tally('theta').hist_total, 'B' )

    #Note that Since A is generated with nthreads=1, we can in principle expect
    #the same results if running again on the same machine. However, there are
    #potential FP issues if changing platform, and B has nthreads=2 which ruins
    #it even on the same platform. Thus, we will only check for statistical
    #compatiblility with the reference files (todo: so we should not plotcmp
    #resA above??).

    print("Verifying resA conversion stability.")
    assert resA.check_compat( resA, threshold=1.0 ) is True
    resA_ref2 = ncmmc.MMCResults(reffileA.read_bytes())
    resA_ref3 = ncmmc.MMCResults(reffileA)
    assert resA_ref._raw_data() == resA_ref2._raw_data()
    assert resA_ref._raw_data() == resA_ref3._raw_data()

    #Next check could in principle fail (but unlikely):
    resA_ref.check_compat( resA, threshold=0.05, check=True )

    print("Verifying resB conversion stability.")
    assert resB.check_compat( resB, threshold=1.0 ) is True
    assert resB_ref.check_compat( resB, threshold=0.01 ) is True

    print("Verifying extreme check_compat.")
    assert resA.check_compat( resA ) is True
    assert resA.check_compat( resB, threshold=1.0 ) is False
    assert resA.check_compat( resB, threshold=0.0 ) is False

    print("Proceed with ref data only.")
    resA = resA_ref
    resB = resB_ref
    del resA_ref
    del resA_ref2
    del resA_ref3
    del resB_ref

    print("More plotting code test")
    resA.tally('theta').plot(rebin_factor=2,logy=True)
    if do_plot:
        #Again, resB with nthreads=2 give irreproducible plot curves:
        resB.tally('theta').plot(rebin_factor=2,logy=True)
    resA.tally('theta').plot(rebin_factor=2,logy=False)
    print("Dump test")
    resA.tally('theta').hist_total.clone(rebin_factor=2).dump()
    print("resA.tally('theta').hist_total.title=%s"%repr(resA.tally('theta').hist_total.title))
    print(resA.tally_names)
    print("Sum test")
    h = resA.tally('theta').histogram_sum(select=['NOSCAT','MULTISCAT_PUREELAS'],
                                       exclude='SINGLESCAT_ELAS')
    h.clone(rebin_factor=2).plot(logy=True)
    h = resA.tally('theta').histogram_sum(select='NOSCAT')
    assert h is resA.tally('theta').histograms['NOSCAT']
    assert resA.tally('theta').histogram_sum() is resA.tally('theta').hist_total
    print( sorted(resA.tally('theta').hist_breakdown.keys()) )

    ncmmc.gen_doc( 'engine' )
    ncmmc.gen_doc( 'src' )
    ncmmc.gen_doc( 'geom' )

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:],
         do_update = '--update' in sys.argv[1:])
