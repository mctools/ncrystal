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
from NCrystalDev.minimc_objects import MMCResults, MMCTallyView
from NCrystalDev.hist import Hist1D
from NCTestUtils.env import ncsetenv
import NCTestUtils.dirs as dirs
from NCrystalDev.exceptions import NCBadInput, NCCalcError
from NCTestUtils.common import ensure_error
import pprint

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
    resB = MMCResults(resB)

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
    resA_ref2 = ncmmc.MMCResults(reffileA.read_text())
    resB_ref = ncmmc.MMCResults(reffileB)
    print("Testing equality")
    assert not (resA_ref==117)
    assert (resA_ref==resA_ref)
    assert (resA_ref==resA_ref2)
    assert (resA_ref!=resB_ref)
    assert not (resA_ref==resB_ref)
    assert resA_ref.tally('theta')==resA_ref.tally('theta')
    assert resA_ref.tally('theta')==resA_ref2.tally('theta')
    assert resA_ref.tally('theta') is not None
    assert resA_ref.tally('theta')!=117
    assert resB_ref.tally('theta')!=resA_ref2.tally('theta')

    #Fake ripout the breakdown hists:
    assert resA_ref == resA_ref2
    d=resA_ref2.tally('theta')._raw_data()
    d['breakdown']=None
    assert resA_ref != resA_ref2

    print("Test plotting code")
    def plotcmp( h, href, title ):
        href.plot(
            do_show=False,error_bands=1.0,
            alpha=0.3,color='blue',label='ref',
            title=False
        )
        plt=h.plot(
            do_show=False,color='none',logy=True,label='new',title=False
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
    #Examples from doc-strings:
    tt = resA.tally('theta')
    tt.histogram_sum(select=['NOSCAT','MULTISCAT_PUREELAS'])
    tt.histogram_sum(exclude='SINGLESCAT_ELAS')

    h.clone(rebin_factor=2).plot(logy=True)
    h = tt.histogram_sum(select='NOSCAT')
    assert h == tt.histograms['NOSCAT']
    assert h is not tt.histograms['NOSCAT']#should be cloned
    assert tt.histogram_sum() == tt.hist_total
    assert tt.histogram_sum() is not tt.hist_total
    print( sorted(tt.hist_breakdown.keys()) )

    ncmmc.gen_doc( 'engine' )
    ncmmc.gen_doc( 'src' )
    ncmmc.gen_doc( 'geom' )


    def pdc( cfgstr, cfgtype ):
        print(f'==> Decoding {cfgtype} "{cfgstr}":')
        d = ncmmc.decode_cfgstr( cfgstr, cfgtype )
        pprint.pp(d)

    pdc('circular;r=0.3;uy=2;uz=0;x=17;ekin=0.025+-0.001', 'src')
    pdc('box;dx=1;dy=1e-3;dz=0.0123', 'geom')
    pdc('tallybins=+;tally=q,de', 'engine')
    with ensure_error(NCBadInput,'cfgstr parameter must be a string'):
        pdc(0.025,'engine')
    with ensure_error(NCBadInput,
                      'cfgtype parameter must be "src", "geom", or "engine"'):
        pdc('foo',0.025)
    with ensure_error(NCBadInput,
                      'cfgtype parameter must be "src", "geom", or "engine"'):
        pdc('foo','bar')


    with ensure_error(NCBadInput,'Data seems to be in an unsupported format'):
        MMCResults(dict( datatype = 'NCrystalMiniMCResults_v17',
                         input = {}, output={} ))
    with ensure_error(NCBadInput,'Data seems to be in an unsupported format'):
        MMCResults(dict( datatype = 'NCrystalMiniMCResults_v1',
                         input = {}, output={}, somethingnew={} ))
    with ensure_error(NCBadInput,'Unsupported data format'):
        MMCResults([1,2,3])
    resV3 = mmc_run( 'void.ncmat',scenario = '1eV on 1cm 1 times',
                     enginecfg = 'nthreads=1;tally=q,de,mu' )
    resV0 = mmc_run( 'void.ncmat',scenario = '1eV on 1cm 1 times',
                     enginecfg = 'nthreads=1;tally=' )
    with ensure_error(NCBadInput,'Tally not available in MiniMC dataset: "e"'
                      ' (available tallies are "de", "mu", "q")'):
        resV3.tally('e')
    with ensure_error(NCBadInput,'Tally not available in MiniMC dataset: "e"'
                      ' (no tallies were enabled!).'):
        resV0.tally('e')
    print('resA long title:')
    print(repr(resA.long_title()))
    print('resB long title:')
    print(repr(resB.long_title()))
    print('resV3 long title:')
    print(repr(resV3.long_title()))
    print('resV0 long title:')
    print(repr(resV0.long_title()))
    with ensure_error(NCCalcError,
                      'Incompatible MMCResults (incompatible input'
                      ' values for "material/cfgstr")'):
        resA.check_compat( resB, check = True )

    with ensure_error(TypeError,'Do not create MMCTallyView objects directly'):
        MMCTallyView( {'foo':'bar'} )
    assert resA.tally('theta').mother is resA
    assert not ( resA == resB )
    assert not ( resV0 == resV3 )
    print("some tally units and short descriptions:")
    for t in resV3.tallies:
        print(repr(t.name),repr(t.unit),repr(t.short_description))

    print("resV3 mu tally dump method:")
    resV3.tally('mu').dump(contents=False)

    from NCrystalDev._mmc_impl import _determine_rebin_factor as drf
    assert drf( current_nbins=100, max_nbins= 200 ) == 1
    assert drf( current_nbins=100, max_nbins= 20 ) == 5
    assert drf( current_nbins=100 ) == 1
    assert drf( current_nbins=100, rebin_factor=5 ) == 5
    assert drf( current_nbins=100, max_nbins=100 ) == 1
    with ensure_error(NCBadInput,
                      'Invalid rebin factor 17 is not a divisor of nbins=100.'):
        drf( current_nbins=100, rebin_factor=17 )
    with ensure_error(NCBadInput,
                      'Can not set both max_nbins and rebin_factor.'):
        drf( current_nbins=100, max_nbins=20, rebin_factor=5 )

    def runerr( errmsg, *a, **kw ):
        with ensure_error(NCBadInput,errmsg):
            mmc_run(*a,**kw)
    runerr( 'Missing required parameter: cfgstr.', None )
    runerr( 'The cfgstr parameter must be a string.',
            1.8 )
    runerr( 'The scenario parameter must be a string.',
            'void.ncmat', scenario=1.8 )
    runerr( 'Invalid enginecfg tallybins entry "e:5:0.0".',
            'void.ncmat', scenario='', enginecfg='tallybins=e:5:0.0,0.02' )
    runerr( 'The srccfg parameter must be a string.',
            'void.ncmat', srccfg=b'1.8' )
    runerr( 'The geomcfg parameter must be a string.',
            'void.ncmat', geomcfg=1.8 )
    runerr( 'The callback_options parameter must be a string.',
            'void.ncmat', callback_options=1.8 )
    runerr( 'Inconsistent parameters. Do not supply callback_options'
            ' without a callback function.',
            'void.ncmat', scenario='', callback_options='' )
    mmc_run('void.ncmat')#just a cfgstring => scenario string is an empty string
    runerr( 'Inconsistent parameters. Do not supply geomcfg or'
            ' srccfg when also supplying a scenario string.',
            'void.ncmat', scenario='2Aa',srccfg='constant;wl=1.8')
    #runerr( 'Missing geomcfg parameter.', 'void.ncmat',srccfg='constant;wl=1.8')
    #runerr( 'Missing srccfg parameter.', 'void.ncmat',geomcfg='sphere;r=1')

    kw = dict(cfgstr='void.ncmat',scenario = '0.01eV on 1cm 1 times',
              enginecfg
              = 'nthreads=1;tally=e;tallybins=e:5:0:02;tallybreakdown=0')
    r = mmc_run( **kw, unpack='dict' )
    h = r['output']['tally']['e']['total']
    assert isinstance(h,Hist1D)
    r['output']['tally']['e']['total'] = h.to_json()
    pprint.pp(r)
    r = mmc_run( **kw, unpack='dict_jsoncompat' )
    pprint.pp(r)
    r = mmc_run( **kw, unpack='json' )
    print(r)
    runerr( 'Invalid value of unpack (must be "dict", "json",'
            ' "dict_jsoncompat", or "object"): '"'foobar'",
            **kw, unpack='foobar' )
    runerr( 'Missing geomcfg parameter.',
            'void.ncmat', srccfg='constant;wl=1.8' )
    assert not MMCResults(r).check_compat(resA)

    #Test plot without breakdown:
    assert 'tallybreakdown=0' in kw['enginecfg']
    r = mmc_run(**kw)
    r.tally('e').plot()
    r.tally('e').plot(logy=None,title=False)
    r.tally('e').plot(title='my title')

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:],
         do_update = '--update' in sys.argv[1:])
