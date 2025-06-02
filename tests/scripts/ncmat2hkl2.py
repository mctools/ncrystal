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

from NCTestUtils.common import print_text_file_with_snipping
import NCrystalDev as NC
import NCrystalDev.cli as nccli
import numpy as np
import pathlib
import sys
do_plot = '--plot' in sys.argv[1:]
if do_plot:
    import matplotlib.pyplot as plt #noqa E402


#default is to encode numbers in produced .laz/.lau files with 14 digits of
#precision (%.14g), but for robustness of the unit test we reduce this in any
#printed content to avoid spurious false positives due to FP instabilities:
test_precision = 10

max_diff_lvl = 1e-7

def cfgstr2hkl(cfgstr,fmt):
    assert fmt in ('lau','laz')
    o = pathlib.Path('./tmp_redirected_output.txt')
    if o.is_file():
        o.unlink()
    nccli.run( 'ncmat2hkl',
               f'--override-prec={test_precision}',#hidden unit test option
               f'--format={fmt}',
               cfgstr,
               '-o', o.name
              )
    v = NC.__version__
    return o.read_text().replace(f'# File created by NCrystal v{v}',
                                 '# File created by NCrystal v<current>')

itest = [0]
def test_cfgstr(cfgstr):
    #NC.registerInMemoryFileData('foo.laz',cfgstr2hkl(cfgstr,'laz'))
    def bragg_only(c):
        return c+';inelas=0;elas=0;bragg=1'
    def make_oriented(c):
        return c + (';dir1=@crys_hkl:0,0,1@lab:0,0,1;'
                    'dir2=@crys_hkl:0,1,0@lab:0,1,0;mos=10deg;dirtol=180deg')
    sc_orig = NC.createScatter(bragg_only(cfgstr))
    sc_orig_oriented = NC.createScatter(make_oriented(bragg_only(cfgstr)))
    hr="="*100+'\n'
    def prtitle(title):
        print('\n'+hr*3 + f"====> {title} <=====\n" +hr*3+'\n')
    for fmt in ('laz','lau'):
        if fmt=='lau':
            itest[0] += 1
            if itest[0] % 3 != 0:
                continue
        fn=f'foo.{fmt}'

        prtitle(f'NCrystal dump of cfgstr={cfgstr}')
        info=NC.createInfo(cfgstr)
        info.dump(verbose=2)
        wls = np.linspace(0.001,info.braggthreshold*1.2,20000)
        wls_sparse = np.linspace(wls[0],wls[-1],2000)

        prtitle(f'Conversion of {cfgstr} to {fmt}')
        NC.registerInMemoryFileData(fn,cfgstr2hkl(cfgstr,fmt))
        rd=NC.createTextData(fn).rawData
        assert rd.endswith('\n')
        print_text_file_with_snipping(rd,nstart=100,nend=40)

        prtitle(f'NCrystal dump of converted {fmt}')
        print_text_file_with_snipping(NC.createInfo(fn).dump_str(verbose=2),
                                      nstart=70,nend=20)

        prtitle(f'Testing bragg component of converted {fmt} vs. original')
        sc=NC.createScatter(bragg_only(fn))
        sc.dump()
        xs=sc.xsect(wl=wls)
        xs_orig=sc_orig.xsect(wl=wls)
        diff=abs(xs-xs_orig).max()
        if do_plot:
            print (f'XS max-abs-diff converted {fmt} vs. original:',diff)
        else:
            print (f'XS max-abs-diff converted {fmt} vs. original < '
                   f'{max_diff_lvl}?: %s'%('yes' if diff<max_diff_lvl else 'no'))
        if do_plot:
            plt.plot(wls/2,xs,label='via laz')
            plt.plot(wls/2,xs_orig,label='orig')
            plt.legend()
            plt.show()
        if not diff < max_diff_lvl:
            raise SystemExit(f'maxdiff test failed with diff: {diff}')

        #Test that can be used for single crystal, with huge mosaicity so we can
        #get non-trivial cross sections without careful aiming:
        sc=NC.createScatter(make_oriented(bragg_only(fn)))
        sc.dump()
        xs = sc.xsect(wl=wls_sparse, direction=(0,0,1))
        xs_orig = sc_orig_oriented.xsect(wl=wls_sparse, direction=(0,0,1))
        diff=abs(xs-xs_orig).max()
        if do_plot:
            print (f'SC XS max-abs-diff converted {fmt} vs. original:',diff)
        else:
            print (f'SC XS max-abs-diff converted {fmt} vs. original < '
                   f'{max_diff_lvl}?: %s'%('yes' if diff<max_diff_lvl else 'no'))
        if do_plot:
            plt.plot(wls_sparse/2,xs,label='via laz')
            plt.plot(wls_sparse/2,xs_orig,label='orig')
            plt.legend()
            plt.show()
        if not diff < max_diff_lvl:
            raise SystemExit(f'SC maxdiff test failed with diff: {diff}')

test_cfgstr( 'stdlib::Al_sg225.ncmat;dcutoff=0.6' )
test_cfgstr( 'stdlib::CaSiO3_sg2_Wollastonite.ncmat;dcutoff=3.2' )
test_cfgstr( 'stdlib::Y2SiO5_sg15_YSO.ncmat;dcutoff=2.0' )
test_cfgstr( 'stdlib::Mg2SiO4_sg62_MagnesiumSilicate.ncmat;dcutoff=1.2' )
test_cfgstr( 'stdlib::PbO-alpha_sg129_Litharge.ncmat;dcutoff=0.9' )
test_cfgstr( 'stdlib::MgCO3_sg167_MagnesiumCarbonate.ncmat;dcutoff=1.0' )
test_cfgstr( 'stdlib::Al2O3_sg167_Corundum.ncmat;dcutoff=1.5')
test_cfgstr( 'stdlib::ZnO_sg186_ZincOxide.ncmat;dcutoff=0.5' )
