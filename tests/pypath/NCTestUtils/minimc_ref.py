
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

# Utilities for detecting changes in given MiniMC output, and making it easy to
# inspect differences and update reference histograms in case of changes.
import NCrystalDev._mmc as ncmmc
from NCrystalDev._hist import Hist1D
from NCrystalDev.core import load as ncload

#Specialisation for std case with which we might want to test many (or even all)
#materials.:
def minimc_unittest_stdsphere( *,
                               cfgstr,
                               neutron_energy,#string like "1eV" or "1Aa"
                               illuminate_uniformly = False,
                               **kwargs #passed to minimc_unittest
                               ):
    mat = ncload( cfgstr )
    ekin = ncmmc._parse_energy( neutron_energy )
    ( sphere_diam_str,
      sphere_diam ) = ncmmc._approx_mfp_as_length_string( mat, ekin = ekin )
    sphere_radius = sphere_diam/2
    n = 1e4 if mat.scatter.isOriented() else 1e5
    srcz = (-sphere_radius)*(1-1e-13)

    if illuminate_uniformly:
        srccfg = f'circular;r={sphere_radius}'
    else:
        srccfg = 'constant'
    srccfg += f';ekin={ekin};z={srcz};n={n}'
    minimc_unittest( cfgstr = cfgstr,
                     srccfg = srccfg,
                     geomcfg = f'sphere;r={sphere_radius}',
                     **kwargs )

def _parse_sysargv():
    import sys
    args = sys.argv[1:]
    if '--help' in args or '-h' in args:
        print("Possible options: --plot or --update")
    return dict( do_plot = '--plot' in args,
                 do_updateref = '--update' in args )

def main_minimc_unittest_stdsphere( *a, **kw ):
    minimc_unittest_stdsphere( *a, **kw, **_parse_sysargv() )

def _detect_caller_filebasename( drop_ext = True ):
    #Find first file in callstack outside the present one.
    import inspect
    import os
    res = None
    for e in inspect.stack():
        if not os.path.samefile(__file__,e.filename):
            res = e.filename
            break
    assert res
    res = os.path.basename( res )
    if drop_ext:
        res = os.path.splitext( res )[0]
    return res

def minimc_unittest( *,
                     key = '<auto>', #<auto> expands to name taken from
                                     #filename of caller.
                     cfgstr,
                     srccfg,
                     geomcfg,
                     do_plot = False,
                     do_updateref = False,
                    ):
    """Run MiniMC with chosen configuration and compare result with reference
    histogram. Will keep reference histogram in <testsdata>mmcrefs/key.json (if
    key contains a slash, a subdir will be used). If do_plot is True, show
    plot vs. ref., and if do_updateref is True, update the reference.

    Runs with nthreads=2.
    """
    from .dirs import get_named_test_data_dir

    if '/' not in key:
        key = f'mmcref/{key}'
    if '<auto>' in key:
        autokey = _detect_caller_filebasename()
        if autokey.startswith('sb_nctest'):
            autokey = autokey.split('_',2)[2]
            assert autokey.startswith('test')
            autokey=autokey[len('test'):]
        key = key.replace('<auto>',autokey)
    assert key.count('/') <= 1, "no more than one slash allowed in key"
    if '/' in key:
        subdir, key = key.split('/')
        assert subdir and key
    else:
        subdir, key = None, key

    reffile = get_named_test_data_dir(
        subdir,
        for_updates = do_updateref
    ).joinpath( key+'.json' )

    res = ncmmc.runsim_diffraction_pattern( cfgstr,
                                            srccfg = srccfg,
                                            geomcfg = geomcfg,
                                            nthreads = 2 )
    h = res.histogram_main.clone(rebin_factor=20)#smaller data files
    if do_plot:
        res.plot_breakdown(rebin_factor=10,logy=True)
        #h.plot(logy=True)

    if do_updateref:
        reffile.write_text(h.to_json())
        print(f"Updated {reffile}")
        return

    h_ref = Hist1D(reffile.read_text())
    pval = h.check_compat( h_ref, return_pval = True )
    print(f"Pvalue for comp. with ref (higher is more compatible): {pval:g}")
    if do_plot:
        h_ref.plot(do_show=False,error_bands=1.0,
                   alpha=0.3,color='blue',label='ref')
        plt=h.plot(do_show=False,color='none',logy=True,label='new')
        plt.legend()
        plt.grid()
        plt.show()
    if pval < 0.001:
        raise SystemExit("ERROR: Possible compatibility issues detected!")
