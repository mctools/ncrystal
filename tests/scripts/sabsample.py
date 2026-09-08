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
from NCTestUtils.env import ncsetenv
from NCrystalDev.misc import evaluate_query as ncquery
from NCTestUtils.stat import kolmogorov_smirnov_pvalue
from NCTestUtils.sabsampleutils import Samples, onedim_projections
from NCTestUtils import sabsampleutils
from dataclasses import dataclass
from typing import Optional
import numpy as np

_seed_cache = [0]
def gen_seedstr():
    _seed_cache[0] += 117
    return str(_seed_cache[0])
def reset_seed():
    _seed_cache[0] = 12345

@dataclass(frozen=True)
class SampleResults:
    method: str
    alpha: np.ndarray
    beta: np.ndarray
    egrid_div_kT : Optional[np.ndarray] = None # noqa FA100

sample_methods =  ['ref','std','legacy','legacy_oversample','vdoslux5']
def dosample_query( method, ekin, nsamples, cfgstr, atomlbl, egrid = None ):
    assert egrid is None, "not implemented"#fixme: implement and use!
    assert method in sample_methods
    seedstr = gen_seedstr()
    if method == 'vdoslux5':
        method = 'std'
        cfgstr += ';vdoslux=5;knllux=6'
        #fixme: egrid with just a few points (one at 2*ekin, the rest very low)
    if method == 'std':
        query = [ 'sab','proc', cfgstr, atomlbl,
                  str(nsamples), seedstr, str(ekin) ]
        reskey = 'sample'
    else:
        assert method in ['ref','legacy','legacy_oversample']
        query = ["sab","refsample", cfgstr, atomlbl,
                 str(nsamples), seedstr, str(ekin), method ]
        reskey = 'refsample' if method=='ref' else 'legacysample'

    r = ncquery( query, huge_arrays=True)
    egrid_div_kT = None
    if 'sabproc' in r and 'egrid' in r['sabproc']:
        egrid_div_kT = r['sabproc']['egrid']

    if reskey:
        r = r[reskey]
    return SampleResults( alpha = r['alpha'],
                          beta = r['beta'],
                          method = method,
                          egrid_div_kT = egrid_div_kT )

def sample(template = None, get_egrid = False, **kwargs):
    sr = dosample_query( **kwargs )
    if template:
        s = template.clone_empty(sr.method)
    else:
        s = Samples(name=sr.method)
    s.add_data( avals = sr.alpha, bvals = sr.beta, keepdata = True )
    return ( s, sr.egrid_div_kT ) if get_egrid else s

def pval_from_samples( samples1, samples2 ):
    ksp = kolmogorov_smirnov_pvalue
    pvals = {}
    pvals_worst = 2.0
    for key, fct in onedim_projections:
        tmp = fct(samples1.alpha,samples1.beta)
        pv = ksp(fct(samples2.alpha,samples2.beta),tmp)
        pvals_worst = min( pv, pvals_worst )
        pvals[key] = pv
    return pvals, pvals_worst

def test( *,cfgstr,ekin, atomlbl='',
          do_plot=False, luxlvl=0, verbose=True ):

    nsamples = int( 1e5 * ( 10**luxlvl ) )
    #def dosample( method, ekin, nsamples, cfgstr, atomlbl, egrid = None ):
    common = { 'cfgstr' : cfgstr, 'atomlbl' : atomlbl,
               'ekin' : ekin, 'nsamples' : nsamples, 'egrid' : None }
    s_ref =  sample( method = 'ref', **common )
    s_std, s_std_egrid_div_kT =  sample( method = 'std', template=s_ref,
                                         get_egrid = True,  **common )

    samples = { 'ref' : s_ref, 'std' : s_std }

    pvals, pvals_worst = pval_from_samples( s_std, s_ref )

    if verbose:
        print("KS unbinned pval (std vs. ref):")
        def print_pvalline(k,p):
            print("   %s : %.4g%%"%(k.rjust(3),p*100))
        for k,p in pvals.items():
            print_pvalline(k,p)
        print_pvalline('worst',pvals_worst)

    if luxlvl!=0:
        for m in sample_methods:
            if m not in samples:
                samples[m] = sample( method = m, template=s_ref, **common )

    if do_plot:
        import matplotlib.pyplot as plt
        import NCrystalDev.cfgstr as nccfgstr
        import NCrystalDev.core as nccore
        from NCrystalDev._numpy import _np_linspace
        from NCrystalDev.plot import PlotContext
        fig, axes = plt.subplots(2, 3, figsize=(12, 6))
        fig.suptitle(cfgstr + (f' ({atomlbl})' if atomlbl else '')+f' {float(ekin):g}eV')

        from NCrystalDev.constants import constant_boltzmann
        info = nccore.createInfo(cfgstr)

        if atomlbl:
            di = [e for e in info.dyninfos if e.lbl==atomlbl]
            assert len(di)==1
            di = di[0]
        else:
            assert len(info.dyninfos)==1
            di = info.dyninfos[0]
        kT = di.temperature*constant_boltzmann
        s_std_emax_div_kT = s_std_egrid_div_kT[-1]
        #s_std_emax = s_std_emax_div_kT * kT
        E_div_kT = float(ekin) / kT
        if hasattr(di,'analyseVDOS'):
            vdoslux = nccfgstr.decodecfg_vdoslux(cfgstr)
            sab = di.loadKernel( vdoslux=vdoslux )
        else:
            sab = di.loadKernel()

        def plot2dcell( ax, avals, bvals, title,
                        amin, amax, bmin, bmax,
                        agrid, bgrid, emax_div_kT = None ):
            nmax = 20000
            if nmax<len(avals):
                #only show up to nmax pts in plot
                idx = np.random.choice(len(avals), size=nmax, replace=False)
                avals = avals[idx]
                bvals = bvals[idx]
            pctx = PlotContext(axis=ax,do_show=False)
            eps = 0.05
            b0 = bmin - eps*(bmax-bmin)
            b1 = bmax + eps*(bmax-bmin)
            a0 = max( 0.0, amin - eps*(amax-amin) )
            a1 = amax + eps*(amax-amin)
            pctx.axis.set_xlim( b0, b1 )
            pctx.axis.set_ylim( a0, a1 )
            ax.set_xlabel('beta')
            ax.set_ylabel('alpha')

            linecommon=dict(color='gray', linestyle=':', linewidth=1, alpha=0.7)
            ax.hlines( [ v for v in agrid
                         #if amin <= v <= amax
                        ],
                       xmin = bgrid[0], xmax=bgrid[-1],
                       **linecommon )
            ax.vlines( [ v for v in bgrid
                         #if bmin <= v <= bmax
                        ],
                       ymin = agrid[0], ymax=agrid[-1],
                       **linecommon )

            #Phase-space edge:
            def print_kb( e, color ):
                assert bmax>=-e#revisit this if it fails
                brangeplot = [max(b0,-e),b1]
                assert brangeplot[1] > brangeplot[0]
                assert e > 0 and not np.isinf(e)
                b = _np_linspace(*brangeplot,5000)
                sbe = np.sqrt(b+e)
                ap = ( sbe + np.sqrt(e) )**2
                am = ( sbe - np.sqrt(e) )**2
                pctx.axis.plot(b,ap,color=color,alpha=0.5)
                pctx.axis.plot(b,am,color=color,alpha=0.5)
            print_kb(E_div_kT,color='red')
            if emax_div_kT is not None:
                print_kb(emax_div_kT,color='green')

            n = len(avals)
            if n>0:
                base_size = 40.0
                s = base_size * max(0.12, (1.0 / np.sqrt(n / 200.0)))
                alpha = min(0.9, max(0.02, 2000.0 / (n + 1)))
                marker = '.' if n > 200000 else 'o'
                pctx.axis.scatter(bvals,avals,
                                  s=s, c='C0', alpha=alpha, marker=marker,
                                  linewidths=0, rasterized=(n >= 200000) )
            pctx.axis.set_title(title)
            pctx.finalise( do_grid = False )

        #same range in 2d plots:
        amin = min( s_std.alpha.min(), s_ref.alpha.min() )
        amax = max( s_std.alpha.max(), s_ref.alpha.max() )
        bmin = min( s_std.beta.min(), s_ref.beta.min() )
        bmax = max( s_std.beta.max(), s_ref.beta.max() )
        common2d = dict( amin = amin, amax = amax,
                         bmin = bmin, bmax = bmax,
                         agrid = sab['alpha'],
                         bgrid = sab['beta'] )

        plot2dcell( axes[0,0],s_std.alpha,s_std.beta,'sampled (std)',
                    emax_div_kT = s_std_emax_div_kT,
                    **common2d,  )
        #NB: The ref plot will also get egrid from s_std:
        plot2dcell( axes[1,0],s_ref.alpha,s_ref.beta,'sampled (ref)',
                    emax_div_kT = s_std_emax_div_kT,
                    **common2d )

        def plothist(ax,key,title):
            pctx = PlotContext(axis=ax,do_show=False)

            hists = dict( (k,s.create_hist(key))
                          for k,s in samples.items() )
            assert hists['ref'].integral>0
            href = hists['ref']
            hstd = hists['std']
            for k,v in hists.items():
                if k!='ref' and v.integral>0.0:
                    v.scale( (href.binwidth*href.integral)
                                    /(v.binwidth*v.integral) )
            href.plot(label='Ref',
                      color='blue',alpha=0.5,
                      error_bands=True,
                      **pctx.kwargs_subcontext())
            colors={ 'ref' : 'blue',
                     'legacy':'orange',
                     'legacy_oversample':'purple',
                     'vdoslux5':'red',
                     'std':'none' }
            assert set(colors.keys()) == set(sample_methods)
            hstd.plot(label='std',
                      color=colors['std'],
                      **pctx.kwargs_subcontext())
            for k,v in hists.items():
                if k not in ('ref','std'):
                    v.plot(label=k,
                           error_bands=True,
                           color=colors[k],
                           alpha=0.5,
                           **pctx.kwargs_subcontext())

            ax.set_title(f'{title} (pval: {pvals[key]*100:.3g}%)')
            assert sum(int(k=='a') for k,v in onedim_projections)==1
            pctx.finalise(do_legend=('draggable' if key=='a' else False),
                          do_grid = True)

        plothist(axes[0,1],'a','alpha')
        plothist(axes[0,2],'b','beta')
        plothist(axes[1,1],'amb','alpha-beta')
        plothist(axes[1,2],'apb','alpha+beta')
        plt.tight_layout()
        plt.show()

    return len(pvals), pvals_worst

def main(do_plot,luxlvl,test_select):
    verbose = do_plot or luxlvl
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')

    testpts = [
        #Fixme: varying vdoslux in CaH2@20K shows that we perhaps do not have an
        #ideal alpha-beta grid generated, since vdoslux 4 or 5 are needed to
        #remove heavy artifacts:

        #dict( cfgstr='solid:H/1gcm3',atomlbl='H',ekin='10.0');,

        dict( cfgstr='stdlib::Polyethylene_CH2.ncmat;knllux=6',atomlbl='C',#fixme: also 'C'
              ekin='15'#fixme: also something extreme, like 10000
             ),

        dict( cfgstr='stdlib::Al_sg225.ncmat;knllux=6',ekin='1000.5'),#fixme: also something normal

        dict( cfgstr='stdlib::CaH2_sg62_CalciumHydride.ncmat;temp=20;vdoslux=2;knllux=1',
              atomlbl='H',
              ekin='0.0016694736654149164'#7Aa
             ),
        dict( cfgstr='stdlib::CaH2_sg62_CalciumHydride.ncmat;temp=20;vdoslux=2',
              atomlbl='H',
              ekin='0.0016694736654149164'#7Aa
             ),
        dict( cfgstr='stdlib::CaH2_sg62_CalciumHydride.ncmat;temp=20;vdoslux=2',
              atomlbl='Ca',
              ekin='0.0016694736654149164'#7A
             ),
        dict( cfgstr='stdlib::Al_sg225.ncmat;vdoslux=0;temp=500', ekin='0.025' ),

        dict( cfgstr='stdlib::Al_sg225.ncmat', ekin='0.025' ),

        dict( cfgstr='stdlib::Al_sg225.ncmat', ekin='0.0025' ),

        dict( cfgstr='stdlib::Al_sg225.ncmat', ekin='25e-10' ),
        #ref is too slow dict( cfgstr='stdlib::Al_sg225.ncmat', ekin='25e-10' ),

        dict( cfgstr='stdlib::Al_sg225.ncmat', ekin='15' ),

        dict( cfgstr='solid::H/1gcm3;vdoslux=5', ekin='15' ),

        dict( cfgstr='solid::H/1gcm3;vdoslux=0', ekin='2' ),

        dict( cfgstr='stdlib::Polyethylene_CH2.ncmat;knllux=6',atomlbl='H',#fixme: also 'C'
              ekin='15'#fixme: also something extreme, like 10000
             ),

        dict( cfgstr='stdlib::Be_sg194.ncmat', ekin='10' ),

        dict( cfgstr='stdlib::Ca_sg229_Calcium-gamma.ncmat;vdoslux=0;temp=800', ekin='0.04' ),

        dict( cfgstr='stdlib::Polyethylene_CH2.ncmat;knllux=6',
              atomlbl='H',#fixme or C?
              ekin='15' ),


        #dict( cfgstr='stdlib::Polyethylene_CH2.ncmat;temp=57.8755;knllux=5',
        #      atomlbl='H',#fixme or C?
        #      ekin='0.09' ),
        #
        #dict( cfgstr='stdlib::Polyethylene_CH2.ncmat;temp=57.8760;knllux=5',
        #      atomlbl='H',#fixme or C?
        #      ekin='0.09' ),
        #
        #dict( cfgstr='stdlib::Polyethylene_CH2.ncmat;temp=57.8755;knllux=5',
        #      atomlbl='C',#fixme or C?
        #      ekin='0.09' ),
        #
        #dict( cfgstr='stdlib::Polyethylene_CH2.ncmat;temp=57.8760;knllux=5',
        #      atomlbl='C',#fixme or C?
        #      ekin='0.09' ),


    ]
    npvals_tot = 0
    pvals = []
    nused = 0
    if test_select:
        test_select = set(test_select)

    for i,data in enumerate(testpts):
        i += 1
        if test_select and i not in test_select:
            print("=============> SKIPPING test %i"%i)
            continue
        nused += 1
        print()
        print("=============>")
        print("=============> Launching test %i"%i)
        print("=============>")
        print()
        npvals, worst_pval = test(**data,
                                  do_plot=do_plot,
                                  luxlvl=luxlvl,
                                  verbose=verbose)
        npvals_tot += npvals
        pvals.append( (worst_pval,i) )

    if not nused:
        raise SystemExit('ERROR: No tests were run!')

    print()
    print("=============>")
    print("=============> Final summary")
    print("=============>")
    print()
    #FIXME: Update next comment and the actual pval_threshold later:
    #I would have expected something like 0.05/npvals_tot to work for
    #pval_threshold, but I actually have to increase it to something like
    #10.0/npvals_tot to see any false positives. So we ad-hoc ensure that the
    #threshold is not below 0.1%:
    pval_threshold = max(1e-3,0.05/npvals_tot)
    test_error = False

    for pval, i in pvals:
        ok = pval>pval_threshold
        okmsg = ( 'OK' if ok else 'FAILED (pval: %.4g%%, not above %.4g%%)'
                  %(pval*100.0,pval_threshold*100.0) )
        print("P-value for test #%i : %s"%(i,okmsg))
        if not ok:
            test_error = True

    if test_error:
        raise SystemExit('ERROR: Some comparisons failed')

    if test_select and len(test_select)!=nused:
        raise SystemExit('ERROR: Test selection selected one or more'
                         ' non-existent tests!')


if __name__ == '__main__':
    import sys
    do_plot = False
    args = sys.argv[1:]
    luxlvl = 0
    if '-h' in args or '--help' in args:
        print("""Arguments:
        <none> : Run as regular test
        --lux  : Increase statistics by factor of 10 and show more curves.
                 This option can be specified multiple times.
        --plot : Show interactive plots
        1 7 18 : Specify digits on the command line to run just those tests.
        """)
        raise SystemExit()
    while '--plot' in args:
        args.remove('--plot')
        do_plot = True
    while '--lux' in args:
        args.remove('--lux')
        sabsampleutils.hist1d_nbins = int(sabsampleutils.hist1d_nbins*1.5)
        luxlvl += 1
    assert not args or all(e.isdigit() for e in args)
    main(do_plot = '--plot' in sys.argv[1:],
         luxlvl = luxlvl,
         test_select = [int(e) for e in args])

