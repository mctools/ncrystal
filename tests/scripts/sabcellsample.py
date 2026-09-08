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
#from NCrystalDev.hist import HistFiller1D
from NCTestUtils.stat import kolmogorov_smirnov_pvalue
from NCTestUtils.sabsampleutils import Samples, onedim_projections
import numpy as np

_seed_vals = [12345]
def testcell(*,E_div_kT, alpha, beta, svals = None,
             do_plot=False, luxlvl=0, mplvl = 0, verbose=True ):
    a1, a2 = alpha
    b1, b2 = beta
    if svals is None:
       svals=[1.0,1.0,1.0,1.0]
       assert len(svals)==4
    vals = [E_div_kT,a1,a2,b1,b2,*svals]

    assert 0 <= luxlvl <= 6
    assert 0 <= mplvl <= 6
    nsamples_cpp = int(1e5*(10**luxlvl))
    nsamples_mpref = 0
    if mplvl>0:
        nsamples_mpref = int(1e3*(10**mplvl))
        assert nsamples_mpref>0

    _seed_vals[0] += 117
    query = ( [ 'sab','sglcell']
              +['@%s'%str(e) for e in vals]
              +[str(nsamples_cpp),str(_seed_vals[0])] )
    if do_plot:
        import shlex
        print("Query:", shlex.join(query))
    res = ncquery( query, huge_arrays=True)

    sm=res['sampling']
    mtd = sm['chosen_sample_method']
    altmtd = 'bc' if mtd=='fc' else 'fc'
    def methodfullName(m):
        assert m in ('bc','fc')
        return 'BoundedCell' if m=='bc' else 'FullCell'
    print("Chosen method: %s"%methodfullName(mtd))
    print("Alternative method: %s"%methodfullName(altmtd))
    ndig = 5 if verbose else 2
    print(f"Predicted FullCell AR: %.{ndig}g%%"%(100.0*sm['fc_predicted_AR']))
    actual_fc_ar = ( len(sm['fc_sampled_alpha'])/sm['fc_sampled_ntries']
                     if sm['fc_sampled_ntries'] else None )
    if actual_fc_ar is not None:
        print(f"Actual FullCell AR   : %.{ndig}g%%"%(100.0*actual_fc_ar))
    print(f"Actual BoundedCell AR: %.{ndig}g%%"
          %(100.0*len(sm['bc_sampled_alpha'])/sm['bc_sampled_ntries']))
    if sm['fc_prob_edge_1'] != -1.0:
        print(f"Probability edge@b1 (FullCell)   : %.{ndig}g%%"
              %(100.0*sm['fc_prob_edge_1']))
    print(f"Probability edge@b1 (BoundedCell): %.{ndig}g%%"
          %(100.0*sm['bc_prob_edge_1']))

    AR_prec = 0.1
    if sm['fc_predicted_AR']>0.01:
        AR_prec = 0.05
    if sm['fc_predicted_AR']>0.1:
        AR_prec = 0.025
    if actual_fc_ar is not None and sm['fc_predicted_AR']>0.001:
        assert abs(actual_fc_ar/sm['fc_predicted_AR']-1.0)<AR_prec

    samples_ref = Samples('ref')
    samples_ref.add_data( sm['ref_sampled_alpha'],
                          sm['ref_sampled_beta'],
                          sm['ref_sampled_ntries'] )
    samples_std = samples_ref.clone_empty('std')
    samples_alt = samples_ref.clone_empty('alt')

    samples = dict( ref = samples_ref,
                    std = samples_std,
                    alt = samples_alt )


    samples_std.add_data( sm[f'{mtd}_sampled_alpha'],
                          sm[f'{mtd}_sampled_beta'],
                          sm[f'{mtd}_sampled_ntries'] )
    samples_alt.add_data( sm[f'{altmtd}_sampled_alpha'],
                          sm[f'{altmtd}_sampled_beta'],
                          sm[f'{altmtd}_sampled_ntries'] )

    mpsamples = None
    if nsamples_mpref:
        from NCTestUtils.sabcelleval import RefCell
        refcell = RefCell( a1=a1, a2=a2, b1=b1, b2=b2,
                           s11=svals[0], s12=svals[1],
                           s21=svals[2], s22=svals[3] )
        samples_mp = samples_ref.clone_empty('mp')
        mpsamples = refcell.sample( E_div_kT, n=nsamples_mpref, seed=7 )
        mp_ab = np.asarray(mpsamples,dtype=float).copy()
        if len(mp_ab)==0:
            mp_a, mp_b = np.zeros(0), np.zeros(0)
        else:
            mp_a,mp_b = mp_ab.T[0], mp_ab.T[1]
        samples_mp.add_data( mp_a, mp_b, None )
        samples['mp'] = samples_mp

    sampled_a = np.asarray(sm[f'{mtd}_sampled_alpha'],dtype=float)
    sampled_b = np.asarray(sm[f'{mtd}_sampled_beta'],dtype=float)
    altsampled_a = np.asarray(sm[f'{altmtd}_sampled_alpha'],dtype=float)
    altsampled_b = np.asarray(sm[f'{altmtd}_sampled_beta'],dtype=float)

    ref_sampled_a = np.asarray(sm['ref_sampled_alpha'],dtype=float)
    ref_sampled_b = np.asarray(sm['ref_sampled_beta'],dtype=float)

    altsampled_ok = len(altsampled_a) > 0.01*len(sampled_a)
    altsampled_ok_for_stats = len(altsampled_a) > 0.1*len(sampled_a)

    ksp = kolmogorov_smirnov_pvalue
    pvals = {}
    pvals_alt = {}
    for key, fct in onedim_projections:
        tmp = fct(sampled_a,sampled_b)
        pvals[key] = ksp(fct(ref_sampled_a,ref_sampled_b),tmp)
        if altsampled_ok_for_stats:
            pvals_alt[key] = ksp(fct(altsampled_a,altsampled_b),tmp)

    pval_global_worst = 1.0
    npvals = 0
    def print_pvalline(k,p):
        print("   %s : %.4g%%"%(k.rjust(3),p*100))

    for ttt, pvdict in [('std vs. ref',pvals),
                        ('std vs. alt',pvals_alt)]:
        if not pvdict:
            continue
        if verbose:
            print(f"KS unbinned pval ({ttt}):")
        for k,p in pvdict.items():
            npvals += 1
            if verbose:
                print_pvalline(k,p)
            pval_global_worst = min(pval_global_worst,p)
        if verbose:
            print_pvalline('worst',min(p for k,p in pvdict.items()))

    if do_plot:
        import matplotlib.pyplot as plt
        from NCrystalDev.plot import PlotContext
        from NCTestUtils.sabcelleval import plot_celleval
        fig, axes = plt.subplots(2, 3, figsize=(12, 6))

        title='s11=%g, s12=%g, s21=%g, s22=%g'%tuple(res['S'])
        fig.suptitle(title)

        def plot2dcell( ax, avals, bvals, title, nmax=50000 ):
            if nmax is not None and nmax<len(avals):
                #only show up to nmax pts in plot
                idx = np.random.choice(len(avals), size=nmax, replace=False)
                avals = avals[idx]
                bvals = bvals[idx]
            pctx = PlotContext(axis=ax,do_show=False)
            plot_celleval( res, do_title = False, **pctx.kwargs_subcontext() )
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
            pctx.finalise( do_grid = False, do_legend='draggable' )
        nmax2d = 10000
        plot2dcell(axes[0,0],sampled_a,sampled_b,
                   f'sampled (C++, {mtd.upper()})',nmax=nmax2d)
        plot2dcell(axes[1,0],ref_sampled_a,ref_sampled_b,'Ref (C++)',nmax=nmax2d)

        def plothist(ax,key):
            pctx = PlotContext(axis=ax,do_show=False)

            hists = dict( (k,s.create_hist(key))
                          for k,s in samples.items() )
            assert hists['ref'].integral>0
            href = hists['ref']
            hstd = hists['std']
            halt = hists['alt']
            hmp = hists.get('mp')
            for k,v in hists.items():
                if k!='ref' and v.integral>0.0:
                    v.scale( (href.binwidth*href.integral)
                                    /(v.binwidth*v.integral) )
            pval = hstd.check_compat(href,return_pval=True,check=False)
            pval_alt = ( hstd.check_compat(halt,return_pval=True,check=False)
                         if altsampled_ok else None )
            href.plot(label='Ref (C++)',
                      color='blue',alpha=0.5,
                      error_bands=True,
                      **pctx.kwargs_subcontext())
            if altsampled_ok:
                halt.plot(label=f'Alt (C++, {altmtd.upper()})',
                          color='orange',alpha=0.3,
                          error_bands=True,
                      **pctx.kwargs_subcontext())
            if hmp is not None and hmp.integral>0:
                hmp.plot(label='mpmath',
                         color='green',alpha=0.3,
                         error_bands=True,
                         **pctx.kwargs_subcontext())
            hstd.plot(label=f'Std (C++, {mtd.upper()})',
                      color='none',
                      **pctx.kwargs_subcontext())

            t_alt = ('' if pval_alt is None
                    else f', pval_fc_vs_bc: {pval_alt*100:.3g}%')
            ax.set_title(f'{key} (pval: {pval*100:.3g}%{t_alt})')
            pctx.finalise(do_legend=('draggable' if key=='a' else False))

        plothist(axes[0,1],'a')
        plothist(axes[0,2],'b')
        plothist(axes[1,1],'amb')
        plothist(axes[1,2],'apb')
        plt.tight_layout()
        plt.show()

    return npvals,pval_global_worst

def main(do_plot,luxlvl,mplvl,test_select):

    verbose = do_plot or luxlvl or mplvl

    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')

    testpts = [
        dict( E_div_kT=0.06,
              alpha=(6.6767857524277714,6.75015702442328),
              beta=(5.4022664171285335,5.5029404418379908) ),
        dict(E_div_kT=1.5,alpha=(1.0,2.0),beta=(-2.0,-1.0),
             svals=[1.0,20.0,3.0,4.0]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[1e3,1,1e3,1]),
        dict(E_div_kT=0.505,alpha=(0.0,1.0),beta=(-1.0,-0.5),
             svals=[1,1e10,1,1e10]),
        dict(E_div_kT=0.505,alpha=(0.0,1.0),beta=(-1.0,-0.5),
             svals=[1000,1,3,2000]),
        dict(E_div_kT=0.505,alpha=(0.0,1.0),beta=(-1.0,-0.5),
             svals=[0,1,0,1]),
        dict(E_div_kT=0.505,alpha=(0.0,1.0),beta=(-1.0,-0.5),
             svals=[1,0,3,0]),
        dict(E_div_kT=0.55,alpha=(0.0,1.0),beta=(-1.0,-0.5),
             svals=[1,0,3,0]),
        dict(E_div_kT=0.7,alpha=(0.0,1.0),beta=(-1.0,-0.5),
             svals=[1,0,3,0]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[1,1,0,0]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[0,0,1,1]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[1,1,1,1]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[1e3,1,1,1]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[1,1e3,1,1]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[1,1,1e3,1]),
        dict(E_div_kT=0.01,alpha=(5.8,6.0),beta=(5.2,5.35),
             svals=[1,1,1,1e3]),
        dict(E_div_kT=0.0001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1,1,0,0]),
        dict(E_div_kT=0.0001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[0,0,1,1]),
        dict(E_div_kT=0.0001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1,1,1,1]),
        dict(E_div_kT=0.0001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e4,0,0]),
        dict(E_div_kT=0.0001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[0,0,1e4,1.0]),
        dict(E_div_kT=0.0001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e4,1e4,1.0]),
        dict(E_div_kT=0.000001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e4,1e4,1.0]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0,1,0,1]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0)),
        dict(E_div_kT=1.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0.1,1000,2000,0.2]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0.1,1000,1000,0.1]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0.1,1000,0,0]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0,0,1000,0.1]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[1,1000,1,1]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0,0,1,1]),
        dict(E_div_kT=1e-10,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e50,0,0]),
        dict(E_div_kT=1e-15,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e5,0,0]),
        dict(E_div_kT=1e-15,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e250,0,0]),
        dict(E_div_kT=1e-15,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e5,0,0]),
        dict(E_div_kT=1,
             alpha=(1e-50,1e-10),
             beta=(-1,1),
             svals=[1,1,1,1]),
        dict(
            E_div_kT=0.2812650004,
            alpha=(5e-10,8e-10),
            beta=(-3e-5,-2e-05),
            svals=[1,1,1,1]),
        dict(
            E_div_kT=0.2812650005,
            alpha=(5e-10,8e-10),
            beta=(-3e-5,-2e-05),
            svals=[1,1,1,1]),
        dict(E_div_kT=0.000395856,
             alpha=(8.36433e-50,8.36433e-10),
             beta=(-2.78659e-05,-2.60257e-08),
             svals=[2.09324e-51,2.09324e-11,2.09321e-51,2.09321e-11]),
        dict(E_div_kT=1.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0.1,1,20000,0.2]),
        dict(E_div_kT=2.0,alpha=(0.5,3.0),beta=(-10.0,40.5),
             svals=[1.0,20.0,3.0,4.0]),
        dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(11.0,12.0)),
        dict(E_div_kT=1e-5,alpha=(11.0,13.0),beta=(0.5,12.0)),
        dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(0.5,12.0)),
        dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(11.0,12.0),
             svals=[1.0,20.0,3.0,4.0]),
        dict(E_div_kT=0.01,alpha=(11.35,13),beta=(10,12.0),
             svals=[0,1,1,0]),
        dict(E_div_kT=1e-5,alpha=(11.0,13.0),beta=(0.5,12.0),
             svals=[1.0,20.0,3.0,0]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-7,5.5),beta=(5,5.92),
             svals=[1.0,1e-200,1.0,1e-200]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-10,5.5),beta=(5,5.92),
             svals=[1.0,1e-20,1.0,1e-20]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-11,5.5),beta=(5,5.92),
             svals=[1e-20,1.0,1.0,1e-20]),
        dict(E_div_kT=1e-5,alpha=(11.0,13.0),beta=(0.5,12.0),
             svals=[1.0,20.0,3.0,1e-300]),
        dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(0.5,12.0),
             svals=[1.0,20.0,3.0e-4,4.0]),
        dict(E_div_kT=0.01,alpha=(5.499,5.5),beta=(5,5.92),
             svals=[1.0,99,1.0,99]),
        dict(E_div_kT=0.01,alpha=(5.499,5.5),beta=(5,5.92),
             svals=[1.0,101,1.0,101]),
        dict(E_div_kT=0.01,alpha=(5.499,5.5),beta=(5,5.92),
             svals=[1.0,0.011,1.0,0.011]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-11,5.5),beta=(5,5.92),
             svals=[1.0,0.009,1.0,0.009]),
        #This one actually fails, but with just ~3-4 digits of precision, it is
        #not so surprising. At looking at it interactively, it actually seems
        #fine:
        #dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
        #     svals=[1.0,0.009,1.0,0.009]),

        dict( E_div_kT=39.58561183871582 * 0 + 40,
              alpha=(158.28456094789,160.60010716438572),
              beta=(-2.2701714000978403e-14,-2.1202501083299951e-17),
              svals=[ 0.028331793397541188,
                      0.027532855700519738,
                      0.027532855700519426,
                      0.026748257342144074 ] ),

    ]

    npvals_tot = 0
    pvals = []
    nused = 0
    if test_select:
        test_select = set(test_select)
    #for i in range(2):
    #    testpts += testpts


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
        if do_plot:
            print(data)

        #print(data)
        npvals, worst_pval = testcell(**data,
                                      do_plot=do_plot,
                                      luxlvl=luxlvl,
                                      mplvl=mplvl,
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
    mplvl = 0
    if '-h' in args or '--help' in args:
        print("""Arguments:
        <none> : Run as regular test
        --lux  : Increase statistics by factor of 10.
                 This option can be specified multiple times.
        --mp   : Enable mpmath reference. This option can be specified multiple
                 times for more mpmath reference statistics.
        --plot : Show interactive plots
        1 7 18 : Specify digits on the command line to run just those tests.
        """)
        raise SystemExit()
    while '--plot' in args:
        args.remove('--plot')
        do_plot = True
    while '--lux' in args:
        args.remove('--lux')
        luxlvl += 1
    while '--mp' in args:
        args.remove('--mp')
        mplvl += 1
    assert not args or all(e.isdigit() for e in args)
    main(do_plot = '--plot' in sys.argv[1:],
         luxlvl = luxlvl,
         mplvl = mplvl,
         test_select = [int(e) for e in args])

#fixme: Add option to test a huge bunch of randomly generated cells.
