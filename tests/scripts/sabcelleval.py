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

# NEEDS: numpy mpmath
#
import mpmath # noqa F401
import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.stabilise_ncpprint # noqa F401
from NCTestUtils.sabcelleval import RefCell, draw_alpha_beta_grid
from NCrystalDev._common import ncpprint
from NCrystalDev._numpy import _np_linspace
from NCTestUtils.env import ncsetenv
from NCrystalDev.misc import evaluate_query as ncquery
import numpy as np

def plot_celleval( data, **kw_plot ):
    from NCrystalDev.plot import PlotContext
    pctx = PlotContext(**kw_plot).check_unused()
    draw_alpha_beta_grid( alphagrid = data['alpha'],
                          betagrid = data['beta'],
                          **pctx.kwargs_subcontext() )
    b = data['beta']
    db = b[1]-b[0]
    a = data['alpha']
    da = a[1]-a[0]
    blim = ( b[0]-0.1*db, b[1]+0.1*db )
    alim = ( max(0.0,a[0]-0.1*da), a[1]+0.1*da )
    pctx.axis.set_xlim( *blim )
    pctx.axis.set_ylim( *alim )

    data_ci = data['cellintegral']
    elist = [ ( data['surveyor']['E_div_kT_touch'], 'touch', 'green', 0.0 ),
              ( data['surveyor']['E_div_kT_cover'], 'cover', 'blue',
                data_ci['full_integral'] ),
              ( data_ci['phasespace_E_div_kT'], 'chosen', 'red',
                list( v for k,v in data_ci['phasespace_integral']
                      if k=='Romberg33' )[0] ) ]
    for e, lbl, col, integral in elist:
        brangeplot = [max(blim[0],-e),blim[1]]
        assert brangeplot[1] > brangeplot[0]
        lble = f'{e:g}kT' if not np.isinf(e) else 'INF'
        lbl = f'{lbl} ({lble}, integral={integral:g})'
        if e > 0 and not np.isinf(e):
            b = _np_linspace(*brangeplot,5000)
            sbe = np.sqrt(b+e)
            ap = ( sbe + np.sqrt(e) )**2
            am = ( sbe - np.sqrt(e) )**2
            pctx.axis.plot(b,ap,color=col,label=lbl)
            pctx.axis.plot(b,am,color=col)
        elif np.isinf(e):
            pctx.axis.plot(brangeplot,[0,0],color=col,label=lbl)
        else:
            pctx.axis.plot(brangeplot,brangeplot,color=col,label=lbl)

    e = data_ci['phasespace_E_div_kT']
    b1, b2 = data['beta']
    a1, a2 = data['alpha']
    from matplotlib import patches
    for ( r_a1, r_a2,clip_betaminus,
          clip_betaplus ) in data_ci['integration_regions']['regions']:
        print(f"AlphaRange [{r_a1},{r_a2}]: clip_betaminus"
              f"={clip_betaminus}, clip_betaplus={clip_betaplus}")

        assert a1 <= r_a1 <= a2
        color=None
        if not ( clip_betaminus or clip_betaplus ):
            #Just a square!
            r = patches.Rectangle((b1, r_a1), b2-b1, r_a2-r_a1,
                                  facecolor=color, edgecolor='k',#'lightblue'
                                  hatch='///', linewidth=1.0)
            pctx.axis.add_patch(r)
            continue
        for aval in _np_linspace(r_a1,r_a2,50):
            #phasespace curve: 4ae=(b-a)^2 <=> |b-a|=sqrt(4ae)
            db = np.sqrt(4*aval*e)
            bm = aval-db if clip_betaminus else b1
            bp = aval+db if clip_betaplus else b2
            _=pctx.axis.plot([bm,bp],[aval,aval],color=color,alpha=0.3)
            if color is None:
                color=_[0].get_color()

    title='s11=%g, s12=%g, s21=%g, s22=%g'%tuple(data['S'])
    pctx.axis.set_title(title)
    return pctx.finalise( do_grid = False, do_legend='draggable' )

def evalcell(*,E_div_kT, alpha, beta, svals = None, do_plot=False ):
    a1, a2 = alpha
    b1, b2 = beta
    if svals is None:
        svals=[1.0,1.0,1.0,1.0]
    assert len(svals)==4
    vals = [E_div_kT,a1,a2,b1,b2,*svals]
    res = ncquery( [ 'sab','sglcell']+['@%s'%str(e) for e in vals])
    refcell = RefCell( a1=a1, a2=a2, b1=b1, b2=b2,
                       s11=svals[0], s12=svals[1],
                       s21=svals[2], s22=svals[3] )
    print("Results from C++:")
    ncpprint(res)
    nbf = 1000
    print(f"Brute force results (n={nbf}):")
    bfres = refcell.phasespace_integral_brute_force( E_div_kT, n=nbf )
    ncpprint(bfres)
    mprefval = refcell.phasespace_integral(E_div_kT)
    mprefval_fullint = refcell.full_integral()
    assert mprefval<=mprefval_fullint
    print(f'mpmath reference phasespace integral: {float(mprefval):.15g}')
    print(f'mpmath reference full integral: {float(mprefval_fullint):.15g}')

    bfprec_fullint = float(abs(bfres['full_integral']/mprefval_fullint-1))
    bfprec_psint = float(abs(bfres['phasespace_integral']/mprefval-1))
    print(f"bruteforce/mpref precision (full integral): {bfprec_fullint:g}")
    print(f"bruteforce/mpref precision (phasespace integral): {bfprec_psint:g}")
    assert bfres['phasespace_integral'] <= bfres['full_integral']

    tgt_bfrec_fullint = 0.005
    tgt_bfrec_psint = 0.01
    if bfres['phasespace_integral']<0.01*bfres['full_integral']:
        tgt_bfrec_psint = 0.06
    if a2<a1*(1+1e-10):
        tgt_bfrec_fullint = 0.2
        tgt_bfrec_psint = 0.2

    assert bfprec_fullint < tgt_bfrec_fullint
    assert bfprec_psint < tgt_bfrec_psint

    vals = sorted( ( float(abs(v/mprefval-1)), v, name )
                   for name, v
                   in res['cellintegral']['phasespace_integral'] )
    for prec, v, name in vals:
        print(f" {name.rjust(10)} : {v:.12g}  [precision lvl {prec:g}]")
    resfullint = res['cellintegral']['full_integral']
    prec = float(abs(resfullint/mprefval_fullint-1))
    print(f" full integral : {resfullint:.12g} [precision lvl {prec:g}]")

    r33 = [e for e in vals if e[2]=='Romberg33'][0]
    if do_plot:
        plot_celleval( res )

    r33prec = 5e-6
    if a2<a1*(1+1e-10):
        r33prec = 1e-3
    assert r33[0] < r33prec, "Romberg33 not suitable as reference"

def main(do_plot):
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')

    testpts = [
        dict( E_div_kT=0.06,
              alpha=(6.6767857524277714,6.75015702442328),
              beta=(5.4022664171285335,5.5029404418379908) ),
        dict(E_div_kT=1.5,alpha=(1.0,2.0),beta=(-2.0,-1.0)),
        dict(E_div_kT=1.5,alpha=(1.0,2.0),beta=(-2.0,-1.0),
             svals=[1.0,20.0,3.0,4.0]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0)),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,-0.5)),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-0.1,6.0)),
        dict(E_div_kT=2.0,alpha=(0.0,12.0),beta=(-3.0,10.5)),
        dict(E_div_kT=2.0,alpha=(0.0,6.0),beta=(-1.0,4.5)),
        dict(E_div_kT=2.0,alpha=(0.5,3.0),beta=(-10.0,40.5)),
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
        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,1e-200,1.0,1e-200]),
        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,1e-20,1.0,1e-20]),
        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
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
        #Note: some rather bad results, attributed to numerical issues when a2
        #~- a1.
        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,99,1.0,99]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,101,1.0,101]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,0.011,1.0,0.011]),
        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,0.009,1.0,0.009]),
    ]

    for i,data in enumerate(testpts):
        i += 1
        #if i not in (13,17):
        #    continue
        print()
        print("=============>")
        print("=============> Launching test %i"%i)
        print("=============>")
        print()
        evalcell(**data,do_plot=do_plot)


if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:])

#fixme: Add option to test a huge bunch of randomly generated cells.
