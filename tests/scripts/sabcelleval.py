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
import NCTestUtils.enable_fpe
import NCTestUtils.stabilise_ncpprint # noqa F401
from NCTestUtils.sabcelleval import RefCell
from NCrystalDev._common import ncpprint
from NCTestUtils.env import ncsetenv
from NCrystalDev.misc import evaluate_query as ncquery

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
    def fmtprec(v):
        if do_plot:
            return '%g'%v
        if v < 1e-14:
            return '<1e-14'
        return '%.1g'%v

    for prec, v, name in vals:
        print(f" {name.rjust(10)} : {v:.11g}  [precision lvl {fmtprec(prec)}]")
    resfullint = res['cellintegral']['full_integral']
    prec = float(abs(resfullint/mprefval_fullint-1))
    print(f" full integral : {resfullint:.11g} [precision lvl {fmtprec(prec)}]")

    f65 = [e for e in vals if e[2]=='Flex65'][0]
    if do_plot:
        from NCTestUtils.sabcelleval import plot_celleval
        plot_celleval( res )

    f65prec = 5e-6
    if a2<a1*(1+1e-10):
        f65prec = 1e-3
    assert f65[0] < f65prec, "Romberg65 not suitable as reference"

def main(do_plot,test_select):
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
        #~- a1 (which is to be expected, since there are not many mantissa bits
        #left to describe values inside the span from e.g. 5.5-1e-12 to 5.5.

        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,99,1.0,99]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,101,1.0,101]),
        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,0.011,1.0,0.011]),
        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
             svals=[1.0,0.009,1.0,0.009]),
        dict( E_div_kT=0.2812650004,
              alpha=(5e-10,8e-10),
              beta=(-3e-5,-2e-05),
              svals=[1,1,1,1] ),
        dict( E_div_kT=0.2812650005,
              alpha=(5e-10,8e-10),
              beta=(-3e-5,-2e-05),
              svals=[1,1,1,1] ),
        dict(E_div_kT=0.000395856,
             alpha=(8.36433e-50,8.36433e-10),
             beta=(-2.78659e-05,-2.60257e-08),
             svals=[2.09324e-51,2.09324e-11,2.09321e-51,2.09321e-11]),
    ]

    if test_select:
        test_select = set(test_select)

    for i,data in enumerate(testpts):
        i += 1
        if test_select and i not in test_select:
            print("=============> SKIPPING test %i"%i)
            continue
        print()
        print("=============>")
        print("=============> Launching test %i"%i)
        print("=============>")
        print()
        evalcell(**data,do_plot=do_plot)


if __name__ == '__main__':
    import sys
    do_plot = False
    args = sys.argv[1:]
    while '--plot' in args:
        args.remove('--plot')
        do_plot = True
    assert not args or all(e.isdigit() for e in args)

    main( do_plot = do_plot,
          test_select = [int(e) for e in args] )

#fixme: Add option to test a huge bunch of randomly generated cells.
