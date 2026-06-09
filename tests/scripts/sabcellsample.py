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

#import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.sabcelleval import RefCell
from NCTestUtils.env import ncsetenv
from NCrystalDev.misc import evaluate_query as ncquery

def testcell(*,E_div_kT, alpha, beta, svals = None, do_plot=False ):
    a1, a2 = alpha
    b1, b2 = beta
    if svals is None:
       svals=[1.0,1.0,1.0,1.0]
       assert len(svals)==4
    vals = [E_div_kT,a1,a2,b1,b2,*svals]
    res = ncquery( [ 'sab','sglcell']+['@%s'%str(e) for e in vals])#fixme: ask for samples also
    refcell = RefCell( a1=a1, a2=a2, b1=b1, b2=b2,
                       s11=svals[0], s12=svals[1],
                       s21=svals[2], s22=svals[3] )
    print("HEJSA")
    mpsamples = refcell.sample( E_div_kT, n=10000, seed=5 )
    print("HEJSA2")

    if do_plot:
        from NCTestUtils.sabcelleval import plot_celleval
        from NCrystalDev.plot import PlotContext
        import numpy as np
        pctx = PlotContext()
        plot_celleval( res, **pctx.kwargs_subcontext() )
        ab = np.asarray(mpsamples,dtype=float).copy()
        print(ab)
        a,b = ab.T[0], ab.T[1]
        pctx.axis.scatter(b,a,marker='.',alpha=0.2)
        pctx.finalise( do_grid = False, do_legend='draggable' )

def main(do_plot):
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')

    #Fixme: work on test points, plots, etc.
    testpts = [
#       dict( E_div_kT=0.06,
#             alpha=(6.6767857524277714,6.75015702442328),
#             beta=(5.4022664171285335,5.5029404418379908) ),
 #       dict(E_div_kT=1.5,alpha=(1.0,2.0),beta=(-2.0,-1.0)),
#       dict(E_div_kT=1.5,alpha=(1.0,2.0),beta=(-2.0,-1.0),
#            svals=[1.0,20.0,3.0,4.0]),
        dict(E_div_kT=0.0001,alpha=(5.5,6.0),beta=(5.2,6.0),
             svals=[1.0,1e4,1e4,1.0]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0,1,0,1]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0)),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0.1,1000,1000,0.1]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[1,1000,1,1]),
        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,4.0),
             svals=[0,0,1,1]),


#        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-3.0,-0.5)),
##        dict(E_div_kT=2.0,alpha=(0.5,12.0),beta=(-0.1,6.0)),
##        dict(E_div_kT=2.0,alpha=(0.0,12.0),beta=(-3.0,10.5)),
##        dict(E_div_kT=2.0,alpha=(0.0,6.0),beta=(-1.0,4.5)),
#        dict(E_div_kT=2.0,alpha=(0.5,3.0),beta=(-10.0,40.5)),
#        dict(E_div_kT=2.0,alpha=(0.5,3.0),beta=(-10.0,40.5),
#             svals=[1.0,20.0,3.0,4.0]),
##        dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(11.0,12.0)),
##        dict(E_div_kT=1e-5,alpha=(11.0,13.0),beta=(0.5,12.0)),
        #dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(0.5,12.0)),#this one?
##        dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(11.0,12.0),
##             svals=[1.0,20.0,3.0,4.0]),
#        dict(E_div_kT=0.01,alpha=(11.35,13),beta=(10,12.0),
#             svals=[0,1,1,0]),
##        dict(E_div_kT=1e-5,alpha=(11.0,13.0),beta=(0.5,12.0),
##             svals=[1.0,20.0,3.0,0]),
##        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
##             svals=[1.0,1e-200,1.0,1e-200]),
##        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
##             svals=[1.0,1e-20,1.0,1e-20]),
##        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
##             svals=[1e-20,1.0,1.0,1e-20]),
##        dict(E_div_kT=1e-5,alpha=(11.0,13.0),beta=(0.5,12.0),
##             svals=[1.0,20.0,3.0,1e-300]),
##        dict(E_div_kT=1e-5,alpha=(0.5,12.0),beta=(0.5,12.0),
##             svals=[1.0,20.0,3.0e-4,4.0]),
##        dict(E_div_kT=0.01,alpha=(5.499,5.5),beta=(5,5.92),
##             svals=[1.0,99,1.0,99]),
##        dict(E_div_kT=0.01,alpha=(5.499,5.5),beta=(5,5.92),
##             svals=[1.0,101,1.0,101]),
##        dict(E_div_kT=0.01,alpha=(5.499,5.5),beta=(5,5.92),
##             svals=[1.0,0.011,1.0,0.011]),
##        #Note: some rather bad results, attributed to numerical issues when a2
##        #~- a1.
##        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
##             svals=[1.0,99,1.0,99]),
##        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
##             svals=[1.0,101,1.0,101]),
##        dict(E_div_kT=0.01,alpha=(5.5-1e-12,5.5),beta=(5,5.92),
##             svals=[1.0,0.011,1.0,0.011]),
##        dict(E_div_kT=0.01,alpha=(5.47*0+5.5-1e-12,5.5),beta=(5,5.92),
##             svals=[1.0,0.009,1.0,0.009]),
##
##
##
###        dict(E_div_kT=1,
###             alpha=(1e-50,1e-10),
###             beta=(-1,1),
###             svals=[1,1,1,1]),
###
##        dict(
##            #E_div_kT=0.2812650004,#FAILS
##            E_div_kT=0.2812650005,#OK
##            alpha=(5e-10,8e-10),
##            beta=(-3e-5,-2e-05),
##            svals=[1,1,1,1]),
##
##
##        dict(E_div_kT=0.000395856,
##             alpha=(8.36433e-50,8.36433e-10),
##             beta=(-2.78659e-05,-2.60257e-08),
##             svals=[2.09324e-51,2.09324e-11,2.09321e-51,2.09321e-11]),
##
    ]

    for i,data in enumerate(testpts):
        i += 1
        #if i not in (29,):
        #    continue
        print()
        print("=============>")
        print("=============> Launching test %i"%i)
        print("=============>")
        print()
        testcell(**data,do_plot=do_plot)


if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:])

#fixme: Add option to test a huge bunch of randomly generated cells.
