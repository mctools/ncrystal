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

from NCrystalDev import NCBadInput, NCCalcError
from NCrystalDev._hist import Hist1D
from NCTestUtils.hists import Hist1Dcpp
from NCTestUtils.randutils import TestRNG
from NCTestUtils.common import ensure_error
import math

def main():

    h = Hist1Dcpp( 5, 0.0, 100.0, clamp_overflows=True )
    h.fill( 12.34 )
    print("JSON:",h.toJSON())
    h.toPyHist().dump()

    h = Hist1Dcpp( 5, 0.0, 100.0, clamp_overflows=False )
    assert h.toPyHist().empty
    h.fill( 10.0 )
    assert not h.toPyHist().empty
    h.fill( 90.0 )
    h.toPyHist().dump()

    h = Hist1Dcpp( 5, 0.0, 100.0, clamp_overflows=False )
    h.fill( -10.0 )
    h.fill( 110.0 )
    h.fill( 1e99 )
    h.toPyHist().dump()

    h = Hist1Dcpp( 5, 0.0, 100.0, clamp_overflows=True )
    h.fill( 110.0 )
    h.toPyHist().dump()

    h = Hist1Dcpp( 5, 0.0, 100.0, clamp_overflows=False )
    h.fill( -10.0, w=0.01 )
    h.fill( 50.0, w=0.01 )
    h.fill( 110.0, w=0.01 )
    h.fill( 1e99, w=0.01 )
    h.toPyHist().dump()

    h1 = Hist1Dcpp( 5, 0.0, 100.0, clamp_overflows=False )
    h1.fill( -10.0, w=0.01 )
    h1.fill( 50.0, w=0.01 )
    h1.fill( 110.0, w=0.01 )
    h1 = h1.toPyHist()

    h2 = Hist1Dcpp( 5, 0.0, 100.0, clamp_overflows=False )
    h2.fill( -10.0, w=0.01 )
    h2.fill( 50.0, w=0.01 )
    h2.fill( 110.0, w=0.01 )
    h2 = h2.toPyHist()
    h2.rebin(1)

    h1.add_contents(h2.clone())
    h1.dump()
    h1.dump(contents = False)

    print( 'h1.bindata: ', h1.bindata(json_compat=True) )
    print( 'h1.stats: ', h1.stats )
    h1.set_title("My Title")
    print( 'h1.to_dict: ', h1.to_dict(json_compat=True) )
    print( 'h1.to_json: ', h1.to_json() )

    pval = h1.check_compat(h1,return_pval=True)
    print("p-value: %g"%pval)
    assert pval==1.0
    h1.check_compat(h1,check=True)
    assert h1.check_compat(h1,return_pval=True,check=True) == 1.0
    assert h1.check_compat(h1) is True

    def testinteg_bins( h, fillws, *a ):
        v,e = h.integrate_bins(*a)
        v_expect = math.fsum(fillws) if fillws else 0.0
        e_expect = math.sqrt(math.fsum([e**2
                                        for e in fillws])) if fillws else 0.0
        print("Integrate (expects %g +- %g): %g +- %g"%(v_expect,e_expect,v,e))
        assert abs(v-v_expect) < 1e-10
        assert abs(e-e_expect) < 1e-10

    def testinteg( h, fillws, *a ):
        v,e = h.integrate(*a)
        v_expect = math.fsum(fillws) if fillws else 0.0
        e_expect = math.sqrt(math.fsum([e**2
                                        for e in fillws])) if fillws else 0.0
        print("Integrate (expects %g +- %g): %g +- %g"%(v_expect,e_expect,v,e))
        assert abs(v-v_expect) < 1e-10
        assert abs(e-e_expect) < 1e-10

    h1 = Hist1Dcpp( 20, -10.0, 10.0, clamp_overflows=True )
    h1.fill( -200.0, w=0.001)
    h1.fill( -9.5, w=0.05 )
    h1.fill( -5.5, w=0.01 )
    h1.fill( -0.5, w=0.1 )
    h1.fill( -0.7, w=0.15 )
    h1.fill( 3.5, w=1.0 )
    h1.fill( 4.5, w=0.5 )
    h1.fill( 8.5, w=0.04)
    h1.fill( 9.5, w=0.03)
    h1.fill( 200.0, w=0.02)
    h1 = h1.toPyHist()
    h1.set_title('Foo title')
    h1.dump()
    h1.clone(rebin_factor=4).dump()

    assert h1.nbins == 20
    testinteg_bins(h1,[0.03,0.02],h1.nbins-1,h1.nbins)
    testinteg_bins(h1,[0.04,0.03,0.02],h1.nbins-2,h1.nbins)
    testinteg_bins(h1,[0.04],h1.nbins-2,h1.nbins-1)
    testinteg_bins(h1,[0.04,0.03,0.02],h1.nbins-2,None)
    testinteg_bins(h1,[0.001,0.05],0,1)
    testinteg_bins(h1,[0.001,0.05],None,1)
    testinteg_bins(h1,[],None,0)
    testinteg_bins(h1,[],0,0)

    testinteg(h1,[1.0,0.5],3.0,5.0)
    testinteg(h1,[1.0],2.0,4.0)
    testinteg(h1,[0.5],4.0,5.0)
    testinteg(h1,[0.03,0.02],9.0,10.0)
    testinteg(h1,[0.001,0.05,0.01],None,-3.0)
    testinteg(h1,[0.01,0.1,0.15],-7.0,3.0)
    testinteg(h1,[0.001,0.05,0.01],-40.0,-3.0)
    testinteg(h1,[1.0,0.5],0.0,7.0)
    testinteg(h1,[1.0,0.5,0.04,0.03,0.02],0.0,None)
    testinteg(h1,[1.0,0.5,0.04,0.03,0.02],0.0,10.0)
    testinteg(h1,[1.0,0.5,0.04,0.03,0.02],0.0,120.0)
    testinteg(h1,[],3.0000000001,3.00000001)
    testinteg(h1,[0.001,0.05,0.01,0.1,0.15,1.0,0.5,0.04,0.03,0.02],None,None)

    with ensure_error(NCBadInput,
                      'Invalid integration range requested.'):
        testinteg(h1,[1.0,0.5],5.0,3.0)

    with ensure_error(NCBadInput,
                      'Invalid bin range requested.'):
        testinteg_bins(h1,[1.0,0.5],3,2)

    with ensure_error(NCBadInput,
                      'Value 3.01 does not correspond exactly to a bin'
                      ' edge within the tolerance.'):
        testinteg(h1,[1.0,0.5],3.01,5.0)

    h1 = Hist1Dcpp( 20, -10.0, 10.0, clamp_overflows=False )
    h1.fill( -200.0, w=0.001)
    h1.fill( -9.5, w=0.05 )
    h1.fill( -5.5, w=0.01 )
    h1.fill( -0.5, w=0.1 )
    h1.fill( -0.7, w=0.15 )
    h1.fill( 3.5, w=1.0 )
    h1.fill( 4.5, w=0.5 )
    h1.fill( 8.5, w=0.04)
    h1.fill( 9.5, w=0.03)
    h1.fill( 200.0, w=0.02)
    h1 = h1.toPyHist()
    h1.set_title('same but not clamping overflows')
    h1.dump()
    h1.clone(rebin_factor=4).dump()

    assert h1.nbins == 20
    testinteg_bins(h1,[0.03],h1.nbins-1,h1.nbins)
    testinteg_bins(h1,[0.04,0.03],h1.nbins-2,h1.nbins)
    testinteg_bins(h1,[0.04],h1.nbins-2,h1.nbins-1)
    testinteg_bins(h1,[0.04,0.03,0.02],h1.nbins-2,None)
    testinteg_bins(h1,[0.05],0,1)
    testinteg_bins(h1,[0.001,0.05],None,1)
    testinteg_bins(h1,[0.001],None,0)
    testinteg_bins(h1,[],0,0)

    testinteg(h1,[1.0,0.5],3.0,5.0)
    testinteg(h1,[1.0],2.0,4.0)
    testinteg(h1,[0.5],4.0,5.0)
    testinteg(h1,[0.03],9.0,10.0)
    testinteg(h1,[0.001,0.05,0.01],None,-3.0)
    testinteg(h1,[0.01,0.1,0.15],-7.0,3.0)
    testinteg(h1,[0.05,0.01],-40.0,-3.0)
    testinteg(h1,[1.0,0.5],0.0,7.0)
    testinteg(h1,[1.0,0.5,0.04,0.03,0.02],0.0,None)
    testinteg(h1,[1.0,0.5,0.04,0.03],0.0,10.0)
    testinteg(h1,[1.0,0.5,0.04,0.03],0.0,120.0)
    testinteg(h1,[],3.0000000001,3.00000001)
    testinteg(h1,[0.001,0.05,0.01,0.1,0.15,1.0,0.5,0.04,0.03,0.02],None,None)

    h1b = h1.clone(rebin_factor=4)
    h1b.dump()
    print( h1b._hist_curve() )
    print( h1b.errorbar_args() )
    print( h1b.bar_args() )

    h2 = Hist1Dcpp( 20, -10.0, 10.0, clamp_overflows=False ).toPyHist()
    h2.set_title("Target of add_contents")
    h2.dump()
    h1.dump()
    h2.add_contents(h1)
    h2.dump()
    assert h2.stats == h1.stats
    assert h2.bindata(json_compat=True) == h1.bindata(json_compat=True)
    h2 = Hist1Dcpp( *h1.binning, clamp_overflows=bool(not h1.has_flow) ).toPyHist()
    h2.set_title("Target of add_contents2")
    h2.adopt_contents(h1,keep_title=False)
    assert h2.to_json() == h1.to_json()

    h3 = Hist1Dcpp( 20, -10.0, 10.0,
                    allow_weights = False, clamp_overflows=False )
    h3.fill( -200.0 )
    h3.fill( -9.5 )
    h3.fill( -5.5 )
    h3.fill( -0.5 )
    h3.fill( -0.7 )
    h3.fill( 3.5 )
    h3.fill( 4.5 )
    h3.fill( 8.5 )
    h3.fill( 9.5 )
    h3.fill( 200.0 )
    h3 = h3.toPyHist()
    h3.dump()

    testinteg_bins(h3,[1,]*1,4,5)
    testinteg_bins(h3,[1,]*2,9,10)
    testinteg_bins(h3,[1,]*3,4,10)
    testinteg_bins(h3,[1,]*1,19,20)
    testinteg_bins(h3,[1,]*1,20,None)
    testinteg(h3,[1,]*10,None,None)
    testinteg_bins(h3,[1,]*1,None,0)
    testinteg_bins(h3,[1,]*2,None,1)
    testinteg_bins(h3,[1,]*2,None,2)
    testinteg_bins(h3,[1,]*5,None,11)
    testinteg_bins(h3,[1,]*7,9,None)
    testinteg(h3,[1,]*10,None,None)
    testinteg(h3,[1,]*8,-10.0,10.0)
    testinteg(h3,[1,]*9,None,10.0)
    testinteg(h3,[1,]*9,-10.0,None)
    testinteg(h3,[1,]*5,-5.0000000001,8.9999999999999)

    tmp = h3.to_json()
    h4 = Hist1Dcpp( *h3.binning,
                    allow_weights = False, clamp_overflows=False ).toPyHist()
    h3.add_contents(h4)
    assert tmp == h3.to_json()

    h4 = Hist1Dcpp( *h3.binning,
                    allow_weights = False, clamp_overflows=True )
    h4.fill( -200.0 )
    h4.fill( -9.5 )
    h4.fill( -5.5 )
    h4.fill( -0.5 )
    h4.fill( -0.7 )
    h4.fill( 3.5 )
    h4.fill( 4.5 )
    h4.fill( 8.5 )
    h4.fill( 9.5 )
    h4.fill( 200.0 )
    h4 = h4.toPyHist()
    tmp = h4.to_json()
    h4b = Hist1Dcpp( *h3.binning,
                     allow_weights = True, clamp_overflows=True ).toPyHist()
    h4.add_contents(h4b)
    assert tmp == h4.to_json()

    class FakeAxis:
        def bar( self, **kwargs ):
            print("SPY axis.bar called with kwargs: %s"%kwargs)
        def errorbar( self, **kwargs ):
            print("SPY axis.errorbar called with kwargs: %s"%str(kwargs))
        def set_xlim( self, *a ):
            print("SPY axis.set_xlim called with args: %s"%str(a))

    class FakePLT:
        def gca(self):
            print("SPY plt.gca() called")
            return FakeAxis()
        def show(self):
            print("SPY plt.show() called")

    h3.plot( plt=FakePLT() )

    h5 = Hist1Dcpp( 6, -3.0, 3.0,
                    allow_weights = True, clamp_overflows=True )
    h5.fill( -0.5 )
    h5.fill( 1.5 )
    h5 = h5.toPyHist()
    h5.dump()
    h6 = Hist1Dcpp( 6, -3.0, 3.0,
                    allow_weights = False, clamp_overflows=True )
    h6.fill( -200.0 )
    h6.fill( -0.5 )
    h6 = h6.toPyHist()
    h6.dump()

    h5.add_contents(h6)
    h5.dump()

    h7 = Hist1Dcpp( 6, -3.0, 3.0,
                    allow_weights = True, clamp_overflows=True )
    h7.fill( -0.5 )
    h7.fill( 1.5 )
    h7 = h7.toPyHist()
    h7.dump()
    h8 = Hist1Dcpp( 6, -3.0, 3.0,
                    allow_weights = False, clamp_overflows=True )
    h8.fill( 200.0 )
    h8 = h8.toPyHist()
    h8.dump()

    h7.add_contents(h8)
    h7.dump()
    h7.dump(highres=True)

    def verify_stat(hist,vals,weights = None):
        if hasattr(hist,'toPyHist'):
            hist = hist.toPyHist()
        if weights is None:
            weights = [1,]*len(vals)
        sumw = math.fsum(weights)
        sumwx = math.fsum( [w*x for x,w in zip(vals,weights)] )
        sumwx2 = math.fsum( [w*x**2 for x,w in zip(vals,weights)] )
        assert sumw > 0.0
        assert sumwx2 > 0.0
        mean = sumwx/sumw
        rms2 = sumwx2 / sumw - (mean)**2
        assert rms2 >= 0.0
        rms = math.sqrt(rms2)

        ref = dict(integral=sumw, mean=mean, rms=rms)
        hvals = dict( integral=hist.integral,
                      mean=hist.mean,
                      rms=hist.rms )
        okall = True
        for k in ref.keys():
            v, vref = hvals[k], ref[k]
            ok = abs(v-vref)<1e-10
            print("  -> Stat: %s = %g (ref) %g (hist) %s"%( k, vref, v,
                                                            'OK' if ok
                                                            else 'BAD' ) )
            okall = ok and okall
        assert okall

    #test RMS merging:
    l1 = [ 1.0, 3.0, 5.0 ]
    l2 = [ -100.0, 3.0, 5.0, 2.0, 1.2 ]
    w1 = [ 0.1, 0.7, 0.3 ]
    w2 = [ 0.6, 0.2, 0.4, 0.5, 2.0 ]

    for clamp in [True,False]:
        for allow_weights in [True,False]:
            histkw = dict( nbins=6,xmin=-3.0,xmax=3.0,
                           allow_weights = allow_weights,
                           clamp_overflows = clamp )
            print()
            print('-'*80)
            print('-'*80)
            print(histkw)
            h_1 = Hist1Dcpp(**histkw)
            w = w1 if allow_weights else None
            if allow_weights:
                for x,ww in zip(l1,w):
                    h_1.fill( x,ww )
            else:
                for x in l1:
                    h_1.fill( x )
            verify_stat( h_1, l1, weights = w)
            h_2 = Hist1Dcpp(**histkw)
            w = w2 if allow_weights else None
            if allow_weights:
                for x,ww in zip(l2,w):
                    h_2.fill( x,ww )
            else:
                for x in l2:
                    h_2.fill( x )
            verify_stat( h_2, l2, weights = w)
            h_3 = h_1.toPyHist()
            h_3.add_contents( h_2.toPyHist() )
            h_3.dump()
            w = w1+w2 if allow_weights else None
            verify_stat( h_3, l1+l2, weights = w)



    #compat:
    rng = TestRNG()
    h = Hist1Dcpp( 100, 0.0, 100.0 )
    h_compat = Hist1Dcpp( 100, 0.0, 100.0 )
    h_slightcompat = Hist1Dcpp( 100, 0.0, 100.0 )
    h_incompat = Hist1Dcpp( 100, 0.0, 100.0 )
    for i in range(10000):
        h.fill( rng.gauss(mu=40.0,sigma=15.0), w=0.001+rng.rand01() )
    for i in range(8000):
        h_compat.fill( rng.gauss(mu=40.0,sigma=15.0), w=0.001+rng.rand01() )
    for i in range(8000):
        h_slightcompat.fill( rng.gauss(mu=39.5,sigma=15.5), w=0.001+rng.rand01() )
    for i in range(10000):
        h_incompat.fill( rng.gauss(mu=50.0,sigma=15.0), w=0.001+rng.rand01() )
    h = h.toPyHist()
    h_compat = h_compat.toPyHist()
    h_slightcompat = h_slightcompat.toPyHist()
    h_incompat = h_incompat.toPyHist()
    if False:
        [_.plot() for _ in (h, h_compat, h_slightcompat, h_incompat) ]
    assert h.check_compat(h_compat)
    pval = h.check_compat(h_compat,return_pval=True,check=True)
    print("h vs. h_compat p-value: %g"%pval)

    assert h.check_compat(h_slightcompat)
    pval = h.check_compat(h_slightcompat,return_pval=True,check=True)
    print("h vs. h_slightcompat p-value: %g"%pval)

    assert not h.check_compat(h_incompat)
    with ensure_error(NCCalcError, 'check_compat failed:'
                      ' p-value is not greater than 0.05.'):
        h.check_compat(h_incompat,check=True)

    pval = h.check_compat(h_incompat,return_pval=True)
    print("h vs. h_incompat p-value: %.10f"%pval)
    assert pval < 1e-13

    with ensure_error(NCBadInput,
                      'chisquare_dist: incompatible histogram integrals.'):
        h.check_compat(h_compat,force_norm=False)

    with ensure_error(NCBadInput,'chisquare_dist: incompatible binnings.'):
        _=Hist1Dcpp( 100, 0.0, 100.0 ).toPyHist()
        _.check_compat(Hist1Dcpp( 100, 0.0, 99.0 ).toPyHist())

    with ensure_error(NCBadInput,'chisquare_dist: '
                      'incompatible overflow settings.'):
        _=Hist1Dcpp(10,0.0,1.0,clamp_overflows=True).toPyHist()
        _.check_compat(Hist1Dcpp(10,0.0,1.0,clamp_overflows=False).toPyHist())


    h1=Hist1Dcpp(10,0.0,1.0,clamp_overflows=False)
    h1.fill(-99.0)
    h1.fill(-99.0)
    h2=Hist1Dcpp(10,0.0,1.0,clamp_overflows=False)
    h2.fill(-9999999999.0)
    h1.toPyHist().check_compat(h2.toPyHist(),check=True)

    with ensure_error(NCBadInput,'histogram scale factor out of range.'):
        Hist1Dcpp(10,0.0,1.0).toPyHist().scale(0.0)
    with ensure_error(NCBadInput,'histogram scale factor out of range.'):
        Hist1Dcpp(10,0.0,1.0).toPyHist().scale(-1.0)
    with ensure_error(NCBadInput,'histogram scale factor out of range.'):
        Hist1Dcpp(10,0.0,1.0).toPyHist().scale(1.0001e150)

    h1=Hist1Dcpp(5,0.0,1.0,allow_weights=False)
    h1.fill(-20.0)
    h1.fill(0.1)
    h1.fill(0.7)
    h1.fill(0.7)
    h1.fill(0.7)
    h1.fill(11.0)
    h1 = h1.toPyHist()
    h1.dump()
    h1.scale(10.0)
    h1.dump()

    with ensure_error(NCBadInput,'incompatible binning'
                      ' (100,0.1,100.1) vs. (100,0.1,99.1).'):
        _=Hist1Dcpp( 100, 0.1, 100.1 ).toPyHist()
        _.add_contents(Hist1Dcpp( 100, 0.1, 99.1 ).toPyHist())

    with ensure_error(NCBadInput,'incompatible under/overflow settings.'):
        _=Hist1Dcpp( 100, 0.1, 100.1, clamp_overflows=False ).toPyHist()
        _.add_contents(Hist1Dcpp( 100, 0.1, 100.1,
                                  clamp_overflows=True ).toPyHist())

    d = h1.to_dict()
    Hist1D(d).dump()

    d = h1.to_dict()
    del d['bindata']['overflow']
    with ensure_error(NCBadInput,'Inconsistent input: both under and '
                      'overflow must be available if one of them is.'):
        Hist1D(d)

    d = h1.to_dict()
    del d['bindata']['overflow_errorsq']
    with ensure_error(NCBadInput,'Inconsistent input: both under and overflow'
                      ' error^2 info must be available if one of them is.'):
        Hist1D(d)

    d = h1.to_dict()
    del d['bindata']['underflow']
    with ensure_error(NCBadInput,'Inconsistent input: both under and '
                      'overflow must be available if one of them is.'):
        Hist1D(d)

    d = h1.to_dict()
    del d['bindata']['underflow_errorsq']
    with ensure_error(NCBadInput,'Inconsistent input: both under and overflow'
                      ' error^2 info must be available if one of them is.'):
        Hist1D(d)

    for n in ['rms','mean','minfilled','maxfilled']:
        d = h1.to_dict()
        del d['stats'][n]
        with ensure_error(NCBadInput,'Inconsistent input: all of '
                          'rms/mean/minfilled/maxfilled stats'
                          ' must be available if any is.'):
            Hist1D(d)

    d = h1.to_dict()
    del d['stats']['integral']
    assert Hist1D(d).to_json() == h1.to_json()

    d = h1.to_dict()
    del d['stats']
    Hist1D(d).add_contents(h1)

    h1=Hist1Dcpp(5,0.0,1.0,allow_weights=True)
    h1.fill(-20.0)
    h1.fill(0.1)
    h1.fill(0.7)
    h1.fill(0.7)
    h1.fill(0.7)
    h1.fill(11.0)
    h1 = h1.toPyHist()
    h1.dump()

    h2=Hist1Dcpp(5,0.0,1.0,allow_weights=True)
    h2.fill(-20.0,w=1.0+1e-14)
    h2.fill(0.1,w=1.0+1e-14)
    h2.fill(0.7)
    h2.fill(0.7)
    h2.fill(0.7)
    h2.fill(11.0,w=1.0-2e-14)
    h2 = h2.toPyHist()
    h2.dump()

    h1.add_contents(h2)
    h1.dump()

    h1=Hist1Dcpp(5,0.0,1.0,allow_weights=False)
    h1.fill(-20.0)
    h1.fill(0.1)
    h1.fill(0.7)
    h1.fill(0.7)
    h1.fill(0.7)
    h1.fill(11.0)
    h1 = h1.toPyHist()
    h1.dump()

    h2=Hist1Dcpp(5,0.0,1.0,allow_weights=True)
    h2.fill(-20.0,w=1.0+1e-14)
    h2.fill(0.1,w=1.0+1e-14)
    h2.fill(0.7)
    h2.fill(0.7)
    h2.fill(0.7)
    h2.fill(11.0,w=1.0-2e-14)
    h2 = h2.toPyHist()
    h2.dump()

    h1.add_contents(h2)
    h1.dump()


if __name__ == '__main__':
    main()
