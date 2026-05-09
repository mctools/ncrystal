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

# FIXME: NOT TESTING MUCH RIGHT NOW!! JUST A PLOT SCRIPT. (also: check if in the
# end we need numpy or not).


def samplepb(x0,y0,x1,y1,nsample):
    from NCrystalDev.misc import evaluate_query as ncquery
    import numpy as np
    res = ncquery( [ 'sab','samplepb',
                     str(x0),str(y0),str(x1),str(y1),str(nsample) ] )
    s = res['samples']
    if not s:
        print("sampling not possible here!")
        return np.ones(0),np.ones(0)
    #print("First few samples: %s"%(res['samples'][0:6]))
    a = np.asarray(res['samples'],dtype=float)
    x=a[:,0]
    y=a[:,1]

    return x,y,dict(rng_per_sample=res['rng_per_sample'],
                    sampler_details=res['sampler_details'])

def plot_samples(axis,x0,y0,x1,y1,nsample, **plotkwargs):
    x, y, details = samplepb(x0=x0,y0=y0,x1=x1,y1=y1,
                             nsample=nsample)
    axis.scatter(x, y, **plotkwargs)
    return details

def plot_pb(axis,xmax,n = 10000, **plotkwargs):
    from NCrystalDev._numpy import _np_linspace
    import numpy as np
    x = _np_linspace( 0.0, xmax, n )
    sx = np.sqrt(x)
    ym = (sx-1.0)**2
    yp = (sx+1.0)**2
    _ = axis.plot( x, ym, **plotkwargs )
    if 'color' not in plotkwargs:
        plotkwargs['color']=_[0].get_color()
    axis.plot( x, yp, **plotkwargs )

def plot_box(axis,x0,y0,x1,y1,**plotkwargs):
    axis.plot( [x0,x0,x1,x1,x0],
               [y0,y1,y1,y0,y0], **plotkwargs )

def show(x0,y0,x1,y1):
    import matplotlib.pyplot as plt
    box = dict(x0=x0,y0=y0,x1=x1,y1=y1)
    fig,axis = plt.subplots()
    xmax = box['x1']*1.2
    plot_pb( axis, xmax = xmax )
    plot_box( axis, **box )
    si = plot_samples( axis, **box, nsample=5000,
                            alpha=0.4, marker='.' )
    rs = si['rng_per_sample']
    approx_ar = 2.0/rs
    axis.set_title("Acceptance rate ~= %.3g%% (RNG/sample=%.3g, %s)"
                   %(approx_ar*100,
                     rs,
                     si['sampler_details']['samplemode']))

    orc_x, orc_y = zip(*si['sampler_details']['overlay_region_curve'])
    #print(si['sampler_details']['overlay_region_curve'])
    axis.plot(orc_x,orc_y,ls='--',alpha=0.5,lw=2)
    axis.set_xlim(0.0)
    axis.set_ylim(min(0.0,min(orc_y)))
    plt.show()

def main():
    return#fixme
    boxes = [ dict(x0=0.9,y0=0.8,x1=1.1,y1=1.1),
              dict(x0=0,y0=0.1,x1=2.5,y1=1.1),
              dict(x0=0,y0=0.0,x1=2.5,y1=1.1),
              dict(x0=0.1,y0=0.1,x1=2.2,y1=60.0),
              dict(x0=0.1,y0=0.1,x1=0.4,y1=0.5),
              dict(x0=0.1,y0=1.01,x1=0.4,y1=20.5),
              dict(x0=10000,y0=9999,x1=1000000,y1=1000000.5),
              dict(x0=0,y0=0,x1=4.9,y1=12),
              dict(x0=0,y0=0,x1=5.1,y1=12),
              dict(x0=0,y0=3,x1=5.1,y1=12),
              dict(x0=0,y0=0,x1=5.1,y1=0.1),
              dict(x0=1.0,y0=0,x1=26,y1=30),
              dict(x0=1.0,y0=0,x1=26,y1=17.3)]
    for b in boxes:
        show(**b)

if __name__=='__main__':
    main()
