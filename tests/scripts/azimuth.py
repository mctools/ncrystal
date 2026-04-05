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
import NCrystalDev.core as nccore
from NCrystalDev.constants import ekin2wl
from NCrystalDev.hist import HistFiller1D
from NCTestUtils.env import ncsetenv
import numpy as np


def load_sglbragg_scatter():
    mat = nccore.load('stdlib::Al_sg225.ncmat;dcutoff=2.2;comp=bragg')
    # Al:
    #  H   K   L  d_hkl[Aa] Mult. FSquared[barn]
    #  1   1   1    2.33803    8     1.77316
    #  2   0   0    2.02479    6     1.73179
    #  2   2   0    1.43174   12     1.57574
    #  3   1   1    1.22099   24     1.46799
    mat.dump()
    hkl = list(mat.info.hklObjects())
    assert len(hkl)==1
    assert mat.scatter.getName() == 'PowderBragg'
    return hkl[0].dspacing, mat.scatter

def calc_mu( indir, outdir_arrays ):
    ix,iy,iz = indir
    assert isinstance(ix,float)
    ox,oy,oz = outdir_arrays
    assert hasattr(ox,'__len__')
    return ox*ix + oy*iy + oz*iz

def calc_dotp( v1, v2 ):
    #Works both for tuple of floats and tuples of arrays
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def normvect( a, b, c ):
    n = np.sqrt(a*a+b*b+c*c)
    assert n>0.0
    m = 1.0/n
    return (a*m,b*m,c*m)

def crossvect_then_norm( a, b, norm=True ):
    a1,a2,a3=a
    b1,b2,b3=b
    v = ( a2*b3-a3*b2,
          a3*b1-a1*b3,
          a1*b2-a2*b1 )
    return normvect(*v)

def calc_azimuthal( indir, outdir_arrays ):
    #We need to define vectors orthogonal to indir.
    indir = normvect(*indir)
    if abs(calc_dotp( indir, (1.,0.,0.)))>0.9:
        u = ( 0.0, 1.0, 0.0 ) #indir ~parallel with x-axis
    else:
        u = ( 1.0, 0.0, 0.0 )
    v = crossvect_then_norm( indir, u )
    u = crossvect_then_norm( indir, v )
    #Ok, we now have orthonormal basis (indir,u,v):
    assert abs(calc_dotp(u,indir))<1e-7
    assert abs(calc_dotp(v,indir))<1e-7
    assert abs(calc_dotp(u,v))<1e-7
    #Project outdir onto the (u,v) plane:
    t = calc_dotp(outdir_arrays,indir)
    ot = ( outdir_arrays[0]-t*indir[0],
           outdir_arrays[1]-t*indir[1],
           outdir_arrays[2]-t*indir[2] )
    #finally get the azimuthal angle:
    return np.arctan2(calc_dotp(ot,v),calc_dotp(ot,u))

def main(do_lux, do_plot):
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')
    dsp, scatterproc = load_sglbragg_scatter()
    print("Loaded material with single dspacing %g Aa"%dsp)
    indirs = [ (0,0,1),(1,0,0),(1,1,1) ]
    if do_lux:
        indirs+= [(0,1,0),(1,1,0)]
    for i in range(len(indirs)):
        ux,uy,uz = indirs[i]
        m = (ux*ux+uy*uy+uz*uz)**0.5
        assert m>0.0
        indirs[i] = (ux/m, uy/m, uz/m)

    nrepeat=int(1e8) if do_lux else int(4e6)
    nbins=200 if do_lux else 100
    for sinthetabragg, indir in zip([ 0.05, 0.5, 0.95 ],indirs):
        print("Checking sinthetabragg=%g indir=(%g,%g,%g)"
              %(sinthetabragg,*indir))
        wl = 2*dsp*sinthetabragg
        mu_expected = 1.0 - 2.0 * sinthetabragg**2


        def sample_and_check_n(thehist,n):
            ekin_final, outdir = scatterproc.scatter(direction=indir, wl=wl,
                                                     repeat=n)
            #Check energy:
            assert (abs(ekin2wl(ekin_final)-wl)<1e-7).all()
            #Check mu:
            assert (abs(calc_mu(indir,outdir)-mu_expected)<1e-7).all()
            #Check azimuthal:
            phi = calc_azimuthal( indir, outdir )
            assert not np.isnan(phi).any()
            assert -np.pi < phi.min() < -0.9*np.pi
            assert np.pi*0.9 < phi.max() < np.pi
            thehist.fill( phi )

        h = HistFiller1D(nbins,-np.pi,np.pi,'phi')
        #keep mem usage low:
        n_at_once = int(1e5)
        assert nrepeat%n_at_once==0
        for i in range(nrepeat//n_at_once):
            sample_and_check_n( h, n_at_once )

        h = h.to_hist1d()

        if do_plot:
            h.plot()

        #Let us check that the distribution is indeed flat:
        bin_av = nrepeat/nbins
        bin_stddev = np.sqrt(bin_av)
        dev = (h.contents-bin_av)/bin_stddev
        assert abs(dev.sum()/len(dev)) < 1e-10#by construction this should be
                                              #zero for any distribution
        print(dev.min(),dev.max())
        assert -5.0 < dev.min() < -1.0
        assert 1.0 < dev.max() < 5.0

if __name__ == '__main__':
    import sys
    main( do_lux = '--lux' in sys.argv[1:],
          do_plot = '--plot' in sys.argv[1:] )
