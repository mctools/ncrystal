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

from NCrystalDev.misc import evaluate_query as ncquery
from NCTestUtils.mpmathctx import get_mpmath_context

def mp_integrate01_kpowx( k, mp ):
    k = mp.mpf(k)
    if k == 1:
        return mp.mpf(1)
    lnk = mp.log(k)
    return (k-mp.mpf(1))/lnk

def mp_integrate01_xkpowx( k, mp ):
    k = mp.mpf(k)
    if k == 1:
        return mp.mpf('1/2')
    lnk = mp.log(k)
    return (k*(lnk-mp.mpf(1))+mp.mpf(1))/(lnk**2)

def mp_sample_xkpowx( k, R, mp ):
    k=mp.mpf(k)
    R=mp.mpf(R)
    one=mp.mpf(1)
    if k==one:
        return R
    return mp.log(one+R*(k-one))/mp.log(k)

def main():
    mp = get_mpmath_context(1000)#fixme: reduce??
    d = ncquery(['util','mathval','kpowx'])
    kvals = d['k']
    nc_kpowx = d['integral01_kpowx']
    nc_xkpowx = d['integral01_xkpowx']
    nc_kpowx_samples = d['sample_kpowx']
    assert len(kvals)>10 and 1.0 in kvals
    assert max(kvals)>1e250 and min(kvals)<1e-250

    cmps = []
    for k, ncval in zip(kvals,nc_kpowx):
        refval = mp_integrate01_kpowx(k,mp)
        cmps.append( (ncval,refval,'int_0^1(k^x)dx [k=%.15g]'%k) )

    for k, ncval in zip(kvals,nc_xkpowx):
        refval = mp_integrate01_xkpowx(k,mp)
        cmps.append( (ncval,refval,'int_0^1(x*k^x)dx [k=%.15g]'%k) )

    worst = None
    for ncval,refval,descr in cmps:
        rd = abs(ncval/refval-mp.mpf(1))
        print('%s = %.14g [precision: %.2g]'%(descr,float(ncval),float(rd)))
        worst = rd if worst is None else max(worst,rd)

    thr = 5e-15
    if not worst < thr:
        print("Worst precision: %.3g"%float(worst))
        raise SystemExit(f'ERROR: Precision not below {thr:g}!')
    else:
        print(f"Precision < {thr:g}? : YES")

    sample_refvals=[]
    for k in nc_kpowx_samples['k']:
        for R in nc_kpowx_samples['R']:
            sample_refvals.append((k,R,mp_sample_xkpowx( k, R, mp )))

    for krrefval, ncval in zip(sample_refvals,nc_kpowx_samples['samples']):
        k,R,refval=krrefval
        cmps.append( (ncval,refval,
                      'sample k^x on [0,1] [k=%.15g,R=%.15g]'%(k,R)) )

    worst = None
    for ncval,refval,descr in cmps:
        rd = abs(ncval/refval-mp.mpf(1))
        print('%s = %.14g [precision: %.2g]'%(descr,float(ncval),float(rd)))
        worst = rd if worst is None else max(worst,rd)

    thr = 1e-14
    if not worst < thr:
        print("Worst precision (samples): %.3g"%float(worst))
        raise SystemExit(f'ERROR: Sampling precision not below {thr:g}!')
    else:
        print(f"Sampling precision < {thr:g}? : YES")

    assert len(nc_kpowx_samples['k']) > 20
    assert len(nc_kpowx_samples['R']) > 7

if __name__=='__main__':
    main()
