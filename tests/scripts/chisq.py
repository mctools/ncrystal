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

# NEEDS: numpy ase

# Note: declaring ase as a dependency is just a workaround to have access to
# scipy (which ase pulls in), without going to the trouble of adding it directly
# as a general test dependency.

# This tests validates the precision of our chisq_cdf function, in the context
# of it's usage for checking histogram compatibility via p-pvalues.

from NCrystalDev._hist import _chisq_cdf as nc_chisq_cdf
from scipy.stats import chi2 as scipy_chi2
import time

def validate(x,k):
    if x<0.1:
        return
    c = nc_chisq_cdf( x, k )
    c_ref = scipy_chi2.cdf( x, k )
    assert 0.0 <= c <= 1.0
    assert 0.0 <= c_ref <= 1.0
    pval = 1.0 - c
    pval_ref = 1.0 - c_ref

    rdiff_cval = abs(c-c_ref)/(1e-9+c_ref)
    rdiff_pval = abs(pval-pval_ref)/(1e-9+pval_ref)
    rdiff = max(rdiff_cval,rdiff_pval)

    target_prec = 0.1
    if (pval > 0.97 and pval_ref > 0.97) and x < 0.5*k:
        #In this case, no histogram will be found to be incompat, so can loosen.
        target_prec = 0.5

    if rdiff>target_prec:
        print(("chisq_cdf(x=%g, k=%i): ncrystal->%g (=1-%g) vs "
               "ref->%g (=1-%g) (rdiff: %.2g)")%( x, k, c, 1.0-c,
                                                  c_ref, 1.0-c_ref, rdiff ))
        raise SystemExit('ERROR: Not fit for purpose p-value calc.')

def main():
    from NCrystalDev._numpy import _np_geomspace
    xvals = _np_geomspace(1e-5,1e5,40)
    kvals = ( [i for i in range(1,10+1)]
              + [i for i in range(12,200,10)]
              + [i for i in range(201,1000,100)]
              + [i for i in range(1001,10000,1000)]
              + [i for i in range(10001,100000*100,300000)] )

    n = 0
    print("Checking x=0")
    for k in kvals:
        validate(0.0,k)#x=0.0
        n+=1

    for k in kvals:
        print(f"Checking k={k} with {len(xvals)} x values.")
        t0 = time.time()
        for x in xvals:
            validate(x,k)
        t0 = time.time() - t0
        print("   (time per validation: %.2g ms"%((t0/len(xvals))*1000.0))

if __name__ == '__main__':
    main()

