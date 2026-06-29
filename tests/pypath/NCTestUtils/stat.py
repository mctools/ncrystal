
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

# A few statistics utilities.

def kolmogorov_smirnov_dist(x, y):
    """Accepts two arrays of sampled values, and constructs the CDF of
    both. Returns the maximum distance between the two CDFs.
    """
    import numpy as np
    x, y = np.asarray(x,dtype=float), np.asarray(y,dtype=float)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        raise RuntimeError("KS distance requires non-empty samples.")

    # Find positions where values cross:
    x, y = np.sort(x), np.sort(y)
    idx_in_y = np.searchsorted(y, x, side='right')
    idx_in_x = np.searchsorted(x, y, side='right')
    del x, y

    # CDF differences at x points:
    tmp = np.arange(1, nx + 1, dtype = float )
    tmp *= (1.0/nx)
    tmp -= idx_in_y * (1.0/ny)
    distx = max(tmp.max(),-tmp.min())#like np.max(np.abs(t)) without temporary arr.
    del idx_in_y

    # CDF differences at y points:
    tmp = np.arange(1, ny + 1, dtype = float )
    tmp *=  ( 1.0 / ny )
    tmp -= idx_in_x * ( 1.0 / nx)
    return max(distx,tmp.max(),-tmp.min())


def kolmogorov_smirnov_pvalue( x, y ):
    """Based on the Kolmogorov-Smirnov distance, find the p-value for the two 1D
       arrays of sampled values to have been sampled from the same underlying
       distribution
    """

    import numpy as np
    x, y = np.asarray(x,dtype=float), np.asarray(y,dtype=float)
    assert len(x.shape)==1
    assert len(y.shape)==1
    n, m = x.size, y.size
    if n == 0 or m == 0:
        raise RuntimeError("KS test requires non-empty samples.")

    #Construct cumulative distributions for the two data sets, evaluated at any
    #pt in either x or y:
    D =  kolmogorov_smirnov_dist(x, y)
    del x, y

    #Let us now convert this distance to a p-value (strictly speaking a
    #asymptotic two-sided Kolmogorov distribution p-value):
    lam = D * np.sqrt(n * m / (n + m))
    if not lam > 0.0:
        return 1.0#no distance, perfect match
    m2lam2, s, sgn = (-2.0 * lam * lam), 0.0, -1
    # p/2 = sum_{k>=1} (-1)^{k-1} exp(-2 k^2 lam^2)
    for k in range(1, 201):
        sgn *= -1
        e = np.exp( m2lam2*k*k )
        s += sgn*e
        if k>20 and e < 1e-10*abs(s):
            break
    return min(1.0,max(0.0,float(2.0 * s)))
