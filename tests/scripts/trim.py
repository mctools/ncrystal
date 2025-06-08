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
#
import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev._sabutils import trim_knl_edges as nc_trim_knl
import numpy as np
import copy

def printknl( knl ):
    a,b,s = knl['alpha'], knl['beta'], knl['sab']
    na,nb = len(a), len(b)
    fmt = { 'float' : lambda x : '%5g'%x }
    print(f"Kernel {na} x {nb}:")
    print(f"   alpha: {np.array2string(a,formatter=fmt)}")
    print(f"   beta : {np.array2string(b,formatter=fmt)}")
    print("   sab  : ",end='')
    print( np.array2string( s,
                            formatter=fmt,
                            prefix='          ' ).lstrip() )

def trimknl( knl ):
    print("trimming...")
    a = knl['alpha']
    b = knl['beta']
    s = knl['sab']
    a,b,s = nc_trim_knl( alphagrid = a, betagrid = b, sab = s )
    knl['alpha'] = a
    knl['beta'] = b
    knl['sab'] = s

def testsimple():
    a = np.asarray( [ 0.0, 0.01, 0.2, 3.0, 40.0 ], dtype = float )
    b = np.asarray( [ -10.0, -2.0, 0.3, 4.0 ], dtype = float )
    orig = dict(
        alpha = a,
        beta = b,
        sab = np.ones( len(a) * len(b) ).reshape( ( len(b), len(a) ) )
    )
    printknl(orig)
    knl2 = copy.deepcopy(orig)
    knl2['sab'][:, 1] = 0.0
    printknl(knl2)
    trimknl(knl2)
    printknl(knl2)
    knl2['sab'][:, 0] = 0.0
    printknl(knl2)
    trimknl(knl2)#noop, low alpha should not get trimmed
    printknl(knl2)

    knl2['sab'][0, :] = 0.0
    printknl(knl2)
    trimknl(knl2)
    printknl(knl2)


    knl2['sab'][-1, :] = 0.0
    knl2['sab'][-2, :] = 0.0
    printknl(knl2)
    trimknl(knl2)
    printknl(knl2)
    del knl2

    knl3 = copy.deepcopy(orig)
    printknl(knl3)
    knl3['sab'][0:2, :] = 0.0
    printknl(knl3)
    trimknl(knl3)
    printknl(knl3)
    del knl3

    knl4 = copy.deepcopy(orig)
    printknl(knl4)
    knl4['sab'][-1, :] = 0.0
    printknl(knl4)
    trimknl(knl4)
    printknl(knl4)
    knl4['sab'][1, :] = 0.0
    printknl(knl4)
    trimknl(knl4)
    printknl(knl4)

def _apply_smin( knl, smin ):
    s = knl['sab']
    s[s<=smin] = 0.0

def _plot_sign( knl, **kwargs ):
    from NCrystalDev.plot import plot_knl
    k = copy.deepcopy(knl)
    s = k['sab']
    s[s>0.0] = 1.0
    s[s<0.0] = -1.0
    plot_knl(k,**kwargs)

def _set_upper_alpha_edge_to_zero( knl,nbins ):
    a = knl['alpha']
    b = knl['beta']
    s = knl['sab']
    s = s.reshape( ( len(b), len(a) ) )
    s[:,-nbins:] = 0.0
    knl['sab'] = s.reshape( (len(a)*len(b),) )

def plot_testsab():
    from NCrystalDev.core import createInfo
    from NCrystalDev.plot import plot_knl
    info = createInfo('stdlib::Al_sg225.ncmat;vdoslux=0')
    assert len(info.dyninfos)==1
    di = info.dyninfos[0]
    knl_orig = di.loadKernel(vdoslux=0)
    def get_knl():
        return copy.deepcopy(knl_orig)

    k = get_knl()
    _apply_smin( k, 1e-3 )
    plot_knl( k )
    _plot_sign( k )
    trimknl( k )
    _plot_sign( k )
    plot_knl( k )

    k = get_knl()
    _apply_smin( k, 1e-3 )
    _plot_sign( k )
    _set_upper_alpha_edge_to_zero( k, 5 )
    _plot_sign( k )
    trimknl( k )
    _plot_sign( k )

def main():
    testsimple()
    if False:
        #nb: uncomment FPE import above if plotting
        plot_testsab()

if __name__ == '__main__':
    main()
