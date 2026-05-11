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
#
import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.stabilise_ncpprint # noqa F401
from NCrystalDev._common import ncpprint
from NCrystalDev.misc import evaluate_query as ncquery
import numpy as np
import math

def plot( alphagrid, betagrid, cell_list, title ):
    import matplotlib.pyplot as plt
    fig,axis = plt.subplots()
    xi = betagrid
    yi = alphagrid

    #Draw grid:
    gridpts = []
    if 0.0 not in xi:
        axis.axvline(0.0, color='red', linestyle=':', linewidth=1, alpha=0.5)
    if 0.0 not in yi:
        axis.axhline(0.0, color='red', linestyle=':', linewidth=1, alpha=0.5)

    for x in xi:
        axis.axvline(x, color='lightgray', linewidth=0.5)
        gridpts += [ (x,y) for y in yi ]
    for y in yi:
        axis.axhline(y, color='lightgray', linewidth=0.5)
    axis.plot(*zip(*gridpts), 'o')

    #Collect energies curves:
    estr2count = {}
    eprev = None
    phasespace_draw_evals = []
    cell_numbering = []
    for e, ialpha, ibeta in cell_list:
        assert eprev is None or e >= eprev
        eprev, estr = e, str(e)
        if estr not in estr2count:
            nbr = len(estr2count)+1
            phasespace_draw_evals.append( ( e, nbr ) )
            estr2count[estr] = nbr
        cell_numbering.append( ( ialpha, ibeta, estr2count[estr], e ) )

    #Add phase-space curves:
    nbr2col = {}
    for e,nbr in phasespace_draw_evals:
        if math.isinf(e):
            nbr2col[nbr]='black'
            print("WARNING: Not showing phasespace curve for E/kT=infinity")
            continue
        bmin = max(-e,betagrid[0])
        bmax = betagrid[-1]
        if not ( bmax > bmin ):
            print(f"WARNING: Not showing phasespace curve for E/kT={e:g}")
        b = np.linspace(bmin,bmax,5000)
        if e > 0:
            sbe = np.sqrt(b+e)
            ap = ( sbe + np.sqrt(e) )**2
            am = ( sbe - np.sqrt(e) )**2
            color = plt.plot(b,ap)[0].get_color()
            plt.plot(b,am,color=color)
        else:
            color = plt.plot(b,b)[0].get_color()
        nbr2col[nbr]=color

    #Draw cell numbers:
    nbr_prev = None
    for ia, ib, nbr, e in cell_numbering:
        if nbr_prev is None or nbr_prev != nbr:
            print("Adding text for %i (E/kT=%g)"%(nbr,e))
            nbr_prev=nbr
        axis.text(0.5*(betagrid[ib-1]+betagrid[ib]),
                  0.5*(alphagrid[ia-1]+alphagrid[ia]),
                  str(nbr),ha='center',va='center',color=nbr2col[nbr])

    da = alphagrid[-1]-alphagrid[0]
    db = betagrid[-1]-betagrid[0]
    axis.set_ylim(-da*0.05,alphagrid[-1]+da*0.05)
    axis.set_xlim(betagrid[0]-db*0.05,betagrid[-1]+db*0.05)

    axis.set_xlabel('beta')
    axis.set_ylabel('alpha')
    axis.set_title(title)
    plt.show()

def surv(alphagrid, betagrid,do_plot=True):
    res = ncquery( [ 'sab','surveyor',
                     '@%s'%('@'.join(str(e) for e in alphagrid)),
                     '@%s'%('@'.join(str(e) for e in betagrid)) ] )
    ncpprint(res)
    if do_plot:
        plot( alphagrid, betagrid, res['touch_list'], 'cells touched' )
        plot( alphagrid, betagrid, res['cover_list'], 'cells covered' )

def main(do_plot):

    def s( alphagrid, betagrid ):
        surv(alphagrid, betagrid, do_plot = do_plot)

    s( [0.1, 0.9, 3.0], [-5.0,-3.0, 0.0, 2.0 ] )
    s( [0.0, 0.9, 3.0], [-5.0,-3.0, 0.0, 2.0 ] )
    s( [0.1,0.9,3.0], [-3.0,-0.5, 1.4, 4.0 ] )
    s( [0.0,1.0,3.0], [-3.0,-0.5, 1.4, 4.0 ] )
    s( np.linspace(0.01,100,5), np.linspace(-10.0,20.0,7) )
    s( [0.0e9, 0.9e9, 3.0e9], [-5.0e9,-3.0e9, 0.0, 2.0e9 ] )

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:])

