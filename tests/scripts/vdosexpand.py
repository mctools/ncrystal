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

from NCrystalDev.misc import evaluate_query as ncquery
from NCTestUtils.common import ( interp1d,
                                 interp1d_loglin,
                                 powspace,
                                 calc_reldiff,
                                 thicken_grid )
from NCrystalDev._numpy import _np_linspace
import shlex
import numpy as np

def standard_tests():
    do_test( 'Al_sg225.ncmat;vdoslux=0' )
    do_test( 'Ca_sg229_Calcium-gamma.ncmat;vdoslux=0' )
    do_test( 'Polyethylene_CH2.ncmat;vdoslux=0', 'H' )
    do_test( 'Al_sg225.ncmat;vdoslux=0', gnmax = 2 )
    #FIXME: Support and verify gnmax < 4!

_default_plot = [False]

def do_test( cfgstr, atomlbl = None, gnmax = None, do_plot = None):
    do_plot = _default_plot[0] if do_plot is None else do_plot
    if not do_plot:
        import NCTestUtils.enable_fpe # noqa F401


    print()
    print("====================================")
    print("===    Testing VDOS expansion    ===")
    print("====================================")
    print()
    print("Reproduce with command: %s"%build_cmdstr( cfgstr = cfgstr,# noqa U031
                                                     atomlbl=atomlbl,
                                                     gnmax=gnmax))
    query = ['vdos','expand',cfgstr]
    if atomlbl is not None or gnmax is not None:
        query.append(atomlbl or '')
    if gnmax is not None:
        query.append(str(gnmax))
    print(f"Query: ncrystal_query {shlex.join(query)}")
    data = ncquery(query,huge_arrays=True)
    assert set(data.keys()) == {'input','output'}
    data_in = data['input']
    data = data['output']
    #import pprint
    #pprint.pprint(data)
    kT = data['kT']

    #FIXME: Add some meaning tests here.
    gns = dict( (nm1+1,gn) for nm1,gn in enumerate(data['Gn']) )
    gnmax_actual = max(gns.keys())
    assert gnmax is None or gnmax_actual==gnmax

    #Prepare Gn(beta) functions:
    nvals = [ nm1+1 for nm1 in range(gnmax_actual) ]
    gnpts = []

    for n in nvals:
        gn_of_beta = gns[n]['values']*kT#*kT preserves unit integral
        emin,emax = gns[n]['energy_range']
        beta = _np_linspace( emin/kT, emax/kT, len(gn_of_beta ) )
        assert 0.999 < np.trapezoid(gn_of_beta, beta) < 1.001
        gnpts.append( (beta, gn_of_beta) )

    title = f'"{data_in["cfgstr"]}"'
    if atomlbl:
        title += f' ({atomlbl})'
    if gnmax:
        title += f' (Gn max forced to G{gnmax})'
    else:
        title += f' (up to G{gnmax_actual})'

    #Grid:
    alphagrid, betagrid = data['alphaGrid'], data['betaGrid']

    #print("ALPHA:",alphagrid)
    #bfine = np.linspace(betagrid[0],betagridbmin,bmax,10000)

    def plot_nselect(n):
        nfreq = 1 if n <=5 else (4 if n<=20 else (20 if n<=100 else 40) )
        return not ( nfreq > 1 and n!=gnmax_actual and n%nfreq!=0 )

    if do_plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, sharex=True)
        fig.subplots_adjust(hspace=0)
        ax, axdiff = axs
        ax.set_title(title)

        beta_highres = thicken_grid( betagrid, 20 )
#            pts_ref = fnalpha(alpha_highres,n)
#            pts_interp = interp1d_loglin(alphagrid,fnalpha(alphagrid,n))(alpha_highres)
#
        for n in nvals:
            if not plot_nselect(n):
                continue
            gn = gnpts[n-1]
            gnfct = interp1d( *gn )
            gnfct_knl = interp1d(betagrid,gnfct(betagrid))
            pts_ref = gnfct( beta_highres )
            pts_interp = gnfct_knl( beta_highres )
            col = ax.plot(*gn,label=f'G{n}')[0].get_color()
            ax.plot(beta_highres,pts_interp,color=col,ls=':')

            #ptsmax = pts_ref.max()
            #rdmask = pts_ref > ptsmax*(1e-6 if n==1 else 1e-3)

            ptsmax = pts_ref.max()
            rdmask = pts_ref > ptsmax*(1e-2 if n==1 else 1e-2)
            axdiff.plot(beta_highres[rdmask],
                        calc_reldiff(pts_ref[rdmask],pts_interp[rdmask]),
                        color=col,alpha=0.5)
            axdiff.axhline(y=1e-1, color="red", linestyle=":")
            axdiff.axhline(y=1e-2, color="green", linestyle=":")

        #intermediate curves as well:
        col_refE0ABGrid='black'
        col_combinedGnFct='grey'
        for pseudocurve, col, showpts in [('combinedGnFct',col_combinedGnFct,True),
                                          ('refE0ABGrid',col_refE0ABGrid,False)]:
            pc_fct_b, pc_fct_vals = data[pseudocurve]
            pc_fct_vals = pc_fct_vals*(gnpts[0][1].max()/pc_fct_vals.max())
            gnfct = interp1d( pc_fct_b, pc_fct_vals )
            gnfct_knl = interp1d(betagrid,gnfct(betagrid))
            pts_ref = gnfct( beta_highres )
            pts_interp = gnfct_knl( beta_highres )
            ax.plot(pc_fct_b, pc_fct_vals,
                    label=pseudocurve,lw=3, alpha=0.6, color=col)
            if showpts:
                ax.plot( betagrid, interp1d(pc_fct_b, pc_fct_vals)(betagrid),
                         'x', color='red' )#this is where we show the points.
            ax.plot(beta_highres, pts_interp, color=col, linestyle=':')
            axdiff.plot(beta_highres,
                        calc_reldiff(pts_ref,pts_interp),
                        color=col,alpha=0.5)

        for a in [ax,axdiff]:
            a.axvline(x=-data['betaMax'], color="red", linestyle=":",
                      label="\u00b1betaMax")
            a.axvline(x=data['betaMax'], color="red", linestyle=":")
            a.vlines( x=betagrid, ymin=0, ymax=1,
                      transform=a.get_xaxis_transform(),
                      color="lightgray", lw=1, alpha=0.7, zorder=0.5 )




        ax.set_title(title)
        axdiff.set_xlabel('beta')
        #plt.semilogy()
        #ax.grid()
        ax.legend(draggable=True,ncol=3)
        axdiff.semilogy()
        axdiff.set_ylim(1e-10,1e1)
        plt.show()

    #def alpha2
    #if do_plot:
    alpha2x = data['alpha2x']
    alphaMax = data['alphaMax']

    def log_fnalpha(alpha, n):
        x = np.asarray(alpha, dtype=float).copy()
        x *= alpha2x
        assert n==int(n) and n>=1
        lognfact = np.log(np.arange(1, n + 1)).sum()
        y = np.zeros_like(x)
        mask = x > 0
        y[mask] = -x[mask] + n * np.log(x[mask]) - lognfact
        y[x == 0] = float('-inf')
        return y.item() if y.ndim == 0 else y

    def fnalpha(alpha, n):
        return np.exp(log_fnalpha(alpha,n))

    if do_plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, sharex=True)
        fig.subplots_adjust(hspace=0)
        ax, axdiff = axs
        ax.set_title(title)

        #print("ALPHAGRID:",alphagrid)
        alpha_highres = thicken_grid( alphagrid, 20 )
#
#powspace( 1e-100, alphaMax*alpha2x,
#                                  len(alphagrid)*50, 2.0 )/alpha2x
        #alpha_highres[-1]=alphaMax
        #assert alpha_highres[-1]>alphaMax*1.0999
        #print("alphaMax=",alphaMax)
        #print("xmax=",alphaMax*alpha2x)

        ptsmax1 = None
        for n in nvals:
            if not plot_nselect(n):
                continue
            pts_ref = fnalpha(alpha_highres,n)
            pts_interp = interp1d_loglin(alphagrid,fnalpha(alphagrid,n))(alpha_highres)
            col = ax.plot(alpha_highres,pts_ref,label=f'n={n}')[0].get_color()
            #print("ARGH alpha_highres[-1]=",alpha_highres[-1],"alphagrid[-1]=",alphagrid[-1])
            ax.plot(alpha_highres,pts_interp,ls=':',color=col)
            ax.plot(alphagrid,fnalpha(alphagrid,n),'x',color=col)#FIXME ARE WE SHOWING THE PTS REPEATEDLY?????

            ptsmax = pts_ref.max()
            if n==1:
                ptsmax1=ptsmax
            rdmask = pts_ref > ptsmax*(1e-6 if n==1 else 1e-3)
            axdiff.plot(alpha_highres[rdmask],
                        calc_reldiff(pts_ref[rdmask],pts_interp[rdmask]),
                        color=col,alpha=0.5)
            axdiff.axhline(y=1e-1, color="red", linestyle=":")
            axdiff.axhline(y=1e-2, color="green", linestyle=":")
            #axdiff.axhline(y=1e-2, color="red", linestyle=":")

        #intermediate curves as well:
        for pseudocurve, col in [('refE0ABGrid',col_refE0ABGrid)]:
            pc_fct_a, pc_fct_vals = data[pseudocurve]
            pc_fct_vals = pc_fct_vals*(ptsmax1/pc_fct_vals.max())
            #rdmask = pc_fct_vals > ptsmax1*1e-6
            #pc_fct_a, pc_fct_vals = pc_fct_a[rdmask], pc_fct_vals[rdmask]
            gnfct = interp1d( pc_fct_a, pc_fct_vals )
            gnfct_knl = interp1d(alphagrid,gnfct(alphagrid))
            pts_ref = gnfct( alpha_highres )
            pts_interp = gnfct_knl( alpha_highres )
            ax.plot(pc_fct_a, pc_fct_vals,
                    label=pseudocurve,lw=3, alpha=0.6, color=col)
            ax.plot(alpha_highres, pts_interp, color=col, linestyle=':')
            axdiff.plot(alpha_highres,
                        calc_reldiff(pts_ref,pts_interp),
                        color=col,alpha=0.5)

        for a in [ax,axdiff]:
            a.axvline(x=alphaMax, color="red", linestyle="--",label="alphaMax")
            a.vlines( x=alphagrid, ymin=0, ymax=1,
                      transform=a.get_xaxis_transform(),
                      color="lightgray", lw=1, alpha=0.7, zorder=0.5 )

        axdiff.semilogy()
        axdiff.set_ylim(1e-10,1e1)
        axdiff.set_xlabel('alpha')
        #plt.semilogy()
        #plt.grid()
        ax.legend(draggable=True,ncol=3,fontsize='small')
        plt.show()

    print(f"nalpha = {len(alphagrid)}")
    print(f"nbeta = {len(betagrid)}")

        #data['betaMax']

#Support interactive usage:

def parse_args():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("cfgstr",help=('Material. Use cfgstr="standardtests" to'
                                  ' ignore other args and go through all'
                                  ' standard tests.'))
    p.add_argument("atomlbl",nargs="?")
    p.add_argument("--gnmax",type=int)
    p.add_argument("--noplot",action="store_true")
    return p.parse_args()

def build_cmdstr(cfgstr,atomlbl=None,gnmax=None):
    import os
    import sys
    prog = os.path.basename(sys.argv[0])
    if prog.startswith('sb_'):
        #simplebuild
        prog = prog.split('_',2)[-1]
    assert prog=='testvdosexpand'
    a=[ prog, cfgstr ]
    if atomlbl is not None:
        a += [atomlbl]
    if gnmax is not None:
        a += ["--gnmax",str(gnmax)]
    return shlex.join(a)

def main():
    import sys
    if len(sys.argv) != 1:
        args = parse_args()
        if args.cfgstr=='standardtests':
            _default_plot[0] = not args.noplot
            standard_tests()
        else:
            do_test( cfgstr = args.cfgstr,
                     atomlbl = args.atomlbl,
                     gnmax = args.gnmax,
                     do_plot = not args.noplot )
    else:
        print("No arguments provided so running standard set"
              " of tests without plotting.")
        standard_tests()

if __name__ == '__main__':
    main()
