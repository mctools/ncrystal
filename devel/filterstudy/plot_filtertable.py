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

"""Plot NCrystal filter tables ("filtertable" query) against the exact cross
sections, for one or more NCrystal cfg-strings.

The top plot shows the exact total cross sections (scattering + absorption)
and the piecewise linear tables. The bottom plot shows the relative difference
|table/exact - 1| on a log scale, together with the requested tolerance.

Examples (run with ./run.sh to use the NCrystal dev build):

  ./run.sh plot_filtertable.py stdlib::Al_sg225.ncmat
  ./run.sh plot_filtertable.py -e -x 1e-4 1 "stdlib::Be_sg194.ncmat;temp=80K" gasmix::air
  ./run.sh plot_filtertable.py --tol 1e-4 -x 1 10 stdlib::Fe_sg229_Iron-alpha.ncmat -o fe.png
"""
import argparse

import numpy as np

try:
    import NCrystalDev as NC
    from NCrystalDev.misc import evaluate_query
except ImportError:
    import NCrystal as NC
    from NCrystal.misc import evaluate_query

WL2EKIN = 0.081804209605330899  # E[eV] = WL2EKIN / wl[Aa]^2


def parse_args():
    p = argparse.ArgumentParser(
        description='Plot filter tables of NCrystal materials against the exact'
                    ' cross sections, and their relative differences.')
    p.add_argument('cfgstr', nargs='+', help='NCrystal cfg-string(s)')
    g = p.add_mutually_exclusive_group()
    g.add_argument('-w', '--wavelength', action='store_true',
                   help='Use wavelength [Aa] as x-axis (default)')
    g.add_argument('-e', '--energy', action='store_true',
                   help='Use neutron energy [eV] as x-axis')
    p.add_argument('-x', '--range', nargs=2, type=float, metavar=('MIN', 'MAX'),
                   help='Range of the x-axis, in the chosen unit (default: '
                        '0.01-100 Aa, or the corresponding energies)')
    p.add_argument('--tol', type=float, default=1e-3,
                   help='Relative tolerance of the tables (default: 1e-3)')
    p.add_argument('--npts', type=int, default=100000,
                   help='Number of points for the exact curves (default: 1e5)')
    p.add_argument('--linx', action='store_true', help='Linear x-axis')
    p.add_argument('--knots', action='store_true',
                   help='Mark the table points on the top plot')
    p.add_argument('-o', '--output', help='Save the plot to this file instead'
                   ' of showing it')
    return p.parse_args()


def bragg_edges(cfgstr):
    """Bragg edge wavelengths (2d) of all crystalline phases."""
    res = []

    def collect(info):
        if info.isMultiPhase():
            for _, ph in info.phases:
                collect(ph)
        elif info.hasHKLInfo():
            res.extend(2.0 * h.d for h in info.hklObjects())
    collect(NC.createInfo(cfgstr))
    return np.unique(res)


def exact_xs(cfgstr, wl):
    return (NC.createScatter(cfgstr).xsect(wl=wl)
            + NC.createAbsorption(cfgstr).xsect(wl=wl))


def table_eval(wl_t, xs_t, wl):
    """Linear interpolation in the table. At a Bragg edge (a pair of points
    with the same wavelength), the value above the edge is used for
    wl >= edge."""
    wl_t, xs_t = np.asarray(wl_t), np.asarray(xs_t)
    i = np.clip(np.searchsorted(wl_t, wl, side='right') - 1, 0, len(wl_t) - 2)
    x0, x1, y0, y1 = wl_t[i], wl_t[i + 1], xs_t[i], xs_t[i + 1]
    dx = np.where(x1 > x0, x1 - x0, 1.0)
    t = np.where(x1 > x0, (wl - x0) / dx, 0.0)
    return y0 + t * (y1 - y0)


def main():
    args = parse_args()
    use_energy = args.energy
    if args.range:
        xmin, xmax = sorted(args.range)
        if not xmin > 0:
            raise SystemExit('The range must be positive')
    else:
        xmin, xmax = ((WL2EKIN / 100.0**2, WL2EKIN / 0.01**2) if use_energy
                      else (0.01, 100.0))
    if use_energy:
        wlmin, wlmax = np.sqrt(WL2EKIN / xmax), np.sqrt(WL2EKIN / xmin)
    else:
        wlmin, wlmax = xmin, xmax

    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 8),
                                   gridspec_kw={'height_ratios': [2, 1]})
    to_x = (lambda wl: WL2EKIN / wl**2) if use_energy else (lambda wl: wl)
    relmin = 1e-12  # floor for the relative differences on the log-scale

    for icfg, cfgstr in enumerate(args.cfgstr):
        col = f'C{icfg}'
        # The table always starts at wavelength 0 (and ends at 500 Aa, unless
        # a longer range is plotted):
        query = ['filtertable', cfgstr, f'tol={args.tol:g}']
        if wlmax > 500.0:
            query.append(f'wlmax={wlmax:.17g}')
        r = evaluate_query(query)
        wl_t, xs_t = np.asarray(r['wl']), np.asarray(r['xs'])
        npts_in_range = int(np.sum((wl_t >= wlmin) & (wl_t <= wlmax)))

        # Exact curve: log-spaced points plus points close to each Bragg edge
        wl = np.geomspace(wlmin, wlmax, args.npts)
        e = bragg_edges(cfgstr)
        e = e[(e > wlmin) & (e < wlmax)]
        wl = np.unique(np.concatenate([wl, e * (1 - 1e-7), e * (1 + 1e-7)]))
        xs = exact_xs(cfgstr, wl)
        xs_tab = table_eval(wl_t, xs_t, wl)
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.where(xs > 0, np.abs(xs_tab / xs - 1.0),
                           np.where(xs_tab == 0, 0.0, np.inf))

        x = to_x(wl)
        label = f'{cfgstr} ({npts_in_range} table points in range)'
        ax1.plot(x, xs, color=col, lw=1.5, alpha=0.5, label=f'{cfgstr} (exact)')
        # The table as actually interpolated (linear in wavelength):
        ax1.plot(x, xs_tab, color=col, lw=0.8, ls='--',
                 label=f'{cfgstr} (table, {npts_in_range} points in range)')
        if args.knots:
            ax1.plot(to_x(wl_t), xs_t, color=col, ls='none', marker='.', ms=3)
        ax2.plot(x, np.maximum(rel, relmin), color=col, lw=0.8, label=label)
        print(f'{cfgstr}: {len(wl_t)} table points ({npts_in_range} in range),'
              f' {r["ndiscontinuities"]} discontinuities, max relative difference {rel.max():.3g}'
              f' (tol={args.tol:g})')

    ax2.axhline(args.tol, color='k', ls=':', lw=1, label=f'tol={args.tol:g}')
    xlabel = 'Neutron energy [eV]' if use_energy else r'Wavelength [$\AA$]'
    ax2.set_xlabel(xlabel)
    ax1.set_ylabel('Cross section [barn/atom]')
    ax2.set_ylabel('|table/exact - 1|')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(relmin, None)
    if not args.linx:
        ax1.set_xscale('log')
    ax1.set_xlim(xmin, xmax)
    ax1.set_title(f'NCrystal filter tables (tol={args.tol:g})')
    ax1.legend(fontsize='small')
    ax2.legend(fontsize='small', loc='lower right')
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)
    fig.tight_layout()
    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f'Saved {args.output}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
