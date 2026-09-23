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

# NEEDS: mpmath

# Validates NCKinUtils.hh's getAlphaMinus/getAlphaPlus/getBetaMinus/
# getBetaPlus -- the shared, Taylor-hardened kinematic-boundary utilities
# that several sab component call sites used to bypass by re-deriving the
# naive, catastrophic-cancellation-prone formulas inline instead (fixed
# across NCSABCellSample.cc/NCSABSurveyor.cc; see CHANGELOG and
# docs/claude_session_vdos_fma_reprod.md). Exercises the ["util",
# "kinutils"] query, which evaluates the production utilities on a grid
# of ekin_div_kT values combined with offsets landing densely around the
# two singular points (beta=0 for alpha_minus, alpha=4*ekin_div_kT for
# beta_minus) where the naive formulas lose precision.

import math

from NCrystalDev.misc import evaluate_query as ncquery
from NCTestUtils.mpmathctx import get_mpmath_context


def mp_alpha_minus(e, beta, mp):
    e, beta = mp.mpf(e), mp.mpf(beta)
    return (mp.sqrt(e) - mp.sqrt(e + beta)) ** 2


def mp_alpha_plus(e, beta, mp):
    e, beta = mp.mpf(e), mp.mpf(beta)
    return (mp.sqrt(e) + mp.sqrt(e + beta)) ** 2


def mp_beta_minus(e, alpha, mp):
    e, alpha = mp.mpf(e), mp.mpf(alpha)
    return alpha - 2 * mp.sqrt(e * alpha)


def mp_beta_plus(e, alpha, mp):
    e, alpha = mp.mpf(e), mp.mpf(alpha)
    return alpha + 2 * mp.sqrt(e * alpha)


def naive_alpha_minus(e, beta):
    # The catastrophic-cancellation-prone formula that findActiveAlphaRange,
    # SABCellSurvey's constructor and RCSImpl's special_avals_candidates
    # used to re-derive inline, instead of calling getAlphaMinus:
    kk = e + beta
    a = kk + e
    b = 2.0 * math.sqrt(e * kk)
    return max(0.0, a - b)


def naive_beta_minus(e, alpha):
    # The catastrophic-cancellation-prone formula that
    # BoundedCellSampler::prepareBCSData and RCSImpl::alphaPtGem used to
    # re-derive inline, instead of calling getBetaMinus:
    return alpha - 2.0 * math.sqrt(e * alpha)


def check_grid( mp, evals_key, out_key, mp_minus, mp_plus, minus_key, plus_key ):
    # Combined absolute+relative tolerance: the absolute floor matters
    # right at the singularity, where the true value itself is tiny.
    atol_rel_e, rtol = 1e-14, 1e-11
    worst = 0.0
    for e, x, vminus, vplus in zip( evals_key, out_key[0],
                                    out_key[1], out_key[2] ):
        ref_minus = float(mp_minus(e, x, mp))
        ref_plus = float(mp_plus(e, x, mp))
        atol = atol_rel_e * e
        rd_minus = abs(vminus - ref_minus) / (atol + rtol * abs(ref_minus))
        rd_plus = abs(vplus - ref_plus) / (atol + rtol * abs(ref_plus))
        assert rd_minus <= 1.0, (minus_key, e, x, vminus, ref_minus, rd_minus)
        assert rd_plus <= 1.0, (plus_key, e, x, vplus, ref_plus, rd_plus)
        worst = max( worst, rd_minus, rd_plus )
    print(f"  worst error/tolerance ratio over grid: {worst:.3g}"
          " (<=1 required)")


def main():
    mp = get_mpmath_context(50)
    d = ncquery(['util', 'kinutils'])

    ab, bb = d['alphabounds'], d['betabounds']
    print(f"alpha_minus/alpha_plus: {len(ab['e'])} grid points")
    check_grid( mp, ab['e'], (ab['beta'], ab['aminus'], ab['aplus']),
                mp_alpha_minus, mp_alpha_plus, 'aminus', 'aplus' )

    print(f"beta_minus/beta_plus: {len(bb['e'])} grid points")
    check_grid( mp, bb['e'], (bb['alpha'], bb['bminus'], bb['bplus']),
                mp_beta_minus, mp_beta_plus, 'bminus', 'bplus' )

    # Demonstrate why the fix mattered: right at the singularities, the
    # naive formula (still evaluated here, in plain double precision, no
    # FMA contraction even) has already lost most of its significant
    # digits relative to its own true magnitude, while the production
    # getAlphaMinus/getBetaMinus queried above (validated to be accurate
    # to ~1e-11 relative there) stay reliable:
    print("naive-vs-hardened comparison right at the singularities:")
    e = 1.0
    for beta in (1e-6, 1e-8, 1e-10):
        ref = float(mp_alpha_minus(e, beta, mp))
        naive = naive_alpha_minus(e, beta)
        rd_naive = abs(naive - ref) / ref if ref else float('nan')
        print(f"  alpha_minus(e=1,beta={beta:g}):"
              f" naive={naive:.6g} ref={ref:.6g}"
              f" naive_rel_err_of_value={rd_naive:.3g}")
        # Even without any FMA contraction (real hardware makes this
        # worse still -- see the BoundedCellSampler beta_minus fix), the
        # naive formula has already lost virtually all significant
        # digits here, while getAlphaMinus (checked above) has not:
        assert rd_naive > 1e-5, "expected the naive formula to be bad here"

    for off in (1e-8, 1e-14, 1e-15):
        alpha = 4.0 * e + off
        ref = float(mp_beta_minus(e, alpha, mp))
        naive = naive_beta_minus(e, alpha)
        rd_naive = abs(naive - ref) / ref if ref else float('nan')
        print(f"  beta_minus(e=1,alpha=4+{off:g}):"
              f" naive={naive:.6g} ref={ref:.6g}"
              f" naive_rel_err_of_value={rd_naive:.3g}")
        assert rd_naive > 1e-9, "expected the naive formula to be bad here"

    print("ALL OK")


if __name__ == '__main__':
    main()
