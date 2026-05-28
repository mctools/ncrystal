
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

class RefCell:

    """Reference evaluation of integrals within a particular S(alpha,beta) grid
    cell, including the usual loglin/linlin interpolation scheme.

    """

    def __init__(self, *, a1,a2,b1,b2,s11,s12,s21,s22,
                 mpmath_dps = 50 ):
        from .mpmathctx import get_mpmath_context
        self.__mp = get_mpmath_context(mpmath_dps)
        mpf = self.__mp.mpf
        self.__alpha = (mpf(a1),mpf(a2))
        self.__beta = (mpf(b1),mpf(b2))
        self.__s = (mpf(s11),mpf(s12),mpf(s21),mpf(s22))
        a1,a2 = self.__alpha
        b1,b2 = self.__beta
        assert 0 <= a1 < a2
        assert b1 < b2
        assert all(s>=0 for s in self.__s)
        self.__fullint = None
        self.__psint = {}
        self.__psintbf = {}
        self.__psranges = {}

    def full_integral( self ):
        if self.__fullint is None:
            a1,a2 = self.__alpha
            b1,b2 = self.__beta
            s11,s12,s21,s22 = self.__s
            integA_at_b1 = _alphaintegral1d(a1,a2,s11,s12,a1,a2,self.__mp)
            integA_at_b2 = _alphaintegral1d(a1,a2,s21,s22,a1,a2,self.__mp)
            self.__fullint = (b2-b1)*(integA_at_b1+integA_at_b2)/2
        return self.__fullint

    def phasespace_integral_ranges( self, E_div_kT ):
        e = self.__mp.mpf(E_div_kT)
        cv = self.__psranges.get(e)
        if cv is not None:
            return cv
        a1,a2 = self.__alpha
        b1,b2 = self.__beta
        v = _find_integration_ranges( e, (a1,a2),(b1,b2) )
        self.__psranges[e] = v
        return v

    def phasespace_integral_brute_force( self, E_div_kT, n=1000 ):
        e = self.__mp.mpf(E_div_kT)
        key=(e,n)
        cv = self.__psintbf.get(key)
        if cv is not None:
            return cv
        v= _brute_force_integral( E_div_kT=e,
                                  alpha=self.__alpha,
                                  beta=self.__beta,
                                  svals=self.__s,
                                  n=n )
        self.__psintbf[key] = v
        return v

    def phasespace_integral( self, E_div_kT ):
        e = self.__mp.mpf(E_div_kT)
        cv = self.__psint.get(e)
        if cv is not None:
            return cv
        ranges = self.phasespace_integral_ranges(e)
        data = dict( alpha = self.__alpha,
                     beta = self.__beta,
                     S = self.__s )
        v = _ref_integral( data, ranges, self.__mp )
        self.__psint[e] = v
        return v

def _mp_integral_k_pow_x(k, a, b,mp):
    #Integrate k^x over x in [a,b].
    assert b>=a
    assert k>0
    mpf = mp.mpf
    k,a,b = mpf(k), mpf(a), mpf(b)
    if k == 1 or b==a:
        return b - a
    lnk = mp.log1p(k - 1) if abs(k-1)<0.1 else mp.log(k)
    eps = (b - a) * lnk
    if abs(eps)>0.01 and abs(lnk)>1e-8:
        #Exact formula (should we use this more often?)
        return mp.exp(a * lnk) * mp.expm1(eps) / lnk

    # Taylor expand exact formula, using lnk=eps/(b-a). Exact formula can be
    # rewritten as (b-a)*exp(a*lnk)*[ (exp(eps)-1)/eps ]
    # The usual taylor expansion of exp(eps) gives:
    #
    # (exp(eps)-1)/eps = sum_{j=1}^{inf}( eps^(j-1)/j! )
    # The j=1 term is 1 and for the other terms term(j)/term(j-1) = eps/j
    #
    assert abs(eps)<0.5
    targetprec = mpf(10)**(-mp.mp.dps)
    term = mpf(1) # j=1
    totsum = term
    nlim = 2000
    for j in range(2,nlim+1):
        term *= (eps/j)
        totsum += term
        if j > 5 and abs(term) < abs(totsum)*targetprec:
            break
    assert j<nlim, "slow convergence"
    return mp.exp(a * lnk) * (b - a) * totsum

def _alphaintegral1d(a1,a2,s1,s2,grid_a1, grid_a2, mp):
    #For S(a) given as a linear interpolation in log(S) (or linear in S if s1 or
    #s2 is zero) between (a1,s1) and (a2,s2), integrate S(a) over the subrange
    #[grid_a1,grid_a2].
    mpf = mp.mpf
    a1,a2,s1,s2 = mpf(a1),mpf(a2),mpf(s1),mpf(s2)
    grid_a1, grid_a2 = mpf(grid_a1), mpf(grid_a2)
    if s1==0 or s2==0:
        #linear interpolation in s [fallback mode]
        # Result is the s-value at the midpoint of [grid_a1,grid_a2] times the
        # width if this interval.
        if grid_a1==a1 and grid_a2==a2:
            return ( (a2-a1) * (s1 + s2 ) ) / 2
        else:
            dsda,dg,avg=(s2-s1)/(a2-a1),grid_a2-grid_a1,(grid_a2+grid_a1)/2
            return dg * (  s1 + (avg-a1)*dsda )
    else:
        #linear interpolation in log(s) [primary mode]
        k, da = s2/s1, a2-a1
        if grid_a1==a1 and grid_a2==a2:
            r1, r2 = mpf(0), mpf(1)
        else:
            invda = mpf(1)/da
            r1 = (grid_a1 - a1)*invda
            r2 = (grid_a2 - a1)*invda if grid_a2!=a2 else mpf(1)
        return s1 * da * _mp_integral_k_pow_x(k, r1, r2, mp)

def _create_s_of_a(a1,a2,s1,s2,mp):
    #Return S(a) which is given as a linear interpolation in log(S) (or linear
    #in S if s1 or s2 is zero) between (a1,s1) and (a2,s2).
    mpf = mp.mpf
    assert 0<=a1<=a2<=1e308
    assert s1>=0 and s2>=0 and s1<=1e308 and s2<=1e308
    s1,s2,a1,a2 = mpf(s1),mpf(s2),mpf(a1),mpf(a2)
    invda = 1/(a2-a1)
    if s1==0 or s2==0:
        #linear interpolation in s [fallback mode]
        def s_of_a_interp(r):
            return (1-r)*s1+r*s2
    else:
        #linear interpolation in log(s) [primary mode]
        k = s2/s1
        def s_of_a_interp(r):
            return s1*(k**r)
    def s_of_a(a):
        if a==a1:
            return s1
        elif a==a2:
            return s2
        else:
            r = (a-a1)*invda
            assert -1e-15 <= r <= 1+1e-15
            r = max(mpf(0),min(mpf(1),r))
            return s_of_a_interp(r)
    return s_of_a

def _ref_integral( data, cellranges, mp ):
    mpf = mp.mpf
    e = cellranges['E_div_kT']
    (a1, a2), (b1, b2) = data['alpha'], data['beta']
    e, a1, a2, b1, b2 = mpf(e), mpf(a1), mpf(a2), mpf(b1), mpf(b2)
    s11, s12, s21, s22 = (mpf(e) for e in data['S'])
    totsum = mpf(0)
    for ( r_a1, r_a2,
          (clip_betaminus, clip_betaplus)) in cellranges['ranges']:
        r_a1, r_a2 = mpf(r_a1), mpf(r_a2)
        if not (clip_betaminus or clip_betaplus):
            #Full box integral of region => no need for numerical quadrature:
            aint_at_b1 = _alphaintegral1d(a1,a2,s11,s12,r_a1, r_a2,mp)
            aint_at_b2 = _alphaintegral1d(a1,a2,s21,s22,r_a1, r_a2,mp)
            contrib = (aint_at_b1+aint_at_b2)*(b2-b1)/2
            totsum += contrib
            continue
        if not clip_betaminus:
            def b_low(a):
                return b1
        else:
            def b_low(a):
                return a - 2*mp.sqrt(e*a)
        if not clip_betaplus:
            def b_up(a):
                return b2
        else:
            def b_up(a):
                return a + 2*mp.sqrt(e*a)
        sofa_at_b1 = _create_s_of_a(a1,a2,s11,s12,mp)
        sofa_at_b2 = _create_s_of_a(a1,a2,s21,s22,mp)
        invb2mb1 = 1/(b2-b1)
        def contrib_at_a(a):
            #always interpolate linearly in b:
            s1, s2 = sofa_at_b1(a), sofa_at_b2(a)
            bl, bu = b_low(a), b_up(a)
            bmiddle = (bu+bl)/2
            smiddle = s1 + (s2-s1)*(bmiddle-b1)*invb2mb1
            return (bu-bl)*smiddle
        contrib = mp.quad(contrib_at_a,[r_a1,r_a2],epsrel=1e-20,epsabs=1e-20)
        totsum += contrib
    return totsum

def _find_integration_ranges( E_div_kT, alpha, beta ):
    e, (a1, a2), (b1, b2) = E_div_kT, alpha, beta
    if b2 <= -e:
        return []
    _ = 2*((e*(b2+e))**0.5)
    ap2 = 2*e+b2+_ #alpha^+(b2)
    if a1 >= ap2:
        return []
    am2 = 2*e+b2-_ #alpha^-(b2)
    if b1 > -e:
        _ = 2*((e*(b1+e))**0.5)
        am1 = 2*e+b1-_ #alpha^-(b1)
        ap1 = 2*e+b1+_#alpha^+(b1)
    else:
        am1,ap1=None,None
    a2 = min(a2,ap2)
    if b2 < 0:
        a1 = max(a1,am2)
    elif b1 > 0:
        a1 = max(a1,am1)
    if not a2 > a1:
        return []
    #Look at lower bounds:
    intervals_lower = []
    if b1 <= -e:
        intervals_lower.append( (a1,a2,True) )#bounded by beta^-(alpha)
    else:
        _a1 = a1
        assert am1 is not None
        assert ap1 is not None
        if _a1 < am1:
            intervals_lower.append( (_a1,am1,True) )#bounded by beta^-(alpha)
            _a1 = am1
        if _a1 < ap1:
            intervals_lower.append( (_a1,min(ap1,a2),False) )#bounded by b1 edge
            _a1 = min(ap1,a2)
        if _a1 < a2:
            intervals_lower.append( (_a1,a2,True) )#bounded by beta^-(alpha)
    #Look at upper bounds:
    intervals_upper = []
    if b2 <= 0.0:
        am2 = 0.0
    _a1 = a1
    if _a1 < am2:
        intervals_upper.append( (_a1,min(am2,a2),True) )#bounded by beta^+(alph)
        _a1 = min(am2,a2)
    if _a1 < a2:
        intervals_upper.append( (_a1,a2,False) )#bounded by b2 edge
    #Combine intervals:
    assert intervals_lower and intervals_upper
    assert intervals_lower[0][0] == intervals_upper[0][0]
    assert intervals_lower[-1][1] == intervals_upper[-1][1]
    res = []
    a1 = intervals_upper[0][0]
    a2 = intervals_upper[-1][1]
    au = a2
    while intervals_lower:
        flags = ( intervals_lower[-1][2], intervals_upper[-1][2])
        if intervals_lower[-1][0] >= intervals_upper[-1][0]:
            al = intervals_lower[-1][0]
            intervals_lower.pop()
        else:
            al = intervals_upper[-1][0]
            intervals_upper.pop()
        res.append( (al,au,flags) )
        au = al
    assert len(intervals_upper)==1
    return dict( ranges = res, E_div_kT = E_div_kT )

def _brute_force_integral_impl( E_div_kT, alpha, beta, svals, n ):
    from NCrystalDev._numpy import _np_linspace
    import numpy as np
    e, (a1, a2), (b1, b2) = E_div_kT, alpha, beta
    assert not np.isinf(e)
    if b2 <= -e:
        return 0.0
    _ = 2*((e*(b2+e))**0.5)
    ap2 = 2*e+b2+_ #alpha^+(b2)
    if a1 >= ap2:
        return 0.0
    da = (a2-a1)/n
    db = (b2-b1)/n
    a = _np_linspace(a1+0.5*da,a2-0.5*da,n)
    b = _np_linspace(b1+0.5*db,b2-0.5*db,n)
    ab = np.column_stack((np.repeat(a, b.size), np.tile(b, a.size)))
    assert len(ab) == n*n
    #assert all(s>0.0 for s in svals), "lin fallback not implemented yet"#fixme
    s11, s12, s21, s22 = svals
    aa, bb = ab.T[0], ab.T[1]
    assert len(aa) == n*n
    assert len(bb) == n*n
    ra = (aa-a1)*(1.0/(a2-a1))
    rb = (bb-b1)*(1.0/(b2-b1))

    if s12==0.0 or s11==0.0:
        sofa_at_b1 = s11*(1.0-ra)+s12*ra
    else:
        sofa_at_b1 = ((s12/s11)**ra)*s11

    if s22==0.0 or s21==0.0:
        sofa_at_b2 = s21*(1.0-ra)+s22*ra
    else:
        sofa_at_b2 = ((s22/s21)**ra)*s21

    ss = sofa_at_b1*(1.0-rb)+sofa_at_b2*rb
    assert len(ss) == n*n
    phasespace_mask = (4*e*aa >= (bb-aa)**2)
    return ss.sum()*da*db, ss[phasespace_mask].sum()*da*db

def _brute_force_integral( E_div_kT, alpha, beta, svals, n ):
    e = float(E_div_kT)
    bfargs = dict( E_div_kT = e,
                   alpha = ( float(alpha[0]), float(alpha[1]) ),
                   beta = ( float(beta[0]), float(beta[1]) ),
                   svals = tuple( float(e) for e in svals ) )
    full, pb = _brute_force_integral_impl(**bfargs, n = n)
    return { 'full_integral': full,
             'phasespace_integral': pb,
             'phasespace_E_div_kT' : e }
