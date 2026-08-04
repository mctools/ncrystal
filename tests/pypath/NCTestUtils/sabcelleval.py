
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

    Also able to sample (alpha,beta) points within the cell and phasespace.

    """

    def __init__(self, *, a1,a2,b1,b2,s11,s12,s21,s22,
                 mpmath_dps = 100 ):
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
        v = _ref_integral( data, ranges, self.__mp, do_sample = False )
        self.__psint[e] = v
        return v

    def sample( self, E_div_kT, n, seed = 123456789123456789 ):
        #Change the global random seed is the only way to reseed
        #mpmath.mp.rand().
        ranges = self.phasespace_integral_ranges(E_div_kT)
        if len(ranges)==0:
            #no contribution!
            return []
        data = dict( alpha = self.__alpha,
                     beta = self.__beta,
                     S = self.__s )
        #fixme: cache this _ref_integral(..) result?
        ri = _ref_integral( data, ranges, self.__mp, do_sample = True )
        from .common import change_random_seed
        with change_random_seed( seed ):
            one = self.__mp.mpf(1)
            def rngfct():
                return one-self.__mp.rand()#we want (0,1] not [0,1)
            return _sample_impl( E_div_kT = E_div_kT, n = n, rng = rngfct,
                                 mp = self.__mp, refintdata = ri )

def draw_alpha_beta_grid(alphagrid,betagrid,**kw_plot):
    from NCrystalDev.plot import PlotContext
    pctx = PlotContext(**kw_plot).check_unused()
    xi, yi, axis = betagrid, alphagrid, pctx.axis
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
    axis.set_xlabel('beta')
    axis.set_ylabel('alpha')
    return pctx.finalise( do_grid = False )

def plot_celleval( data, do_title=True, **kw_plot ):
    import numpy as np
    from NCrystalDev._numpy import _np_linspace
    from NCrystalDev.plot import PlotContext

    pctx = PlotContext(**kw_plot).check_unused()
    draw_alpha_beta_grid( alphagrid = data['alpha'],
                          betagrid = data['beta'],
                          **pctx.kwargs_subcontext() )
    b = data['beta']
    db = b[1]-b[0]
    a = data['alpha']
    da = a[1]-a[0]
    blim = ( b[0]-0.1*db, b[1]+0.1*db )
    alim = ( max(0.0,a[0]-0.1*da), a[1]+0.1*da )
    pctx.axis.set_xlim( *blim )
    pctx.axis.set_ylim( *alim )

    data_ci = data['cellintegral']
    elist = [ ( data['surveyor']['E_div_kT_touch'], 'touch', 'green', 0.0 ),
              ( data['surveyor']['E_div_kT_cover'], 'cover', 'blue',
                data_ci['full_integral'] ),
              ( data_ci['phasespace_E_div_kT'], 'chosen', 'red',
                list( v for k,v in data_ci['phasespace_integral']
                      if k=='Romberg33' )[0] ) ]
    for e, lbl, col, integral in elist:
        assert blim[1]>=-e#revisit this if it fails
        brangeplot = [max(blim[0],-e),blim[1]]
        assert brangeplot[1] > brangeplot[0]
        lble = f'{e:g}kT' if not np.isinf(e) else 'INF'
        lbl = f'{lbl} ({lble}, integral={integral:g})'
        if e > 0 and not np.isinf(e):
            b = _np_linspace(*brangeplot,5000)
            sbe = np.sqrt(b+e)
            ap = ( sbe + np.sqrt(e) )**2
            am = ( sbe - np.sqrt(e) )**2
            pctx.axis.plot(b,ap,color=col,label=lbl)
            pctx.axis.plot(b,am,color=col)
        elif np.isinf(e):
            pctx.axis.plot(brangeplot,[0,0],color=col,label=lbl)
        else:
            pctx.axis.plot(brangeplot,brangeplot,color=col,label=lbl)

    e = data_ci['phasespace_E_div_kT']
    b1, b2 = data['beta']
    a1, a2 = data['alpha']
    from matplotlib import patches
    for ( r_a1, r_a2,clip_betaminus,
          clip_betaplus ) in data_ci['integration_regions']['regions']:
        print(f"AlphaRange [{r_a1},{r_a2}]: clip_betaminus"
              f"={clip_betaminus}, clip_betaplus={clip_betaplus}")

        assert a1 <= r_a1 <= a2
        color=None
        if not ( clip_betaminus or clip_betaplus ):
            #Just a square!
            r = patches.Rectangle((b1, r_a1), b2-b1, r_a2-r_a1,
                                  facecolor=color, edgecolor='k',#'lightblue'
                                  hatch='///', linewidth=1.0)
            pctx.axis.add_patch(r)
            continue
        for aval in _np_linspace(r_a1,r_a2,50):
            #phasespace curve: 4ae=(b-a)^2 <=> |b-a|=sqrt(4ae)
            db = np.sqrt(4*aval*e)
            bm = aval-db if clip_betaminus else b1
            bp = aval+db if clip_betaplus else b2
            _=pctx.axis.plot([bm,bp],[aval,aval],color=color,alpha=0.3)
            if color is None:
                color=_[0].get_color()

    if do_title:
        title='s11=%g, s12=%g, s21=%g, s22=%g'%tuple(data['S'])
        pctx.axis.set_title(title)
    return pctx.finalise( do_grid = False, do_legend='draggable' )

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
    targetprec = mpf(10)**(-mp.dps)
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
    assert a2>a1
    grid_a1, grid_a2 = mpf(grid_a1), mpf(grid_a2)
    assert grid_a2>=grid_a1
    if grid_a2==grid_a1:
        return mp.mpf(0)
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

class ContribAtA:
    def __init__(self, mp, e, a1,a2,b1,b2,
                 r_a1, r_a2, clip_betaminus, clip_betaplus,
                 sofa_at_b1, sofa_at_b2):
        if not clip_betaminus:
            self.__b_low = lambda a : b1
        else:
            self.__b_low = lambda a : a - 2*mp.sqrt(e*a)
        if not clip_betaplus:
            self.__b_up = lambda a : b2
        else:
            self.__b_up = lambda a : a + 2*mp.sqrt(e*a)
        self.__invb2mb1 = 1/(b2-b1)
        self.__b1 = b1
        self.__r_a1 = r_a1
        self.__r_a2 = r_a2
        self.__sofa_at_b1 = sofa_at_b1
        self.__sofa_at_b2 = sofa_at_b2
    def __call__(self, a ):
        assert a>=self.__r_a1
        assert a<=self.__r_a2
        s1 = self.__sofa_at_b1(a)
        s2 = self.__sofa_at_b2(a)
        bl, bu = self.__b_low(a), self.__b_up(a)
        bmiddle = (bu+bl)/2
        smiddle = s1 + (s2-s1)*(bmiddle-self.__b1)*self.__invb2mb1
        return (bu-bl)*smiddle

def _ref_integral( data, cellranges, mp, do_sample ):
    mpf = mp.mpf
    e = cellranges['E_div_kT']
    (a1, a2), (b1, b2) = data['alpha'], data['beta']
    e, a1, a2, b1, b2 = mpf(e), mpf(a1), mpf(a2), mpf(b1), mpf(b2)
    s11, s12, s21, s22 = (mpf(e) for e in data['S'])
    sofa_at_b1 = _create_s_of_a(a1,a2,s11,s12,mp)
    sofa_at_b2 = _create_s_of_a(a1,a2,s21,s22,mp)

    totsum = mpf(0)
    if do_sample:
        sampleinfo = dict( e=e, a1=a1, a2=a2, b1=b1, b2=b2,
                           s11=s11, s12=s12, s21=s21, s22=s22)
        sampleregions = []
        sampleinfo['regions'] = sampleregions

    for ( r_a1, r_a2, (clip_betaminus, clip_betaplus)) in cellranges['ranges']:
        r_a1, r_a2 = mpf(r_a1), mpf(r_a2)
        assert r_a2 >= r_a1
        if do_sample or not (clip_betaminus or clip_betaplus):
            aint_at_b1 = _alphaintegral1d(a1,a2,s11,s12,r_a1,r_a2,mp)
            aint_at_b2 = _alphaintegral1d(a1,a2,s21,s22,r_a1,r_a2,mp)
        if do_sample:
            rangeinfo = ( r_a1, r_a2, (clip_betaminus,clip_betaplus) )
            sampleregions.append( dict( rangeinfo = rangeinfo,
                                        aint_at_b1 = aint_at_b1,
                                        aint_at_b2 = aint_at_b2 ) )

        if not (clip_betaminus or clip_betaplus):
            #Full box integral of region => no need for numerical quadrature:
            contrib = (aint_at_b1+aint_at_b2)*(b2-b1)/2
            totsum += contrib
            if do_sample:
                sampleregions[-1]['contrib']=contrib
            continue

        contrib_at_a = ContribAtA(mp=mp,e=e, a1=a1,a2=a2,b1=b1,b2=b2,
                                  r_a1=r_a1, r_a2=r_a2,
                                  clip_betaminus=clip_betaminus,
                                  clip_betaplus=clip_betaplus,
                                  sofa_at_b1=sofa_at_b1,
                                  sofa_at_b2=sofa_at_b2)

        contrib, err = mp.quad(contrib_at_a,[r_a1,r_a2],
                               epsrel=1e-20,epsabs=1e-20,error=True)
        assert abs(err)<=abs(contrib*1e-20)
        totsum += contrib
        if do_sample:
            _inv_db = (b2-b1)**(-1)
            sampleregions[-1]['contrib']=contrib
            sampleregions[-1]['contrib_at_a']=contrib_at_a
    if do_sample:
        sampleinfo['contrib_total'] = totsum
        if totsum>=0.0:
            cc = mpf(0.0)
            for i in range(len(sampleinfo['regions'])):
                cc += sampleinfo['regions'][i]['contrib']
                sampleinfo['regions'][i]['R_select_cumul'] = cc/totsum
            assert mp.almosteq( sampleinfo['regions'][-1]['R_select_cumul'],
                                mpf(1) )
            return sampleinfo
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
        #bounded by beta^-(alpha)
        if a1<e<a2:
            #potentially better quadrature if breaking at phasespace endpoint
            intervals_lower.append( (a1,e,True) )
            intervals_lower.append( (e,a2,True) )
        else:
            intervals_lower.append( (a1,a2,True) )
    else:
        _a1 = a1
        assert am1 is not None
        assert ap1 is not None
        if _a1 < min(am1,a2):
            intervals_lower.append( (_a1,min(am1,a2),True) )#bounded by b^-(a)
            _a1 = min(am1,a2)
        if _a1 < min(ap1,a2):
            intervals_lower.append( (_a1,min(ap1,a2),False) )#bounded by b1 edge
            _a1 = min(ap1,a2)
        if _a1 < a2:
            intervals_lower.append( (_a1,a2,True) )#bounded by beta^-(alpha)

    assert intervals_lower[-1][1]==a2
    #Look at upper bounds:
    intervals_upper = []
    if b2 <= 0.0:
        am2 = 0.0
    _a1 = a1
    if _a1 < min(am2,a2):
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
    assert a2 > a1
    au = a2
    while intervals_lower:
        flags = ( intervals_lower[-1][2], intervals_upper[-1][2])
        if intervals_lower[-1][0] >= intervals_upper[-1][0]:
            al = intervals_lower[-1][0]
            intervals_lower.pop()
        else:
            al = intervals_upper[-1][0]
            intervals_upper.pop()
        assert au > al
        res.append( (al,au,flags) )
        au = al
    assert len(intervals_upper)==1
    return dict( ranges = res, E_div_kT = E_div_kT )

def _brute_force_integral_impl( E_div_kT, alpha, beta, svals, n ):
    import numpy as np
    from NCrystalDev._numpy import _np_linspace
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

def _sample_impl( E_div_kT, n, rng, mp, refintdata ):
    e = mp.mpf(E_div_kT)
    regions = refintdata['regions']
    pthresholds = [r['R_select_cumul'] for r in regions]

    #Determine number of samples in each region:
    if len(pthresholds)==1:
        n_per_region = [n]
    else:
        n_per_region=[]
        #mpmath does not have multinomial distribution directly.
        p = [rng() for i in range(n)]
        p.sort()
        nnprev=0
        for pthr in pthresholds:
            nn = next((i for (i,pval) in enumerate(p) if pval>pthr),
                      len(p))
            n_per_region.append( nn-nnprev )
            nnprev = nn
    navg = [float(e*n) for e in pthresholds]
    for i in range(len(navg)-1,0,-1):
        navg[i] -= navg[i-1]

    a1 = refintdata['a1']
    a2 = refintdata['a2']
    b1 = refintdata['b1']
    b2 = refintdata['b2']
    s11 = refintdata['s11']
    s12 = refintdata['s12']
    s21 = refintdata['s21']
    s22 = refintdata['s22']
    lin_s_at_b1 = not min(s11,s12)>0
    lin_s_at_b2 = not min(s21,s22)>0

    sofa_at_b1 = _create_s_of_a(a1,a2,s11,s12,mp)
    sofa_at_b2 = _create_s_of_a(a1,a2,s21,s22,mp)

    db = b2 - b1
    invdb = 1/db
    e = refintdata['e']
    result = []
    for nreg, reg in zip(n_per_region,regions):
        if nreg==0:
            continue
        assert reg['contrib'] > 0
        r_a1, r_a2, (clip_betaminus,clip_betaplus) = reg['rangeinfo']
        W1 = reg['aint_at_b1']
        W2 = reg['aint_at_b2']
        fullboxcontrib = (W1+W2)*db/2
        fullboxAR = reg['contrib'] / fullboxcontrib
        if fullboxAR > 0.2:
            p_w1 = W1 / ( W1 + W2 )
            rasampler1 = _gen_ra_sampler(mp,
                                         sofa_at_b1(r_a1),sofa_at_b1(r_a2),
                                         lin_s_at_b1)
            rasampler2 = _gen_ra_sampler(mp,
                                         sofa_at_b2(r_a1),sofa_at_b2(r_a2),
                                         lin_s_at_b2)
            ireg = 0
            while ireg < nreg:
                is_b1edge = rng()<=p_w1
                rasampler = rasampler1 if is_b1edge else rasampler2
                ra = rasampler(rng)
                a = r_a1*(1-ra)+r_a2*ra
                #Triangle which has non-zero height at chosen side and 0 at
                #the far side:
                rb = min(rng(),rng()) if is_b1edge else max(rng(),rng())
                b = b1*(1-rb)+b2*rb
                if (a-b)**2 <= 4*a*e:
                    ireg += 1
                    result.append( (a,b) )
        else:

            #A single box overlay not efficient enough! Sample a based on
            #contrib_at_a function, using uniform overlay sampling in a number of
            #cells. Cells granularity is adaptive to high level of change.
            asampler = AlphaRangeSampler(contrib_at_alpha = reg['contrib_at_a'],
                                         mp = mp, a1 = r_a1, a2 = r_a2,
                                         E_div_kT = e, b1=b1, b2=b2)
            bl,bu = b1, b2
            for i in range(nreg):
                a = asampler.sample(rng)
                assert a>=r_a1
                assert a<=r_a2
                if clip_betaminus or clip_betaplus:
                    bwidth = 2*mp.sqrt(e*a)
                    if clip_betaminus:
                        bl = max(b1,a-bwidth)
                    if clip_betaplus:
                        bu = min(b2,a+bwidth)
                sb1 = sofa_at_b1(a)
                sb2 = sofa_at_b2(a)
                rbl = (bl-b1)*invdb
                rbu = (bu-b1)*invdb
                sbl = (1-rbl)*sb1+rbl*sb2
                sbu = (1-rbu)*sb1+rbu*sb2
                while True:
                    b = bl+rng()*(bu-bl)
                    rb = (b-b1)*invdb
                    s = (1-rb)*sb1+rb*sb2
                    if max(sbl,sbu)*rng()<s:
                        break
                assert (a-b)**2 <= 4*a*e*(1+mp.mpf(1e-40))
                result.append( (a,b) )
    return result

def _gen_ra_sampler(mp, s1, s2, lin_s):
    if s1==s2:
        def ra(rng):
            return rng()
    if not lin_s:
        assert min(s1,s2)>0.0
        #log-lin
        k = s2/s1
        lnk = mp.log(k)
        if abs(lnk)<1e-10:#FIXME: Depend on dps?
            c1 = (lnk**2 + 3*lnk + 6)/6
            c2 = - (lnk**2 + lnk)/2
            c3 = lnk**2/3
            def ra(rng):
                R=rng()
                return R*(c1+R*(c2+R*c3))
            return ra
        else:
            invlnk = 1/lnk
            km1 = k-1
            def ra(rng):
                return mp.log(1+rng()*km1)*invlnk
            return ra
    else:
        #lin:
        return _gen_triangle01_sampler(mp,s1,s2)

def _gen_triangle01_sampler(mp,s1,s2):
    #sample according to pdf which is linear function between (0,s1) and (1,s2)
    if s1==s2:
        return lambda rng : rng()
    mins = min(s1,s2)
    smid = (s1+s2)/2
    punif = mins/smid
    s1largest = s1 > s2
    def sampler(rng):
        if rng()<punif:
            #uniform base:
            return rng()
        #triangle part:
        r = mp.sqrt(rng())
        return 1-r if s1largest else r
    return sampler

class AlphaRangeSampler:

    def __init__(self, mp, contrib_at_alpha, a1, a2, E_div_kT, b1, b2 ):
        #Initial subdivision:
        a1 = mp.mpf(a1)
        a2 = mp.mpf(a2)
        b1 = mp.mpf(b1)
        b2 = mp.mpf(b2)
        e = mp.mpf(E_div_kT)
        avals = mp.linspace(a1,a2,65)
        f = contrib_at_alpha
        def dev( _f1, _f2 ):
            if min(_f1,_f2)==0:
                return mp.mpf(0)#linear triangle => don't divide cell further
            return mp.mpf(0) if _f1 == _f2 else 2*abs(_f1-_f2)/(_f1+_f2)
        cells = []
        for i in range(len(avals)-1):
            ra1 = avals[i]
            ra2 = avals[i+1]
            f1, f2 = f(ra1), f(ra2)
            assert f1>=0.0
            assert f2>=0.0
            cells.append( (dev(f1,f2),f1,f2,ra1,ra2) )
        #Keep subdividing until largest dev is less than devthr:
        devthr = 0.6#tested with 0.4...0.9, all seems ok.
        ncellslimit = 2000
        while True:
            cells.sort()
            c=cells[-1]
            if c[0]<devthr:
                break
            #divide this cell:
            _,f1,f2,ra1,ra2 = c
            ramid = (ra1+ra2)/2
            fmid = f(ramid)
            assert fmid>=0.0
            cells.pop()
            cells.append( (dev(f1,fmid),f1,fmid,ra1,ramid) )
            cells.append( (dev(fmid,f2),fmid,f2,ramid,ra2) )
            assert len(cells)<ncellslimit

        #special points with possible extrema:
        special_avals = []
        mpf = mp.mpf
        sqrte=mp.sqrt(e)
        da1 = 2*mp.sqrt(abs(e*(b1+e)))
        da2 = 2*mp.sqrt(abs(e*(b2+e)))
        for sa in [ e, b1/mpf(3), b2/mpf(3),
                    ( mp.sqrt( abs(e + b1) ) - sqrte )**2,
                    ( mp.sqrt( abs(e + b2) ) - sqrte )**2,
                    2*e+b1-da1,
                    2*e+b1+da1,
                    2*e+b2-da2,
                    2*e+b2+da2 ]:
            if a1 < sa < a2:
                special_avals.append(sa)

        #Now, prepare uniform overlay sampler for these cells:
        overlay_contrib = []
        totsum = mp.mpf(0)
        safety = mpf('1.1')
        finalcells = []
        for c in cells:
            _,f1,f2,ra1,ra2 = c
            assert ra2-ra1 > 0.0
            overlayf = max(f1,f2)
            assert overlayf >= 0.0
            for sa in special_avals:
                if ra1 < sa < ra2:
                    overlayf = max(overlayf,f(sa))
            overlayf *= safety
            cntb = (ra2-ra1)*overlayf
            assert cntb >= 0.0
            finalcells.append( (overlayf,ra1,ra2) )
            totsum += cntb
            overlay_contrib.append( totsum )
        self.__cells = finalcells
        self.__overlay_contrib = overlay_contrib
        self.__f = f

    def sample( self, rng ):
        #Pick cell (bisect_left returns index of first entry not below Rselect):
        from bisect import bisect_left
        while True:
            Rselect = rng() * self.__overlay_contrib[-1]
            idx = bisect_left(self.__overlay_contrib, Rselect)
            c = self.__cells[idx]
            overlayf, ra1, ra2 = c
            da = ra2-ra1
            a = min(ra2,ra1 + rng()*da)
            fa = self.__f(a)
            assert fa <= overlayf
            if overlayf*rng() <= fa:
                return a
