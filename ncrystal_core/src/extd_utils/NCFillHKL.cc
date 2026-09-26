////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/extd_utils/NCFillHKL.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include <bitset>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    using SmallVectD = SmallVector<double,64>;//64 atomic positions is usually (but not always) enough.

    inline void fillHKL_getWhkl(SmallVectD& out_whkl, const double ksq, const SmallVectD & msd)
    {
      nc_assert( msd.size() == out_whkl.size() );
      //Sears, Acta Cryst. (1997). A53, 35-45
      double kk2 = 0.5*ksq;

      //NB: Do not call out_whkl.clear() followed by push_back, as the usage of
      //SmallVector means we will discard the large allocation and have to do
      //constant reallocations!!

      auto it = msd.begin();
      auto itE = msd.end();
      auto itOut = out_whkl.begin();
      for( ; it!=itE; ++it )
        *itOut++ = kk2*(*it);
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    constexpr const double fsquarecut_lowest_possible_value = 1.0e-300;

    [[noreturn]] void throwCombinatoricsTooGreat( double dcutoff )
    {
      NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach"
                      " dcutoff = "<<dcutoff<<" Aa (you can try"
                      " to increase the target value with the dcutoff"
                      " parameter)");
    }

    //For now we allow selection of a particular hkl value via an env var (a
    //hacky workarond required for certain validation plots - we should support
    //this in NCMatCfg instead):
    Optional<HKL> selectedHKLFromEnv()
    {
      Optional<HKL> res;
      std::string selecthklcfg = ncgetenv("FILLHKL_SELECTHKL");
      if ( !selecthklcfg.empty() ) {
        VectS parts;
        split(parts,selecthklcfg,0,',');
        nc_assert_always(parts.size()==3);
        res = HKL( str2int(parts.at(0)),
                   str2int(parts.at(1)),
                   str2int(parts.at(2)) );
      }
      return res;
    }

    struct PreCalc {
      SmallVector<SmallVector<Vector,32>,4> atomic_pos;//atomic coordinates
      SmallVectD csl;//coherent scattering length
      SmallVectD msd;//mean squared displacement
      int max_h, max_k, max_l;
      SmallVectD whkl_thresholds;
      PairDD ksq_preselect_interval;
      PairDD dcut_interval;
    };

    PreCalc fillHKLPreCalc( const StructureInfo& si,
                            const AtomInfoList& atomList,
                            const FillHKLCfg& cfg)
    {
      PreCalc res;
      for ( auto& ai : atomList ) {
        nc_assert( ai.msd().has_value() );
        if ( ! ( ncabs( ai.atomData().coherentScatLen() ) > 0.0 ) )
          continue;//ignore "sterile" species
        res.msd.push_back( ai.msd().value() );
        res.csl.push_back( ai.atomData().coherentScatLen() );
        SmallVector<Vector,32> pos;
        pos.reserve_hint( ai.unitCellPositions().size() );
        for ( const auto& p : ai.unitCellPositions() )
          pos.push_back( p.as<Vector>() );
        res.atomic_pos.push_back( std::move(pos) );
      }

      {
        auto max_hkl = estimateHKLRange( cfg.dcutoff,
                                         si.lattice_a, si.lattice_b, si.lattice_c,
                                         si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
        res.max_h = max_hkl.h;
        res.max_k = max_hkl.k;
        res.max_l = max_hkl.l;
      }

      nc_assert_always(res.msd.size()==res.atomic_pos.size());
      nc_assert_always(res.msd.size()==res.csl.size());

      //cache some thresholds for efficiency (see locations where it is used
      //for more comments):
      res.whkl_thresholds.reserve_hint(res.csl.size());
      for ( auto i : ncrange( res.csl.size() ) ) {
        if ( cfg.fsquarecut < 0.01 && cfg.fsquarecut > fsquarecut_lowest_possible_value )
          res.whkl_thresholds.push_back(std::log(ncabs(res.csl.at(i)) / cfg.fsquarecut ) );
        else
          res.whkl_thresholds.push_back(kInfinity);//use inf when not true that fsqcut^2 << fsq
      }

      auto clampNormal = [](double x)
      {
        //valueInInterval might trigger FPE if used with infinity
        return ncclamp( x, std::numeric_limits<double>::min(), std::numeric_limits<double>::max() );
      };

      //Acceptable range of ksq=(2pi/dspacing)^2 and dspacing (ksq range expanded
      //slightly to avoid removing too much - the real check is on dspacing and is
      //performed later):
      res.ksq_preselect_interval = { clampNormal( (k4PiSq*(1.0-1e-14)) / ncsquare(cfg.dcutoffup) ),
                                     clampNormal( (k4PiSq*(1.0+1e-14)) / ncsquare(cfg.dcutoff) ) };
      res.dcut_interval = { clampNormal(cfg.dcutoff), clampNormal(cfg.dcutoffup) };
      return res;
    }

    class FSquaredCalc {
      //Calculates |F|^2 for HKL points. First call setKSq with the squared
      //length of the wave vector (i.e. the d-spacing) of the HKL point(s),
      //which updates the Debye-Waller factors, and then call calc(..) for any
      //number of HKL points with that ksq value (e.g. symmetry-equivalent
      //ones).
    public:
      FSquaredCalc( const PreCalc& pc, double fsquarecut,
                    bool no_forceunitdebyewallerfactor )
        : m_pc(pc),
          m_fsquarecut(fsquarecut),
          m_use_dw(no_forceunitdebyewallerfactor)
      {
        m_cache_factors.resize(m_pc.csl.size(),0.0);
        //init with unit factors in case of forceunitdebyewallerfactor:
        m_whkl.resize(m_pc.msd.size(),1.0);
      }

      bool empty() const { return m_whkl.empty(); }//all elements have bcoh=0?

      //Returns false if |F|^2 is guaranteed to be below fsquarecut:
      bool setKSq( double ksq )
      {
        if (m_use_dw) nclikely {
          fillHKL_getWhkl(m_whkl, ksq, m_pc.msd);
        }
        double real_or_imag_upper_limit(0.0);
        for( unsigned i=0; i < m_whkl.size(); ++i ) {
          if ( m_whkl[i] > m_pc.whkl_thresholds[i]) {
            m_cache_factors[i] = 0.0;
            continue;//Abort early to save exp/cos/sin calls. Note that
                     //O(fsquarecut) here corresponds to O(fsquarecut^2)
                     //contributions to final FSquared - for which we demand
                     //>fsquarecut below. We only do this when fsquarecut<1e-2
                     //(see calculations for whkl_thresholds above).
          } else {
            double factor = m_pc.csl[i]*std::exp(-m_whkl[i]);
            m_cache_factors[i] = factor;
            //Assuming cos(phase)*factor=sin(phase)*factor=|factor| gives us a
            //cheap upper limit on fsquared:
            real_or_imag_upper_limit += m_pc.atomic_pos[i].size()*ncabs( factor );
          }
        }
        //If the upper limit on fsq is below fsquarecut, we can skip already and
        //avoid needless calculations further down:
        return !(real_or_imag_upper_limit*real_or_imag_upper_limit*2.0<m_fsquarecut);
      }

      double calc( const Vector& hkl ) const
      {
        //Time to calculate phases and sum up contributions. Use numerically
        //stable summation, for better results on low-symmetry crystals (the
        //main cost here is anyway the phase calculations, not the summation):
        StableSum real, imag;
        for( unsigned i=0 ; i < m_whkl.size(); ++i ) {
          double factor = m_cache_factors[i];
          if (!factor)
            continue;
          StableSum cpsum, spsum;
          for ( auto& pos : m_pc.atomic_pos[i] ) {
            //Phase is hkl.dot(pos)*2pi. We speed up the expensive calculation
            //of sin+cos by a factor of 3 by shifting the phase to [0,2pi]
            //(easily done by simply NOT multiplying with 2pi) and using our
            //own fast sincos_02pi through sincos_2pix. Since typically 99% of
            //the hkl initialisation time is spent calculating sin+cos here,
            //that actually translates into an overall speedup of a factor of
            //3 (measured in NCrystal v2.7.0)!
            const double phase_div2pi = hkl.dot(pos);
            auto spcp = sincos_2pix(phase_div2pi);
            cpsum.add(spcp.cos);
            spsum.add(spcp.sin);
          }
          real.add(cpsum.sum() * factor);
          imag.add(spsum.sum() * factor);
        }
        return ncsquare( real.sum() ) + ncsquare( imag.sum() );
      }

      //Upper limit on |F|^2 (all contributions in phase, no Debye-Waller
      //factors):
      double fsquaredUpperLimit() const
      {
        double f = 0.0;
        for ( auto i : ncrange( m_pc.csl.size() ) )
          f += ncabs( m_pc.csl[i] ) * m_pc.atomic_pos[i].size();
        return ncsquare( f );
      }

    private:
      const PreCalc& m_pc;
      double m_fsquarecut;
      bool m_use_dw;
      SmallVectD m_cache_factors;
      SmallVectD m_whkl;
    };

    //Estimated number of hkl points (in the half-space) with dspacing in
    //[dcutoff,dcutoffup]: The reciprocal lattice has a density of V points
    //per Aa^-3 (d*=1/d convention), so this is V times half the volume of the
    //spherical shell with radii 1/dcutoffup and 1/dcutoff:
    double estimateNHKLPoints( const StructureInfo& si, const FillHKLCfg& cfg )
    {
      const double kmin = 1.0 / cfg.dcutoffup;
      const double kmax = 1.0 / cfg.dcutoff;
      return ( 2.0 * kPi / 3.0 ) * si.volume * ( kmax*kmax*kmax
                                                   - kmin*kmin*kmin );
    }

    HKLList calculateHKLPlanesNoSym( const StructureInfo& structureInfo,
                                     const AtomInfoList& atomList,
                                     const FillHKLCfg& cfg,
                                     bool no_forceunitdebyewallerfactor,
                                     bool env_ignorefsqcut )
    {
      //Without space group, planes are grouped into families by (FSquared,
      //dspacing) values. Accepted hkl points are buffered, and whenever the
      //buffer is full (and at the end) sorted by dspacing and merged into the
      //families in a single sweep. The grouping thus does not depend on the
      //loop order or any binning of values, and the memory used for buffering
      //is bounded (in practice a single buffer suffices for most crystals):
      const Optional<HKL> do_select = selectedHKLFromEnv();
      const RotMatrix rec_lat = getReciprocalLatticeRot( structureInfo );
      const auto precalc = fillHKLPreCalc( structureInfo, atomList, cfg );
      FSquaredCalc fsqcalc( precalc, cfg.fsquarecut,
                            no_forceunitdebyewallerfactor );

      HKLList hkllist;
      if ( fsqcalc.empty() )
        return hkllist;//all elements have bcoh=0?

      struct Pt { double d; double fsq; HKL hkl; };
      const std::size_t nbufmax = ( cfg.max_buffered_points
                                    ? cfg.max_buffered_points
                                    : (16u<<20) / sizeof(Pt) );
      std::vector<Pt> pts;
      {
        //Reserve for the estimated number of points (the F2 cut typically
        //removes some), or 0.5MB if the estimate is not available:
        const double est = estimateNHKLPoints( structureInfo, cfg );
        pts.reserve( static_cast<std::size_t>(
                       est > 0.0
                       ? ncmin( 1.05 * est + 64.0, double(nbufmax) )
                       : double( std::min<std::size_t>( nbufmax,
                                                   (1u<<19) / sizeof(Pt) ) ) ) );
      }

      //(reference dspacing, index in hkllist) of families from previous
      //buffers, sorted by dspacing:
      std::vector<std::pair<double,std::size_t>> famidx;

      //Families are matched against the values of their first member, but
      //finally get the average values of all members, calculated as first +
      //mean(member-first). The differences are exact (Sterbenz lemma), so the
      //result is in practice independent of the order of members, is exactly
      //the common value if all members agree, and never leaves their range:
      std::vector<PairDD> famdevsums;//(d,fsq) sums of member-first

      const double tol = cfg.merge_tolerance;
      auto isCompatible = [tol]( const Pt& p, const HKLInfo& hi )
      {
        return ( ncabs( p.d - hi.dspacing ) < tol * ( p.d + hi.dspacing )
                 && ncabs( p.fsq - hi.fsquared ) < tol * ( p.fsq + hi.fsquared ) );
      };

      auto mergeBuffer = [&]( bool last )
      {
        std::sort( pts.begin(), pts.end(),
                   []( const Pt& a, const Pt& b )
                   {
                     if ( a.d != b.d )
                       return a.d < b.d;
                     if ( a.fsq != b.fsq )
                       return a.fsq < b.fsq;
                     return a.hkl < b.hkl;
                   } );
        //Families created in this sweep have increasing dspacing, so those
        //compatible in dspacing with a given point are always a suffix of
        //hkllist, starting at iwin (which never decreases):
        const std::size_t nfam_prev = hkllist.size();
        std::size_t iwin = nfam_prev;
        for ( auto& p : pts ) {
          const std::size_t nofam = std::numeric_limits<std::size_t>::max();
          std::size_t ifam = nofam;
          if ( !famidx.empty() ) {
            //Families from previous buffers first (conservative search
            //window, isCompatible has the final say):
            auto it = std::lower_bound( famidx.begin(), famidx.end(),
                                        p.d * ( 1.0 - 3.0 * tol ),
                                        []( const std::pair<double,
                                            std::size_t>& e, double d )
                                        { return e.first < d; } );
            const double dmax = p.d * ( 1.0 + 3.0 * tol );
            for ( ; it != famidx.end() && it->first <= dmax; ++it ) {
              if ( isCompatible( p, hkllist[it->second] ) ) {
                ifam = it->second;
                break;
              }
            }
          }
          if ( ifam == nofam ) {
            while ( iwin < hkllist.size()
                    && !( ncabs( p.d - hkllist[iwin].dspacing )
                          < tol * ( p.d + hkllist[iwin].dspacing ) ) )
              ++iwin;
            for ( auto i : ncrange( iwin, hkllist.size() ) ) {
              if ( isCompatible( p, hkllist[i] ) ) {
                ifam = i;
                break;
              }
            }
          }
          if ( ifam != nofam ) {
            //Compatible with existing family, simply add HKL point to it.
            HKLInfo& hi = hkllist[ifam];
            hi.multiplicity += 2;
            hi.explicitValues->list.get<std::vector<HKL>>().push_back( p.hkl );
            famdevsums[ifam].first += p.d - hi.dspacing;
            famdevsums[ifam].second += p.fsq - hi.fsquared;
            continue;
          }
          //Not fitting in existing group, set up new.
          if ( hkllist.size()>1000000 && !env_ignorefsqcut )//guard against crazy setups
            throwCombinatoricsTooGreat( cfg.dcutoff );
          hkllist.emplace_back();
          HKLInfo& hi = hkllist.back();
          hi.multiplicity = 2;
          hi.fsquared = p.fsq;
          hi.dspacing = p.d;
          hi.explicitValues = ncmake_unique<HKLInfo::ExplicitVals>();
          hi.explicitValues->list.emplace<std::vector<HKL>>();
          hi.explicitValues->list.get<std::vector<HKL>>().push_back( p.hkl );
          famdevsums.emplace_back( 0.0, 0.0 );
        }
        pts.clear();
        if ( last )
          return;
        //Update famidx with the new families (already sorted by dspacing):
        const auto nidx_prev = famidx.size();
        for ( auto i : ncrange( nfam_prev, hkllist.size() ) )
          famidx.emplace_back( hkllist[i].dspacing, i );
        std::inplace_merge( famidx.begin(), famidx.begin() + nidx_prev,
                            famidx.end() );
      };

      //Brute-force loop over h,k,l indices (skipping half, since (h,k,l) and
      //(-h,-k,-l) are always in the same family):
      for( int loop_h=0;loop_h<=precalc.max_h;++loop_h ) {
        for( int loop_k=(loop_h?-precalc.max_k:0);loop_k<=precalc.max_k;++loop_k ) {
          for( int loop_l=-precalc.max_l;loop_l<=precalc.max_l;++loop_l ) {
            if ( loop_h==0 && loop_k==0 && loop_l<=0)
              continue;

            const Vector hkl(loop_h,loop_k,loop_l);

            //calculate waveVector, wave number and dspacing:
            Vector waveVector = rec_lat*hkl;
            const double ksq = waveVector.mag2();
            if ( !valueInInterval(precalc.ksq_preselect_interval,ksq))
              continue;

            if ( !fsqcalc.setKSq( ksq ) )
              continue;

            if ( do_select.has_value()
                 && !( HKL(loop_h,loop_k,loop_l) == do_select.value() ) )
              continue;

            const double FSquared = fsqcalc.calc( hkl );

            //skip weak or impossible reflections:
            if(FSquared<cfg.fsquarecut)
              continue;

            //Calculate d-spacing and recheck cut:
            const double dspacing = k2Pi / std::sqrt( ksq );

            if ( !valueInInterval( precalc.dcut_interval, dspacing ) )
              continue;

            pts.push_back( Pt{ dspacing, FSquared,
                               HKL( loop_h, loop_k, loop_l ) } );
            if ( pts.size() >= nbufmax )
              mergeBuffer( false );
          }//loop_l
        }//loop_k
      }//loop_h
      mergeBuffer( true );

      //Use average values:
      nc_assert_always( famdevsums.size() == hkllist.size() );
      for ( auto i : ncrange( hkllist.size() ) ) {
        HKLInfo& hi = hkllist[i];
        const double n = 0.5 * hi.multiplicity;
        hi.dspacing += famdevsums[i].first / n;
        hi.fsquared += famdevsums[i].second / n;
      }
      famdevsums.clear();
      famdevsums.shrink_to_fit();

      //Sort explicit HKL entries and use first as representative index:
      for ( auto& hi : hkllist ) {
        auto& v = hi.explicitValues->list.get<std::vector<HKL>>();
        std::sort(v.begin(),v.end());
        v.shrink_to_fit();
        hi.hkl = v.front();
      }

      //NB: Not sorting by dspace (InfoBuilder will anyway do it and it is
      //slightly complicated to do consistently).
      hkllist.shrink_to_fit();
      return hkllist;
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace {

    class SymHKLSeenTracker {
    private:
      static constexpr unsigned fast_small_C = 128;//128;//always enabled, uses 4*C^3 bits [C=128 gives 1.04MB]
      static constexpr unsigned fast_large_C = 512;//512;//rarely used, on-demand usage only [C=512 gives 67MB]
      static constexpr unsigned n_small = 4*fast_small_C*fast_small_C*fast_small_C;
      static constexpr unsigned n_large = 4*fast_large_C*fast_large_C*fast_large_C;
      using FastArraySmall = std::bitset<n_small>;
      using FastArrayLarge = std::bitset<n_large>;
      //Both bitsets on the stack (to prevent stack overflow), but the smaller one is always set up.
      std::unique_ptr<FastArraySmall> m_seen;//<--- this is the workhorse which is almost always used exclusively. Fast and not too big.
      std::unique_ptr<FastArrayLarge> m_seenLarge;
      std::set<HKL> m_seenFallBack;//<--- ultimate fallback, bad performance but always works.
      double m_dcutoff;//for err msg (-1 means err disabled)
    public:
      SymHKLSeenTracker( double dcutoff ) : m_seen(ncmake_unique<FastArraySmall>()), m_dcutoff(dcutoff) {}
      bool isFirstCheck( const HKL& hkl ) {
        auto idx = calcFastIdx<fast_small_C>(hkl);
        if ( idx.has_value() ) nclikely {
          auto e = (*m_seen)[idx.value()];
          if ( (bool)e )
            return false;
          e = true;
          return true;
        } else {
          return isFirstCheckFallBack(hkl);
        }
      }

      //Pretend that *some* of the values != v where already seen (not all, due
      //to internal storage being dynamic).
      void optimiseForSelection( const HKL& v )
      {
        m_seen->set();//sets all to true, pretending they were already processed
        auto idx = calcFastIdx<fast_small_C>(v);
        if ( idx.has_value() )
          m_seen->set(idx.value(),false);
      }
    private:
      bool isFirstCheckFallBack( const HKL& v ) {
        auto idx = calcFastIdx<fast_large_C>(v);
        if ( idx.has_value() ) {
          if (!m_seenLarge) ncunlikely {
            m_seenLarge = ncmake_unique<FastArrayLarge>();
          }
          auto e = (*m_seenLarge)[idx.value()];
          if ( (bool)e )
            return false;
          e = true;
          return true;
        }
        //Ultimate fallback:
        auto it_and_inserted = m_seenFallBack.insert(v);
        if ( m_seenFallBack.size() == 100000000 && m_dcutoff != -1.0 )
          throwCombinatoricsTooGreat( m_dcutoff );
        return it_and_inserted.second;
      }

      template<int C>
      Optional<unsigned> calcFastIdx( const HKL&v ) const
      {
        //NOTE: l varies most frequently in the calling loop, then k, then h. So
        //for cache-locality we should make sure that indices close in l are
        //close, etc. (this is particularly important if overspilling to the
        //m_seenLarge cache). Note on this note: The EqRefl remapping of HKL
        //values screws this up, but benchmarking still showed the code below to
        //be fastest.

        //Works if h in range 0..C-1 (C values), and k,l in range -(C-1)..C (2C values)
        static_assert(C>=2&&C<=10000,"");
        nc_assert( v.h >= 0 );
        constexpr int TwoC = 2*C;
        constexpr int Cm1 = (C-1);
        constexpr int mCm1 = -(C-1);
        Optional<unsigned> res;
        if ( v.h < C && std::min(v.k,v.l) >= mCm1 && std::max(v.k,v.l) <= C ) {
          nc_assert( Cm1 + v.k >= 0 && Cm1 + v.k < TwoC );
          nc_assert( Cm1 + v.l >= 0 && Cm1 + v.l < TwoC );
          //res = static_cast<unsigned>(v.h + C * ( ( Cm1 + v.k) +  TwoC * ( Cm1 + v.l) ));This way would be very slow
          res = static_cast<unsigned>( (Cm1 + v.l) + TwoC * ( ( Cm1 + v.k) + TwoC * v.h ) );//And this way much better
          nc_assert( res < 4*C*C*C );
        }
        return res;
      }
    };
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    HKLList calculateHKLPlanesWithSymEqRefl( const StructureInfo& structureInfo,
                                             const AtomInfoList& atomList,
                                             const FillHKLCfg& cfg,
                                             bool no_forceunitdebyewallerfactor,
                                             bool env_ignorefsqcut )
    {
      nc_assert_always(structureInfo.spacegroup!=0);
      //Caller sets fsquarecut to 0.0 and then clamps it:
      nc_assert( !env_ignorefsqcut
                 || cfg.fsquarecut <= fsquarecut_lowest_possible_value );

      const RotMatrix rec_lat = getReciprocalLatticeRot( structureInfo );
      EqRefl sym( structureInfo.spacegroup,
                  usesRhombohedralAxes( structureInfo.spacegroup,
                                        structureInfo.alpha ) );
      auto sym_findrepval = [&sym]( int hh, int kk, int ll )
      {
        //NB: Tried to get eqv hkl with smallest min(|h|,|k|,|l|) instead to
        //avoid using the large cache in symSeenTracker, but profiling showed
        //this to cause a slowdown of 50% over the entire data library (in both
        //rel and dbg builds!).
        return sym.getEquivalentReflectionsRepresentativeValue(hh,kk,ll);
      };

      SymHKLSeenTracker symSeenTracker( env_ignorefsqcut ? -1.0 : cfg.dcutoff );

      //Make sure to always skip the (0,0,0) group:
      symSeenTracker.isFirstCheck(sym_findrepval(0,0,0));

      Optional<HKL> do_select = selectedHKLFromEnv();
      if ( do_select.has_value() ) {
        do_select = sym_findrepval( do_select.value().h,
                                    do_select.value().k,
                                    do_select.value().l );
        //Pure efficiency improvement, mark *some* of the other values as
        //already seen (not all, due to the std::set fallback in
        //symSeenTrackar):
        symSeenTracker.optimiseForSelection( do_select.value() );
      }

      const auto precalc = fillHKLPreCalc( structureInfo, atomList, cfg );
      FSquaredCalc fsqcalc( precalc, cfg.fsquarecut,
                            no_forceunitdebyewallerfactor );

      HKLList hkllist;
      if ( fsqcalc.empty() )
        return hkllist;//all elements have bcoh=0?

      //Self-check that the structure has the symmetry assumed by EqRefl for
      //the space group, by calculating F2 for all planes of groups with low
      //hkl indices. Differences above 1% (or tiny F2) indicate e.g. a wrong
      //space group:
      const double f2max = fsqcalc.fsquaredUpperLimit();
      constexpr int symcheck_maxhkl = 4;

      //We now conduct a brute-force loop over h,k,l indices, adding calculated
      //info in the following containers along the way. For reasons of symmetry
      //we ignore roughly half (but not all since the sym_key's might have sign
      //flips).

      //Loose ksq preselection on the loop hkl, to skip the costly EqRefl
      //lookup for most hkl points outside the range. Equivalent hkl points
      //have identical ksq up to rounding errors, so the slack guarantees
      //identical results to preselecting on the representative hkl only (the
      //upper value may be DBL_MAX, avoid overflow):
      const double ksq_up = precalc.ksq_preselect_interval.second;
      const PairDD ksq_loose_interval
        = { precalc.ksq_preselect_interval.first * ( 1.0 - 1e-5 ),
            ksq_up > 1e300 ? ksq_up : ksq_up * ( 1.0 + 1e-5 ) };

      for( int loop_h = 0 ; loop_h <= precalc.max_h; ++loop_h ) {
        for( int loop_k = (loop_h?-precalc.max_k:0); loop_k <= precalc.max_k; ++loop_k ) {
          for( int loop_l = -precalc.max_l; loop_l <= precalc.max_l; ++loop_l ) {

            if ( ! valueInInterval( ksq_loose_interval,
                                    ( rec_lat * Vector( loop_h, loop_k,
                                                        loop_l ) ).mag2() ) )
              continue;

            auto sym_key = sym_findrepval( loop_h, loop_k, loop_l );
            if (!symSeenTracker.isFirstCheck(sym_key))
              continue;//Already seen this sym_key once.

            //calculate waveVector at the cost of a matrix multiplication, and
            //preselect on its squared magnitude:
            const Vector hkl(sym_key.h,sym_key.k,sym_key.l);
            Vector waveVector = rec_lat*hkl;
            const double ksq = waveVector.mag2();
            if ( ! valueInInterval( precalc.ksq_preselect_interval , ksq ) )
              continue;

            if ( !fsqcalc.setKSq( ksq ) )
              continue;

            const double FSquared = fsqcalc.calc( hkl );

            //skip weak or impossible reflections:
            if(FSquared<cfg.fsquarecut)
              continue;

            //Calculate d-spacing and recheck cut:
            const double dspacing = k2Pi / std::sqrt( ksq );

            if ( !valueInInterval( precalc.dcut_interval, dspacing ) )
              continue;

            if ( do_select.has_value() && !(sym_key == do_select.value()) )
              continue;

            if ( hkllist.size()> 1000000 && !env_ignorefsqcut )//guard against crazy setups
              throwCombinatoricsTooGreat( cfg.dcutoff );

            if ( hkllist.size() == decltype(hkllist)::nsmall+1 )
              hkllist.reserve_hint( 4096 );

            hkllist.emplace_back();
            auto& entry = hkllist.back();
            entry.dspacing = dspacing;
            entry.fsquared = FSquared;
            auto sym_list = sym.getEquivalentReflections( sym_key );
            entry.hkl = sym_list.front();
            entry.multiplicity = sym_list.size() * 2;
            if ( std::abs(sym_key.h) <= symcheck_maxhkl
                 && std::abs(sym_key.k) <= symcheck_maxhkl
                 && std::abs(sym_key.l) <= symcheck_maxhkl ) {
              for ( auto& e : sym_list ) {
                const double f2 = fsqcalc.calc( Vector( e.h, e.k, e.l ) );
                if ( ncabs( f2 - FSquared ) > 0.01*ncmax( FSquared, 1e-4*f2max ) )
                  NCRYSTAL_THROW2(BadInput,"Crystal structure is not consistent"
                                  " with space group "<<structureInfo.spacegroup
                                  <<" (the planes "<<sym_key.h<<","<<sym_key.k
                                  <<","<<sym_key.l<<" and "<<e.h<<","<<e.k<<","
                                  <<e.l<<" should be symmetry-equivalent but"
                                  " have |F|^2 values of "<<FSquared<<" and "
                                  <<f2<<" barn)");
              }
            }
          }//loop_l
        }//loop_k
      }//loop_h

      //NB: Not sorting by dspace (InfoBuilder will anyway do it and it is
      //slightly complicated to do consistently).

      hkllist.shrink_to_fit();
      return hkllist;
    }
  }
}

NC::HKLList NC::calculateHKLPlanes( const StructureInfo& structureInfo,
                                    const AtomInfoList& atomList,
                                    FillHKLCfg cfg )
{
  if ( atomList.empty() )
    NCRYSTAL_THROW(BadInput,"calculateHKLPlanes needs a non-empty AtomInfoList");
  for ( auto& ai : atomList ) {
    if ( !ai.msd().has_value() ) {
      //NB: strictly not needed if coherent scat len of that entry is vanishing,
      //but for now we keep the requirement of always needing msd just to be
      //consistent (we could reconsider this):
      NCRYSTAL_THROW(BadInput,"calculateHKLPlanes needs an AtomInfoList"
                     " which includes mean-squared-displacements of all atoms");
    }
  }

  nc_assert_always(cfg.dcutoff>0.0&&cfg.dcutoff<cfg.dcutoffup);

  const bool env_ignorefsqcut = ncgetenv_bool("FILLHKL_IGNOREFSQCUT");
  if (env_ignorefsqcut)
    cfg.fsquarecut = 0.0;

  if ( cfg.fsquarecut>=0.0 )
    cfg.fsquarecut = ncmax(cfg.fsquarecut,fsquarecut_lowest_possible_value);

  bool no_forceunitdebyewallerfactor;
  if ( cfg.use_unit_debye_waller_factor.has_value() ) {
    //Caller requested behaviour:
    no_forceunitdebyewallerfactor = ! cfg.use_unit_debye_waller_factor.value();
  } else {
    //Fall-back to global default behaviour (which can be modified with env
    //var for historic reasons):
    no_forceunitdebyewallerfactor = !(ncgetenv_bool("FILLHKL_FORCEUNITDEBYEWALLERFACTOR"));
  }

  if ( structureInfo.spacegroup != 0 )
    return calculateHKLPlanesWithSymEqRefl( structureInfo, atomList, cfg,
                                            no_forceunitdebyewallerfactor,
                                            env_ignorefsqcut );
  return calculateHKLPlanesNoSym( structureInfo, atomList, cfg,
                                  no_forceunitdebyewallerfactor,
                                  env_ignorefsqcut );
}
