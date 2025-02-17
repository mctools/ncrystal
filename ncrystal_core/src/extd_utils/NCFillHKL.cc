////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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
  //map keys used during search for hkl families.
  typedef unsigned FamKeyType;
  FamKeyType keygen(double fsq, double dsp) {
    //for most typical fsquared and dsp numbers, this will encode 3 significant
    //digits of fsq and dsp and their exponents into 4 separate bit areas. In
    //princple, the exponents could "overflow" internally for some obscure cases,
    //but that should be so rare as to not ruin our performance (the important
    //point is that (fsq,dsp) pairs that differ more than O(1e-3) should get
    //different keys, but that otherwise they should almost always not):
    nc_assert(fsq>0.0&&dsp>0.0);
    const int exponent_fsq = std::ceil(std::log10(fsq));
    const double mantissa_fsq = fsq * pow(10,-exponent_fsq);
    const int exponent_dsp = std::ceil(std::log10(dsp));
    const double mantissa_dsp = dsp * pow(10,-exponent_dsp);
    return unsigned(mantissa_fsq*1000+0.5)*4000000
      + unsigned(mantissa_dsp*1000+0.5)*4000
      + ncmax(3000+exponent_fsq*30 + exponent_dsp,0) ;
  }
}

#if ( defined(__clang__) || !defined(__GNUC__) || __GNUC__ >= 5 )
//Use mempool in C++11 (but not gcc 4.x.y)
#  define NCRYSTAL_NCMAT_USE_MEMPOOL
#endif

#ifdef NCRYSTAL_NCMAT_USE_MEMPOOL
#include <stdexcept>
#include <type_traits>
#include <scoped_allocator>
namespace NCRYSTAL_NAMESPACE {
  //We use a simple expanding-only memory pool for the temporary multimap used
  //to detect hkl families. This results in better cache locality and should
  //hopefully reduce memory fragmentation. TODO: Consider moving this pool
  //infrastructure to common utilities.
  class MemPool {
  public:
    explicit MemPool(std::size_t s) : m_size(s), m_offset(s+1) { m_chunks.reserve(64); }
    MemPool(MemPool const &) = delete;
    MemPool & operator=(MemPool const &) = delete;
    ~MemPool() { for (auto& e: m_chunks) ::operator delete(e); }
    void deallocate(void *, std::size_t) {}//ignore
    void * allocate(std::size_t n, std::size_t alignment)
    {
      nc_assert(n&&n<m_size);
      m_offset = ((m_offset + alignment - 1) / alignment ) * alignment;//move up to alignment
      if (m_offset + n > m_size) {//must grow
        m_chunks.push_back(m_data = static_cast<unsigned char *>(::operator new(m_size)));
        m_offset = 0;
      }
      void * result = m_data + m_offset;
      m_offset += n;
      return result;
    }
  private:
    unsigned char * m_data;
    std::size_t const m_size;
    std::size_t m_offset;
    std::vector<unsigned char*> m_chunks;
  };

  template <typename T>
  struct MemPoolAllocator {
    //Boiler-plate needed to make stl container use our MemPool:
    template <typename U> friend struct MemPoolAllocator;
    using value_type = T;
    using pointer = T *;
    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
    explicit MemPoolAllocator(MemPool * a) : m_pool(a) {}
    template <typename U> MemPoolAllocator(MemPoolAllocator<U> const & rhs) : m_pool(rhs.m_pool) {}
    pointer allocate(std::size_t n)
    {
      return static_cast<pointer>(m_pool->allocate(n * sizeof(T), alignof(T)));
    }
    void deallocate(pointer /*p*/, std::size_t /*n*/)
    {
      //if MemPool deallocate was not no-op, we should call it here: m_pool->deallocate(p, n * sizeof(T));
    }
    template <typename U> bool operator==(MemPoolAllocator<U> const & rhs) const { return m_pool == rhs.m_pool; }
    template <typename U> bool operator!=(MemPoolAllocator<U> const & rhs) const { return m_pool != rhs.m_pool; }
  private:
    MemPool * m_pool;
  };
  typedef std::multimap<FamKeyType, size_t, std::less<const FamKeyType>,
                        std::scoped_allocator_adaptor<MemPoolAllocator<std::pair<const FamKeyType, size_t>>>> FamMap;
}
#else
namespace NCRYSTAL_NAMESPACE {
  typedef std::multimap<FamKeyType,size_t> FamMap;
}
#endif

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
  }
  namespace detail {
    HKLList calculateHKLPlanesWithSymEqRefl( const StructureInfo&,
                                             const AtomInfoList&,
                                             FillHKLCfg,
                                             bool no_forceunitdebyewallerfactor );

    struct PreCalc {
      SmallVector<SmallVector<Vector,32>,4> atomic_pos;//atomic coordinates
      SmallVectD csl;//coherent scattering length
      SmallVectD msd;//mean squared displacement
      SmallVectD cache_factors;
      int max_h, max_k, max_l;
      SmallVectD whkl_thresholds;
      SmallVectD whkl;
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
      res.cache_factors.resize(res.csl.size(),0.0);

      //cache some thresholds for efficiency (see locations where it is used
      //for more comments):
      res.whkl_thresholds.reserve_hint(res.csl.size());
      for ( auto i : ncrange( res.csl.size() ) ) {
        if ( cfg.fsquarecut < 0.01 && cfg.fsquarecut > fsquarecut_lowest_possible_value )
          res.whkl_thresholds.push_back(std::log(ncabs(res.csl.at(i)) / cfg.fsquarecut ) );
        else
          res.whkl_thresholds.push_back(kInfinity);//use inf when not true that fsqcut^2 << fsq
      }

      res.whkl.resize(res.msd.size(),1.0);//init with unit factors in case of forceunitdebyewallerfactor

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
    return detail::calculateHKLPlanesWithSymEqRefl( structureInfo,
                                                    atomList,
                                                    std::move(cfg),
                                                    no_forceunitdebyewallerfactor );

  //For now we allow selection of a particular hkl value via an env var (a hacky
  //workarond required for certain validation plots - we should support this in
  //NCMatCfg instead).
  bool do_select = false;
  int select_h(0),select_k(0),select_l(0);
  std::string selecthklcfg = ncgetenv("FILLHKL_SELECTHKL");
  if (!selecthklcfg.empty()) {
    do_select = true;
    VectS parts;
    split(parts,selecthklcfg,0,',');
    nc_assert_always(parts.size()==3);
    select_h = str2int(parts.at(0));
    select_k = str2int(parts.at(1));
    select_l = str2int(parts.at(2));
  }

  const RotMatrix rec_lat = getReciprocalLatticeRot( structureInfo );

  auto cache = detail::fillHKLPreCalc( structureInfo, atomList, cfg);

  //We now conduct a brute-force loop over h,k,l indices, adding calculated info
  //in the following containers along the way:

  //Breaking O(N^2) complexity in compatibility searches by using map (the key
  //is an integer composed from Fsquared and d-spacing, and although clashes are
  //allowed, it should only clash rarely or efficiency is compromised):

  HKLList hkllist;
  if ( cache.whkl.empty() )
    return hkllist;//all elements have bcoh=0?

#ifdef NCRYSTAL_NCMAT_USE_MEMPOOL
  MemPool pool(10000000);
  MemPoolAllocator<void> poolalloc(&pool);
  FamMap fsq2hklidx(poolalloc);
#else
  FamMap fsq2hklidx;
#endif

  for( int loop_h=0;loop_h<=cache.max_h;++loop_h ) {
    for( int loop_k=(loop_h?-cache.max_k:0);loop_k<=cache.max_k;++loop_k ) {
      for( int loop_l=-cache.max_l;loop_l<=cache.max_l;++loop_l ) {
        const Vector hkl(loop_h,loop_k,loop_l);

        //calculate waveVector, wave number and dspacing:
        Vector waveVector = rec_lat*hkl;
        const double ksq = waveVector.mag2();
        if ( !valueInInterval(cache.ksq_preselect_interval,ksq))
          continue;

        if (no_forceunitdebyewallerfactor) nclikely {
          fillHKL_getWhkl(cache.whkl, ksq, cache.msd);
        }

        //calculate |F|^2
        double real_or_imag_upper_limit(0.0);
        for( unsigned i=0; i < cache.whkl.size(); ++i ) {
          if ( cache.whkl[i] > cache.whkl_thresholds[i]) {
            cache.cache_factors[i] = 0.0;
            continue;//Abort early to save exp/cos/sin calls. Note that
                     //O(fsquarecut) here corresponds to O(fsquarecut^2)
                     //contributions to final FSquared - for which we demand
                     //>fsquarecut below. We only do this when fsquarecut<1e-2
                     //(see calculations for whkl_thresholds above).
          } else {
            double factor = cache.csl[i]*std::exp(-cache.whkl[i]);
            cache.cache_factors[i] = factor;
            //Assuming cos(phase)*factor=sin(phase)*factor=|factor| gives us a cheap upper limit on
            //fsquared:
            real_or_imag_upper_limit += cache.atomic_pos[i].size()*ncabs( factor );
          }
        }

        //If the upper limit on fsq is below fsquarecut, we can skip already and
        //avoid needless calculations further down:
        if(real_or_imag_upper_limit*real_or_imag_upper_limit*2.0<cfg.fsquarecut)
          continue;

        if ( loop_h==0 && loop_k==0 && loop_l<=0)
          continue;

        if ( do_select && (loop_h!=select_h||loop_k!=select_k||loop_l!=select_l) )
            continue;

        //Time to calculate phases and sum up contributions. Use numerically
        //stable summation, for better results on low-symmetry crystals (the
        //main cost here is anyway the phase calculations, not the summation):
        StableSum real, imag;
        for( unsigned i=0 ; i < cache.whkl.size(); ++i ) {
          double factor = cache.cache_factors[i];
          if (!factor)
            continue;
          StableSum cpsum, spsum;
          for ( auto& pos : cache.atomic_pos[i] ) {
            //Phase is hkl.dot(pos)*2pi. We speed up the expensive
            //calculation of sin+cos by a factor of 3 by shifting the phase to
            //[0,2pi] (easily done by simply NOT multiplying with 2pi) and using
            //our own fast sincos_02pi through sincos_2pix. Since typically 99%
            //of the hkl initialisation time is spent calculating sin+cos here,
            //that actually translates into an overall speedup of a factor of 3
            //(measured in NCrystal v2.7.0)!
            const double phase_div2pi = hkl.dot(pos);
            auto spcp = sincos_2pix(phase_div2pi);
            cpsum.add(spcp.cos);
            spsum.add(spcp.sin);
          }
          real.add(cpsum.sum() * factor);
          imag.add(spsum.sum() * factor);
        }

        const double FSquared = ncsquare( real.sum() ) + ncsquare( imag.sum() );

        //skip weak or impossible reflections:
        if(FSquared<cfg.fsquarecut)
          continue;

        //Calculate d-spacing and recheck cut:
        const double kval = std::sqrt( ksq );
        const double invkval = 1.0 / kval;
        const double dspacing = k2Pi * invkval;

        if ( !valueInInterval( cache.dcut_interval, dspacing ) )
          continue;

        //Key for our fsq2hklidx multimap:
        FamKeyType searchkey(keygen(FSquared,dspacing));

        FamMap::iterator itSearchLB = fsq2hklidx.lower_bound(searchkey);
        FamMap::iterator itSearch(itSearchLB), itSearchE(fsq2hklidx.end());
        bool isnewfamily = true;
        for ( ; itSearch!=itSearchE && itSearch->first == searchkey; ++itSearch ) {
          nc_assert(itSearch->second<hkllist.size());
          HKLInfo& hi = hkllist[itSearch->second];
          if ( ncabs(FSquared-hi.fsquared) < cfg.merge_tolerance*(FSquared+hi.fsquared )
               && ncabs(dspacing-hi.dspacing) < cfg.merge_tolerance*(dspacing+hi.dspacing ) )
            {
              //Compatible with existing family, simply add HKL point to it.
              hi.multiplicity += 2;
              nc_assert(hi.explicitValues->list.has_value<std::vector<HKL>>());
              hi.explicitValues->list.get<std::vector<HKL>>().emplace_back(loop_h,loop_k,loop_l);
              isnewfamily = false;
              break;
            }
        }
        if (isnewfamily) {
          //Not fitting in existing group, set up new.
          if ( hkllist.size()>1000000 && !env_ignorefsqcut )//guard against crazy setups
            NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach"
                            " dcutoff = "<<cfg.dcutoff<<" Aa (you can try"
                            " to increase the target value with the dcutoff"
                            " parameter)");
          HKLInfo hi;
          hi.hkl = HKL{ loop_h, loop_k, loop_l };
          hi.multiplicity = 2;
          hi.fsquared = FSquared;
          hi.dspacing = dspacing;
          hi.explicitValues = std::make_unique<HKLInfo::ExplicitVals>();
          hi.explicitValues->list.emplace<std::vector<HKL>>();
          hi.explicitValues->list.get<std::vector<HKL>>().reserve(24);//shrinked below
          hi.explicitValues->list.get<std::vector<HKL>>().emplace_back(loop_h,loop_k,loop_l);
          fsq2hklidx.insert(itSearchLB,FamMap::value_type(searchkey,hkllist.size()));
          hkllist.emplace_back(std::move(hi));
        }
      }//loop_l
    }//loop_k
  }//loop_h

  //Sort explicit HKL entries and use first as representative index:
  for ( auto& hi : hkllist ) {
    auto& v = hi.explicitValues->list.get<std::vector<HKL>>();
    std::sort(v.begin(),v.end());
    v.shrink_to_fit();
    hi.hkl = v.front();
  }

  //NB: Not sorting by dspace (InfoBuilder will anyway do it and it is slightly
  //complicated to do consistently).
  hkllist.shrink_to_fit();
  return hkllist;
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
      SymHKLSeenTracker( double dcutoff ) : m_seen(std::make_unique<FastArraySmall>()), m_dcutoff(dcutoff) {}
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
            m_seenLarge = std::make_unique<FastArrayLarge>();
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
          NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach"
                          " dcutoff = "<<m_dcutoff<<" Aa (you can try"
                          " to increase the target value with the dcutoff"
                          " parameter)");
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

NC::HKLList NC::detail::calculateHKLPlanesWithSymEqRefl( const StructureInfo& structureInfo,
                                                         const AtomInfoList& atomList,
                                                         FillHKLCfg cfg,
                                                         bool no_forceunitdebyewallerfactor )
{
  nc_assert_always(structureInfo.spacegroup!=0);

  const bool env_ignorefsqcut = ncgetenv_bool("FILLHKL_IGNOREFSQCUT");
  nc_assert( !env_ignorefsqcut || cfg.fsquarecut == 0.0 );//due to logic in calling function

  const RotMatrix rec_lat = getReciprocalLatticeRot( structureInfo );
  EqRefl sym(structureInfo.spacegroup);
  auto sym_findrepval = [&sym]( int hh, int kk, int ll )
  {
    //NB: Tried to get eqv hkl with smallest min(|h|,|k|,|l|) instead to avoid
    //using the large cache in symSeenTracker, but profiling showed this to
    //cause a slowdown of 50% over the entire data library (in both rel and dbg
    //builds!).
    return sym.getEquivalentReflectionsRepresentativeValue(hh,kk,ll);
  };

  SymHKLSeenTracker symSeenTracker( env_ignorefsqcut ? -1.0 : cfg.dcutoff );

  //Make sure to always skip the (0,0,0) group:
  symSeenTracker.isFirstCheck(sym_findrepval(0,0,0));

  //For now we allow selection of a particular hkl value via an env var (a hacky
  //workarond required for certain validation plots - we should support this in
  //NCMatCfg instead).
  Optional<HKL> do_select;
  std::string selecthklcfg = ncgetenv("FILLHKL_SELECTHKL");
  if ( !selecthklcfg.empty() ) {
    VectS parts;
    split(parts,selecthklcfg,0,',');
    nc_assert_always(parts.size()==3);
    do_select = sym_findrepval( str2int(parts.at(0)),
                                str2int(parts.at(1)),
                                str2int(parts.at(2)) );
    //Pure efficiency improvement, mark *some* of the other values as already
    //seen (not all, due to the std::set fallback in symSeenTrackar):
    symSeenTracker.optimiseForSelection( do_select.value() );
  }

  auto cache = detail::fillHKLPreCalc( structureInfo, atomList, cfg);

  HKLList hkllist;
  if ( cache.whkl.empty() )
    return hkllist;//all elements have bcoh=0?
  //hkllist.reserve( 4096 );

  //We now conduct a brute-force loop over h,k,l indices, adding calculated info
  //in the following containers along the way. For reasons of symmetry we ignore
  //roughly half (but not all since the sym_key's might have sign flips).

  for( int loop_h = 0 ; loop_h <= cache.max_h; ++loop_h ) {
    for( int loop_k = (loop_h?-cache.max_k:0); loop_k <= cache.max_k; ++loop_k ) {
      for( int loop_l = -cache.max_l; loop_l <= cache.max_l; ++loop_l ) {

        auto sym_key = sym_findrepval( loop_h, loop_k, loop_l );
        if (!symSeenTracker.isFirstCheck(sym_key))
          continue;//Already seen this sym_key once.

        //calculate waveVector at the cost of a matrix multiplication, and
        //preselect on its squared magnitude:
        const Vector hkl(sym_key.h,sym_key.k,sym_key.l);
        Vector waveVector = rec_lat*hkl;
        const double ksq = waveVector.mag2();
        if ( ! valueInInterval( cache.ksq_preselect_interval , ksq ) )
          continue;

        if (no_forceunitdebyewallerfactor) nclikely {
          fillHKL_getWhkl( cache.whkl, ksq, cache.msd);
        }

        //calculate |F|^2
        double real_or_imag_upper_limit(0.0);
        for( unsigned i=0; i < cache.whkl.size(); ++i ) {
          if ( cache.whkl[i] > cache.whkl_thresholds[i]) {
            cache.cache_factors[i] = 0.0;
            continue;//Abort early to save exp/cos/sin calls. Note that
                     //O(fsquarecut) here corresponds to O(fsquarecut^2)
                     //contributions to final FSquared - for which we demand
                     //>fsquarecut below. We only do this when fsquarecut<1e-2
                     //(see calculations for whkl_thresholds above).
          } else {
            double factor = cache.csl[i]*std::exp(-cache.whkl[i]);
            cache.cache_factors[i] = factor;
            //Assuming cos(phase)*factor=sin(phase)*factor=|factor| gives us a cheap upper limit on
            //fsquared:
            real_or_imag_upper_limit += cache.atomic_pos[i].size()*ncabs( factor );
          }
        }

        //If the upper limit on fsq is below fsquarecut, we can skip already and
        //avoid needless calculations further down:
        if(real_or_imag_upper_limit*real_or_imag_upper_limit*2.0<cfg.fsquarecut)
          continue;

        //Time to calculate phases and sum up contributions. Use numerically
        //stable summation, for better results on low-symmetry crystals (the
        //main cost here is anyway the phase calculations, not the summation):
        StableSum real, imag;
        for( unsigned i=0 ; i < cache.whkl.size(); ++i ) {
          double factor = cache.cache_factors[i];
          if (!factor)
            continue;
          StableSum cpsum, spsum;
          for ( auto& pos : cache.atomic_pos[i] ) {
            //Phase is hkl.dot(pos)*2pi. We speed up the expensive
            //calculation of sin+cos by a factor of 3 by shifting the phase to
            //[0,2pi] (easily done by simply NOT multiplying with 2pi) and using
            //our own fast sincos_02pi through sincos_2pix. Since typically 99%
            //of the hkl initialisation time is spent calculating sin+cos here,
            //that actually translates into an overall speedup of a factor of 3
            //(measured in NCrystal v2.7.0)!
            const double phase_div2pi = hkl.dot(pos);
            auto spcp = sincos_2pix(phase_div2pi);
            cpsum.add(spcp.cos);
            spsum.add(spcp.sin);
          }
          real.add(cpsum.sum() * factor);
          imag.add(spsum.sum() * factor);
        }

        const double FSquared = ncsquare( real.sum() ) + ncsquare( imag.sum() );

        //skip weak or impossible reflections:
        if(FSquared<cfg.fsquarecut)
          continue;

        //Calculate d-spacing and recheck cut:
        const double dspacing = k2Pi / std::sqrt( ksq );

        if ( !valueInInterval( cache.dcut_interval, dspacing ) )
          continue;

        if ( do_select.has_value() && !(sym_key == do_select.value()) )
            continue;

        if ( hkllist.size()> 1000000 && !env_ignorefsqcut )//guard against crazy setups
          NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach"
                          " dcutoff = "<<cfg.dcutoff<<" Aa (you can try"
                          " to increase the target value with the dcutoff"
                          " parameter)");

        if ( hkllist.size() == decltype(hkllist)::nsmall+1 )
          hkllist.reserve_hint( 4096 );

        hkllist.emplace_back();
        auto& entry = hkllist.back();
        entry.dspacing = dspacing;
        entry.fsquared = FSquared;
        auto sym_list = sym.getEquivalentReflections( sym_key );
        entry.hkl = sym_list.front();
        entry.multiplicity = sym_list.size() * 2;
      }//loop_l
    }//loop_k
  }//loop_h

  //NB: Not sorting by dspace (InfoBuilder will anyway do it and it is slightly
  //complicated to do consistently).

  hkllist.shrink_to_fit();
  return hkllist;
}
