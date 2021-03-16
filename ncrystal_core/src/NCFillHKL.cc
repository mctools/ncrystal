////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/internal/NCFillHKL.hh"
#include "NCrystal/internal/NCOrientUtils.hh"
#include "NCrystal/internal/NCRotMatrix.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/NCDefs.hh"
#include <cstdlib>

namespace NC = NCrystal;

namespace NCrystal {
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
namespace NCrystal {
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
namespace NCrystal {
  typedef std::multimap<FamKeyType,size_t> FamMap;
}
#endif

namespace NCrystal {
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

void NC::fillHKL( NC::Info& info, FillHKLCfg cfg )
{
  const bool env_ignorefsqcut = std::getenv("NCRYSTAL_FILLHKL_IGNOREFSQCUT");
  if (env_ignorefsqcut)
    cfg.fsquarecut = 0.0;

  constexpr const double fsquarecut_lowest_possible_value = 1.0e-300;
  if ( cfg.fsquarecut>=0.0 )
    cfg.fsquarecut = ncmax(cfg.fsquarecut,fsquarecut_lowest_possible_value);

  bool no_forceunitdebyewallerfactor;
  if ( cfg.use_unit_debye_waller_factor.has_value() ) {
    //Caller requested behaviour:
    no_forceunitdebyewallerfactor = ! cfg.use_unit_debye_waller_factor.value();
  } else {
    //Fall-back to global default behaviour (which can be modified with env
    //var for historic reasons):
    no_forceunitdebyewallerfactor = !(std::getenv("NCRYSTAL_FILLHKL_FORCEUNITDEBYEWALLERFACTOR"));
  }

  //For now we allow selection of a particular hkl value via an env var (a hacky
  //workarond required for certain validation plots - we should support this in
  //NCMatCfg instead).
  bool do_select = false;
  int select_h(0),select_k(0),select_l(0);
  const char * selecthklcfg = std::getenv("NCRYSTAL_FILLHKL_SELECTHKL");
  if (selecthklcfg) {
    do_select = true;
    VectS parts;
    split(parts,selecthklcfg,0,',');
    nc_assert_always(parts.size()==3);
    select_h = str2int(parts.at(0));
    select_k = str2int(parts.at(1));
    select_l = str2int(parts.at(2));
  }

  nc_assert_always(!info.isLocked());
  nc_assert_always(info.hasAtomInfo());
  nc_assert_always(info.hasAtomMSD());
  nc_assert_always(info.hasStructureInfo());
  nc_assert_always(!info.hasHKLInfo());
  nc_assert_always(cfg.dcutoff>0.0&&cfg.dcutoff<cfg.dcutoffup);

  const RotMatrix rec_lat = getReciprocalLatticeRot( info );

  const double min_ds_sq(cfg.dcutoff*cfg.dcutoff);
  const double max_ds_sq(cfg.dcutoffup*cfg.dcutoffup);

  //Collect info for each atom in suitable format for use for calculations below:
  SmallVector<SmallVector<Vector,16>,4> atomic_pos;//atomic coordinates
  SmallVectD csl;//coherent scattering length
  SmallVectD msd;//mean squared displacement
  SmallVectD cache_factors;

  AtomList::const_iterator it (info.atomInfoBegin()), itE(info.atomInfoEnd());
  for (;it!=itE;++it) {
    nc_assert( it->msd().has_value() );
    msd.push_back( it->msd().value() );
    csl.push_back( it->atomData().coherentScatLen() );
    SmallVector<Vector,16> pos;
    pos.reserve_hint( it->unitCellPositions().size() );
    for ( const auto& p : it->unitCellPositions() )
      pos.push_back( p.as<Vector>() );
    atomic_pos.push_back( std::move(pos) );
  }

  int max_h, max_k, max_l;
  estimateHKLRange(cfg.dcutoff,rec_lat,max_h, max_k, max_l);

  nc_assert_always(msd.size()==atomic_pos.size());
  nc_assert_always(msd.size()==csl.size());
  cache_factors.resize(csl.size(),0.0);

  //cache some thresholds for efficiency (see below where it is used for more
  //comments):
  SmallVectD whkl_thresholds;
  whkl_thresholds.reserve_hint(csl.size());
  for (size_t i = 0; i<csl.size(); ++i) {
    if ( cfg.fsquarecut < 0.01 && cfg.fsquarecut > fsquarecut_lowest_possible_value )
      whkl_thresholds.push_back(std::log(ncabs(csl.at(i)) / cfg.fsquarecut ) );
    else
      whkl_thresholds.push_back(kInfinity);//use inf when not true that fsqcut^2 << fsq
  }

  //We now conduct a brute-force loop over h,k,l indices, adding calculated info
  //in the following containers along the way:
  HKLList hkllist;
  std::vector<std::vector<short> > eqv_hkl_short;

  //Breaking O(N^2) complexity in compatibility searches by using map (the key
  //is an integer composed from Fsquared and d-spacing, and although clashes are
  //allowed, it should only clash rarely or efficiency is compromised):

#ifdef NCRYSTAL_NCMAT_USE_MEMPOOL
  MemPool pool(10000000);
  MemPoolAllocator<void> poolalloc(&pool);
  FamMap fsq2hklidx(poolalloc);
#else
  FamMap fsq2hklidx;
#endif

  SmallVectD whkl;//outside loop for reusage
  while ( whkl.size() < msd.size() )
    whkl.push_back(1.0);//init with unit factors in case of forceunitdebyewallerfactor

  //NB, for reasons of symmetry we ignore half of the hkl vectors (ignoring
  //h,k,l->-h,-k,-l and 000). This means, half a space, and half a plane and
  //half an axis,  hence the loop limits:

  for( int loop_h=0;loop_h<=max_h;++loop_h ) {
    for( int loop_k=(loop_h?-max_k:0);loop_k<=max_k;++loop_k ) {
      for( int loop_l=-max_l;loop_l<=max_l;++loop_l ) {
        if(loop_h==0 && loop_k==0 && loop_l<=0)
          continue;
        if ( do_select && (loop_h!=select_h||loop_k!=select_k||loop_l!=select_l) )
            continue;
        const Vector hkl(loop_h,loop_k,loop_l);

        //calculate waveVector, wave number and dspacing:
        Vector waveVector = rec_lat*hkl;
        const double ksq = waveVector.mag2();
        const double dspacingsq = (k2Pi*k2Pi)/ksq;
        if( dspacingsq < min_ds_sq || dspacingsq > max_ds_sq )
          continue;

        if (no_forceunitdebyewallerfactor) {
          nclikely fillHKL_getWhkl(whkl, ksq, msd);
        }

        //calculate |F|^2
        double real_or_imag_upper_limit(0.0);
        for( unsigned i=0; i < whkl.size(); ++i ) {
          if ( whkl[i] > whkl_thresholds[i]) {
            cache_factors[i] = 0.0;
            continue;//Abort early to save exp/cos/sin calls. Note that
                     //O(fsquarecut) here corresponds to O(fsquarecut^2)
                     //contributions to final FSquared - for which we demand
                     //>fsquarecut below. We only do this when fsquarecut<1e-2
                     //(see calculations for whkl_thresholds above).
          } else {
            double factor = csl[i]*std::exp(-whkl[i]);
            cache_factors[i] = factor;
            //Assuming cos(phase)=sin(phase)=1 gives us a cheap upper limit on
            //fsquared:
            real_or_imag_upper_limit += atomic_pos[i].size()*factor;
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
        for( unsigned i=0 ; i < whkl.size(); ++i ) {
          double factor = cache_factors[i];
          if (!factor)
            continue;
          auto itAtomPos = atomic_pos[i].begin();
          auto itAtomPosEnd = atomic_pos[i].end();
          StableSum cpsum, spsum;
          for(;itAtomPos!=itAtomPosEnd;++itAtomPos) {
            double phase = hkl.dot(*itAtomPos) * k2Pi;
            double cp,sp;
            sincos(phase,cp,sp);
            cpsum.add(cp);
            spsum.add(sp);
          }
          real.add(cpsum.sum() * factor);
          imag.add(spsum.sum() * factor);
        }
        double realsum = real.sum();
        double imagsum = imag.sum();
        double FSquared = (realsum*realsum+imagsum*imagsum);

        //skip weak or impossible reflections:
        if(FSquared<cfg.fsquarecut)
          continue;

        //normalise waveVector so we can use it below as a demi_normal:
        waveVector *= 1.0 / std::sqrt(ksq);

        const double dspacing = std::sqrt(dspacingsq);//TODO: store dspacingsquared in multimap and avoid some sqrt calls.
        FamKeyType searchkey(keygen(FSquared,dspacing));//key for our fsq2hklidx multimap

        FamMap::iterator itSearchLB = fsq2hklidx.lower_bound(searchkey);
        FamMap::iterator itSearch(itSearchLB), itSearchE(fsq2hklidx.end());
        bool isnewfamily = true;
        for ( ; itSearch!=itSearchE && itSearch->first == searchkey; ++itSearch ) {
          nc_assert(itSearch->second<hkllist.size());
          HKLInfo * hklinfo = &hkllist[itSearch->second];
          if ( ncabs(FSquared-hklinfo->fsquared) < cfg.merge_tolerance*(FSquared+hklinfo->fsquared )
              && ncabs(dspacing-hklinfo->dspacing) < cfg.merge_tolerance*(dspacing+hklinfo->dspacing ) )
            {
              //Compatible with existing family, simply add normals to it.
              hklinfo->demi_normals.push_back(waveVector.as<HKLInfo::Normal>());
              if (cfg.expandhkl) {
                nc_assert(itSearch->second<eqv_hkl_short.size());
                eqv_hkl_short[itSearch->second].push_back(loop_h);
                eqv_hkl_short[itSearch->second].push_back(loop_k);
                eqv_hkl_short[itSearch->second].push_back(loop_l);
              }
              isnewfamily = false;
              break;
            }
        }
        if (isnewfamily) {

          if ( hkllist.size()>1000000 && !env_ignorefsqcut )//guard against crazy setups
            NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach requested dcutoff = "<<cfg.dcutoff<<" Aa");

          HKLInfo hi;
          hi.h=loop_h;
          hi.k=loop_k;
          hi.l=loop_l;
          hi.fsquared = FSquared;
          hi.dspacing = dspacing;
          hi.demi_normals.push_back(waveVector.as<HKLInfo::Normal>());
          fsq2hklidx.insert(itSearchLB,FamMap::value_type(searchkey,hkllist.size()));
          hkllist.emplace_back(std::move(hi));
          if (cfg.expandhkl) {
            eqv_hkl_short.push_back(std::vector<short>());
            std::vector<short>& last = eqv_hkl_short.back();
            last.reserve(3);
            last.push_back(loop_h);
            last.push_back(loop_k);
            last.push_back(loop_l);
          }
        }
      }//loop_l
    }//loop_k
  }//loop_h

  //update HKLlist and copy to info
  info.enableHKLInfo(cfg.dcutoff,cfg.dcutoffup);

  HKLList::iterator itHKL, itHKLB(hkllist.begin()), itHKLE(hkllist.end());
  for(itHKL=itHKLB;itHKL!=itHKLE;++itHKL) {
    unsigned deminorm_size = itHKL->demi_normals.size();
    itHKL->multiplicity=deminorm_size*2;
    if(cfg.expandhkl) {
      std::vector<short>& eh = eqv_hkl_short.at(itHKL-itHKLB);
#if __cplusplus >= 201402L
      //Our make_unique for c++11 seems to have problems with arrays
      itHKL->eqv_hkl = std::make_unique<short[]>(deminorm_size*3);
#else
      itHKL->eqv_hkl = decltype(itHKL->eqv_hkl)(new short[deminorm_size*3]());
#endif
      std::copy(eh.begin(), eh.end(), &itHKL->eqv_hkl[0]);
    }
  }
  info.setHKLList(std::move(hkllist));;

}
