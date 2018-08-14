////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCFillHKL.hh"
#include "NCOrientUtils.hh"
#include "NCRotMatrix.hh"
#include "NCLatticeUtils.hh"
#include "NCNeutronSCL.hh"
#include "NCrystal/NCDefs.hh"

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

#if __cplusplus >= 201103L && ( defined(__clang__) || !defined(__GNUC__) || __GNUC__ >= 5 )
//Use mempool in C++11 (but not gcc 4.x.y)
#  define NCRYSTAL_NCMAT_USE_MEMPOOL
#endif

#ifdef NCRYSTAL_NCMAT_USE_MEMPOOL
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <scoped_allocator>
namespace NCrystal {
  //We use a simple expanding-only memory pool for the temporary multimap used
  //to detect hkl families. This results better cache locality and should
  //prevent memory fragmentation. TODO for NC2: Consider moving this pool
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
  void fillHKL_getWhkl(std::vector<double>& out_whkl, const double ksq, const std::vector<double> & msd)
  {
    //Sears, Acta Cryst. (1997). A53, 35-45

    if (out_whkl.size()!=msd.size())
      out_whkl.resize(msd.size());//should only happen first time called

    double kk2 = 0.5*ksq;
    std::vector<double>::const_iterator it, itE(msd.end());
    std::vector<double>::iterator itOut(out_whkl.begin());
    for(it=msd.begin();it!=itE;++it)
      *itOut++ = kk2*(*it);
  }
}

void NCrystal::fillHKL( NCrystal::Info &info,
                        double dcutoff, double dcutoffup, bool expandhkl,
                        double fsquarecut, double merge_tolerance )
{
  nc_assert_always(!info.isLocked());
  nc_assert_always(info.hasAtomInfo());
  nc_assert_always(info.hasAtomPositions());
  nc_assert_always(info.hasAtomMSD());
  nc_assert_always(info.hasStructureInfo());
  nc_assert_always(!info.hasHKLInfo());
  nc_assert_always(dcutoff>0.0&&dcutoff<dcutoffup);

  const RotMatrix rec_lat = getReciprocalLatticeRot( info );

  const double min_ds_sq(dcutoff*dcutoff);
  const double max_ds_sq(dcutoffup*dcutoffup);

  //Collect info for each atom in suitable format for use for calculations below:
  std::vector<std::vector<Vector> > atomic_pos;//wyckoff positions
  std::vector<double> csl;////coherent scattering length
  std::vector<double> msd;//mean squared displacement

  AtomList::const_iterator it (info.atomInfoBegin()), itE(info.atomInfoEnd());
  for (;it!=itE;++it) {
    const std::string& elem_name = NeutronSCL::instance()->getAtomName(it->atomic_number);//TODO for NC2: lookup directly with atomic_number
    msd.push_back(it->mean_square_displacement);
    csl.push_back(NeutronSCL::instance()->getCoherentSL(elem_name));
    std::vector<Vector> pos;
    pos.reserve(it->positions.size());
    for (size_t i = 0; i<it->positions.size();++i) {
      pos.push_back( Vector(it->positions.at(i).x,
                            it->positions.at(i).y,
                            it->positions.at(i).z) );
    }
    atomic_pos.push_back(pos);

  }

  int max_h, max_k, max_l;
  estimateHKLRange(dcutoff,rec_lat,max_h, max_k, max_l);

  nc_assert_always(msd.size()==atomic_pos.size());
  nc_assert_always(msd.size()==csl.size());

  //cache some thresholds for efficiency (see below where it is used for more
  //comments):
  std::vector<double> whkl_thresholds;
  whkl_thresholds.reserve(csl.size());
  for (size_t i = 0; i<csl.size(); ++i) {
    if (fsquarecut<0.01)
      whkl_thresholds.push_back(std::log(ncabs(csl.at(i)) / fsquarecut ) );
    else
      whkl_thresholds.push_back(std::numeric_limits<double>::infinity());//use inf when not true that fsqcut^2 << fsq
  }

  //We now conduct a brute-force loop over h,k,l indices, adding calculating
  //info in the following containers along the way:
  NCrystal::HKLList hkllist;
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

  std::vector<double> whkl;//outside loop for reusage

  //NB, for reasons of symmetry we ignore half of the hkl vectors (ignoring
  //h,k,l->-h,-k,-l and 000). This means, half a space, and half a plane and
  //half an axis,  hence the loop limits:

  for( int loop_h=0;loop_h<=max_h;++loop_h ) {
    for( int loop_k=(loop_h?-max_k:0);loop_k<=max_k;++loop_k ) {
      for( int loop_l=-max_l;loop_l<=max_l;++loop_l ) {
        if(loop_h==0 && loop_k==0 && loop_l<=0)
          continue;
        const Vector hkl(loop_h,loop_k,loop_l);

        //calculate waveVector, wave number and dspacing:
        Vector waveVector = rec_lat*hkl;
        const double ksq = waveVector.mag2();
        const double dspacingsq = (4*M_PI*M_PI)/ksq;
        if( dspacingsq<min_ds_sq || dspacingsq>max_ds_sq )
          continue;

        fillHKL_getWhkl(whkl, ksq, msd);
        nc_assert(msd.size()==whkl.size());

        //normalise waveVector so we can use it below as a demi_normal:
        waveVector *= 1.0/std::sqrt(ksq);

        //calculate |F|^2
        double real=0., img=0.;
        for(unsigned i=0;i<whkl.size();i++) {
          if ( whkl[i] > whkl_thresholds[i])
            continue;//Abort early to save exp/cos/sin calls. Note that
                     //O(fsquarecut) here corresponds to O(fsquarecut^2)
                     //contributions to final FSquared - for which we demand
                     //>fsquarecut below. We only do this when fsquarecut<1e-2
                     //(see calculations for whkl_thresholds above).
          double factor = csl[i]*std::exp(-whkl[i]);
          std::vector<Vector>::const_iterator itAtomPos(atomic_pos[i].begin()), itAtomPosEnd(atomic_pos[i].end());
          for(;itAtomPos!=itAtomPosEnd;++itAtomPos) {
            double phase=hkl.dot(*itAtomPos)*(2.0*M_PI);
            real += std::cos(phase)*factor;
            img += std::sin(phase)*factor;
          }
        }
        double FSquared = (real*real+img*img);

        //skip weak or impossible reflections:
        if(FSquared<fsquarecut)
          continue;

        const double dspacing = sqrt(dspacingsq);
        FamKeyType searchkey(keygen(FSquared,dspacing));//key for our fsq2hklidx multimap

        FamMap::iterator itSearchLB = fsq2hklidx.lower_bound(searchkey);
        FamMap::iterator itSearch(itSearchLB), itSearchE(fsq2hklidx.end());
        bool isnewfamily = true;
        for ( ; itSearch!=itSearchE && itSearch->first == searchkey; ++itSearch ) {
          nc_assert(itSearch->second<hkllist.size());
          HKLInfo * hklinfo = &hkllist[itSearch->second];
          if ( ncabs(FSquared-hklinfo->fsquared) < merge_tolerance*(FSquared+hklinfo->fsquared )
              && ncabs(dspacing-hklinfo->dspacing) < merge_tolerance*(dspacing+hklinfo->dspacing ) )
            {
              //Compatible with existing family, simply add normals to it.
              hklinfo->demi_normals.push_back(HKLInfo::Normal(waveVector.x(),waveVector.y(),waveVector.z()));
              if (expandhkl) {
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

          if (hkllist.size()>1000000)//guard against crazy setups
            NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach requested dcutoff = "<<dcutoff<<" Aa");

          NCrystal::HKLInfo hi;
          hi.h=loop_h;
          hi.k=loop_k;
          hi.l=loop_l;
          hi.fsquared = FSquared;
          hi.dspacing = dspacing;
          hi.demi_normals.push_back(HKLInfo::Normal(waveVector.x(),waveVector.y(),waveVector.z()));
          fsq2hklidx.insert(itSearchLB,FamMap::value_type(searchkey,hkllist.size()));
          hkllist.push_back(hi);
          if (expandhkl) {
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
  info.enableHKLInfo(dcutoff,dcutoffup);

  HKLList::iterator itHKL, itHKLB(hkllist.begin()), itHKLE(hkllist.end());
  for(itHKL=itHKLB;itHKL!=itHKLE;++itHKL) {
    unsigned deminorm_size = itHKL->demi_normals.size();
    itHKL->multiplicity=deminorm_size*2;
    if(expandhkl) {
      std::vector<short>& eh = eqv_hkl_short.at(itHKL-itHKLB);
      itHKL->eqv_hkl = new short[deminorm_size*3];
      std::copy(eh.begin(), eh.end(), itHKL->eqv_hkl);
    }
    info.addHKL(*itHKL);
  }

}
