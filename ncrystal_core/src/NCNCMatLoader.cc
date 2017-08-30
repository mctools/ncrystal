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

#include "NCNCMatLoader.hh"
#include "NCVector.hh"
#include "NCrystal/NCInfo.hh"

#include "NCFile.hh"
#include "NCLatticeUtils.hh"
#include "NCrystal/NCException.hh"
#include "NCNeutronSCL.hh"
#include <string>
#include <sstream>
#include <iterator>
#include <limits>
#include <streambuf>
#include <algorithm>
#include <cmath>
#include <iostream>

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


NCrystal::NCMatLoader::NCMatLoader(const char* fn)
  : m_debye_temp(-1.), m_sg(0)
{
  nc_assert(NCrystal::file_exists(fn));
  std::ifstream infile;
  infile.open(fn);
  if (!infile.good())
    NCRYSTAL_THROW2(DataLoadError,"Problems opening file: "<<fn);

  std::string tmpstr;
  getline(infile, tmpstr);
  if (tmpstr !="NCMAT v1")//Todo for NC2: update check here (but still support v1)
    NCRYSTAL_THROW(DataLoadError,"File is in an unsupported NCMAT version" );
  m_version.swap(tmpstr);

  if(ignoreCharNTimes(infile,1,'@').eof())
    NCRYSTAL_THROW2(DataLoadError,"data file does not contain any @" << fn );

  infile.seekg(-1, infile.cur);

  while(infile >> tmpstr)
    {
      if(tmpstr=="@CELL")
        {
          ignoreCharNTimes(infile,1,'\n');

          double a=0., b=0.,c=0.,alpha=0.,beta=0.,gamma=0. ;
          bool hasLen=false; bool hasAng=false;
          std::streampos oldpos;
          while(infile >> tmpstr)
            {
              if (tmpstr=="lengths")
                {
                  infile >> tmpstr;
                  a =  atof(tmpstr.c_str());
                  m_a=a;
                  infile >> tmpstr;
                  b =  atof(tmpstr.c_str());
                  m_b=b;
                  infile >> tmpstr;
                  c =  atof(tmpstr.c_str());
                  m_c=c;
                  hasLen=true;
                }
              else if (tmpstr=="angles")
                {
                  infile >> tmpstr;
                  alpha =  atof(tmpstr.c_str());
                  m_alpha=alpha;
                  infile >> tmpstr;
                  beta =  atof(tmpstr.c_str());
                  m_beta=beta;
                  infile >> tmpstr;
                  gamma =  atof(tmpstr.c_str());
                  m_gamma=gamma;
                  hasAng=true;
                }
              else if(tmpstr[0]=='@')
                {
                  infile.seekg(-tmpstr.size(), infile.cur);
                  break;
                }
            }

          if(!hasLen)
            NCRYSTAL_THROW2(DataLoadError,"NCMatLoader::NCMatLoader data file does not contain unitcell lengths " << fn );

          if(!hasAng)
            NCRYSTAL_THROW2(DataLoadError,"NCMatLoader::NCMatLoader data file does not contain unitcell angles " << fn );

          alpha *= M_PI/180.;
          beta  *= M_PI/180.;
          gamma *= M_PI/180.;

          m_lattice.m_cell = getLatticeRot(a,b,c,alpha,beta,gamma);
          m_lattice.m_reciprocal = getReciprocalLatticeRot(a,b,c,alpha,beta,gamma);
        }
      else if (tmpstr=="@ATOMPOSITIONS")
        {
          m_atomnum=0;
          ignoreCharNTimes(infile,1,'\n');
          while(infile >> tmpstr)
            {
              if (tmpstr[0]!='#' && tmpstr[0]!='@')
                {
                  std::string ele = tmpstr;
                  infile >> tmpstr;
                  double px =  atof(tmpstr.c_str());
                  infile >> tmpstr;
                  double py =  atof(tmpstr.c_str());
                  infile >> tmpstr;
                  double pz =  atof(tmpstr.c_str());

                  Vector pos(px,py,pz);
                  //pos.print();
                  std::map<std::string, std::vector<Vector> >::iterator it = m_lattice.m_atomic_pos.find(ele);
                  if ( it == m_lattice.m_atomic_pos.end() )
                    {
                      std::vector<Vector> vv;
                      vv.push_back(pos);
                      m_lattice.m_atomic_pos[ele]=vv;
                    }
                  else
                    {
                      it->second.push_back(pos);
                    }
                  m_atomnum++;
                  ignoreCharNTimes(infile,1,'\n');

                }
              if(tmpstr[0]=='@')
                {
                  infile.seekg(-tmpstr.size(), infile.cur);
                  break;
                }
            }
        }
      else if (tmpstr=="@SPACEGROUP")
        {
          ignoreCharNTimes(infile,1,'\n');
          while(infile >> tmpstr)
            {
              if (tmpstr[0]!='#' && tmpstr[0]!='@')
                {
                  m_sg = atoi(tmpstr.c_str());
                }
              if(tmpstr[0]=='@')
                {
                  infile.seekg(-tmpstr.size(), infile.cur);
                  break;
                }
            }
        }
      else if (tmpstr=="@DEBYETEMPERATURE")
        {
          ignoreCharNTimes(infile,1,'\n');
          while(infile >> tmpstr)
            {
              if (tmpstr[0]!='#' && tmpstr[0]!='@')
                {
                  if (std::isdigit (tmpstr[0]) )
                    m_debye_temp = atof(tmpstr.c_str());
                  else
                    {
                      std::string element = tmpstr;
                      infile >> tmpstr ;
                      m_debye_map[element] =  atof(tmpstr.c_str());
                    }
                }
              if(tmpstr[0]=='@')
                {
                  infile.seekg(-tmpstr.size(), infile.cur);
                  break;
                }
            }
        }

    }
  infile.close();

  if(m_lattice.m_atomic_pos.size()==0)
    NCRYSTAL_THROW2(DataLoadError,"failed to load lattice info from file " << fn);

  if(!m_sg)
    NCRYSTAL_THROW2(DataLoadError,"unit cell space group is not provided in the file " << fn );

  //TODO for NC2: Debye temperature will become optional when other inelastic models become available
  if(m_debye_temp==-1. && m_debye_map.empty())
    NCRYSTAL_THROW2(DataLoadError,"Debye temperature is not provided in the file " << fn );

  if (!m_debye_map.empty()) {
    std::map<std::string, std::vector<Vector> >::const_iterator itAtom, itAtomE(m_lattice.m_atomic_pos.end());
    for (itAtom=m_lattice.m_atomic_pos.begin();itAtom!=itAtomE;++itAtom) {
      if (m_debye_map.find(itAtom->first)==m_debye_map.end())
        NCRYSTAL_THROW2(DataLoadError,"Some per-elment Debye temperatures are specified but one is missing for \""<<itAtom->first<<"\" in the file: "<<fn);
    }
    if (m_debye_map.size()!=m_lattice.m_atomic_pos.size())
      NCRYSTAL_THROW2(DataLoadError,"Per-elment Debye temperatures for elements not in the lattice are specified in the file: "<<fn);
  }
}

NCrystal::NCMatLoader::~NCMatLoader()
{
}

const NCrystal::NCMatLoader::Lattice& NCrystal::NCMatLoader::getLattice() const
{
  return m_lattice;
}


std::map<std::string, unsigned> NCrystal::NCMatLoader::getAtomMap()
{
  std::map<std::string, unsigned> an;
  for(std::map<std::string, std::vector<Vector> >::const_iterator it=m_lattice.m_atomic_pos.begin();it!=m_lattice.m_atomic_pos.end();++it)
    {
      an[it->first]=it->second.size();
    }
  return an;
}

unsigned NCrystal::NCMatLoader::getSpacegroupNum() const
{
  return m_sg;
}

unsigned NCrystal::NCMatLoader::getNumAtom() const
{
  unsigned tot=0;
  for(std::map<std::string, std::vector<Vector> >::const_iterator it=m_lattice.m_atomic_pos.begin();it!=m_lattice.m_atomic_pos.end();++it)
    {
      tot += it->second.size();
    }
  return tot;
}


std::ifstream& NCrystal::NCMatLoader::ignoreCharNTimes(std::ifstream& file, unsigned num, const char& c)
{
  for(unsigned i=0; i < num ; ++i){
    file.ignore(std::numeric_limits<std::streamsize>::max(),c);
  }
  return file;
}

void NCrystal::NCMatLoader::fillHKL( NCrystal::Info &info,  const std::map<std::string, double>& msdmap,
                                     double min_ds, double max_ds, bool expandhkl ) const
{
  nc_assert_always( min_ds>=0.0 && max_ds > min_ds );
  if (min_ds==0) {
    //very simple heuristics here for now, can likely be improved.
    if (m_lattice.m_atomic_pos.size()>40)
      min_ds = 0.3;
    else
      min_ds = 0.2;
    const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
    if (verbose)
      std::cout<<"NCrystal::NCMATFactory::automatically selected dcutoff level "<< min_ds<< " Aa"<<std::endl;
  }

  nc_assert_always( !msdmap.empty() );

  const double min_ds_sq(min_ds*min_ds);
  const double max_ds_sq(max_ds*max_ds);
  //Reciprocal lattice
  const RotMatrix& rec_lat = m_lattice.m_reciprocal;

  //coherent scattering length
  std::vector<double> csl;//csl per element
  std::vector<double> msd;//msd per element
  std::vector<std::vector<Vector> > atomic_pos;//list of atom positions for each element

  std::map<std::string, std::vector<Vector> >::const_iterator it, itE(m_lattice.m_atomic_pos.end());
  for(it = m_lattice.m_atomic_pos.begin();it != itE; ++it)
    {
      atomic_pos.push_back(it->second);
      csl.push_back(NeutronSCL::instance()->getCoherentSL(it->first));
      std::map<std::string,double>::const_iterator itmsd = msdmap.find(it->first);
      if(itmsd == msdmap.end())
        {
          NCRYSTAL_THROW2(BadInput,"Input mean squared displacement does not contain data for element " << it->first);
        }
      else
        msd.push_back(itmsd->second);
    }

  const double tolerance=1e-6;//tolerance for Fsquare & dspacing comparisons when composing hkl families

  int max_h, max_k, max_l;
  estimateHKLRange(min_ds,rec_lat,max_h, max_k, max_l);

  nc_assert_always(msd.size()==atomic_pos.size());
  nc_assert_always(msd.size()==csl.size());

  //cache some thresholds for efficiency (see below where it is used for more
  //comments):
  std::vector<double> whkl_thresholds;
  whkl_thresholds.reserve(csl.size());
  for (size_t i = 0; i<csl.size(); ++i)
    whkl_thresholds.push_back(log(csl.at(i) / 1e-5 ) );

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
        const Vector waveVector = rec_lat*hkl;
        const double ksq = waveVector.mag2();
        //const double k = sqrt(ksq);
        //const double dspacing = (2*M_PI)/k;
        const double dspacingsq = (4*M_PI*M_PI)/ksq;
        if( dspacingsq<min_ds_sq || dspacingsq>max_ds_sq )
          continue;

        getWhkl(whkl, ksq, msd);
        nc_assert(msd.size()==whkl.size());

        //calculate |F|^2
        double real=0., img=0.;
        for(unsigned i=0;i<whkl.size();i++)
          {
            if ( whkl[i] > whkl_thresholds[i])
              continue;//This is the same as calculating the factor and
                       //requiring factor>1e-5 [and note that O(1e-5) here
                       //corresponds to O(1e-10) contributions to final FSquared
                       //- for which we demand >1e-5 below]. Aborting early
                       //saves expensive exp/cos/sin calls.
            double factor = csl[i]*exp(-whkl[i]);
            std::vector<Vector>::const_iterator itAtomPos(atomic_pos[i].begin()), itAtomPosEnd(atomic_pos[i].end());
            for(;itAtomPos!=itAtomPosEnd;++itAtomPos)
              {
                double phase=hkl.dot(*itAtomPos)*(2.0*M_PI);
                real += cos(phase)*factor;
                img += sin(phase)*factor;
              }
          }
        double FSquared = (real*real+img*img);

        //skip weak or impossible reflections:
        if(FSquared<1e-5)//NB: Hardcoded to same value as in .nxs factory (also note the threshold cut above for early abort)
          continue;

        const double dspacing = sqrt(dspacingsq);
        FamKeyType searchkey(keygen(FSquared,dspacing));//key for our fsq2hklidx multimap

        FamMap::iterator itSearchLB = fsq2hklidx.lower_bound(searchkey);
        FamMap::iterator itSearch(itSearchLB), itSearchE(fsq2hklidx.end());
        bool isnewfamily = true;
        for ( ; itSearch!=itSearchE && itSearch->first == searchkey; ++itSearch ) {
          nc_assert(itSearch->second<hkllist.size());
          HKLInfo * hklinfo = &hkllist[itSearch->second];
          if ( ncabs(FSquared-hklinfo->fsquared) < tolerance*(FSquared+hklinfo->fsquared )
               && ncabs(dspacing-hklinfo->dspacing) < tolerance*(dspacing+hklinfo->dspacing ) )
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
            NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach requested dcutoff = "<<min_ds<<" Aa");

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
  HKLList::iterator itHKL, itHKLB(hkllist.begin()), itHKLE(hkllist.end());
  for(itHKL=itHKLB;itHKL!=itHKLE;++itHKL)
    {
      unsigned deminorm_size = itHKL->demi_normals.size();
      itHKL->multiplicity=deminorm_size*2;
      if(expandhkl)
        {
          std::vector<short>& eh = eqv_hkl_short.at(itHKL-itHKLB);
          itHKL->eqv_hkl = new short[deminorm_size*3];
          std::copy(eh.begin(), eh.end(), itHKL->eqv_hkl);
        }
      info.addHKL(*itHKL);
    }
}

void NCrystal::NCMatLoader::getWhkl(std::vector<double>& out_whkl, const double ksq, const std::vector<double> & msd) const
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
