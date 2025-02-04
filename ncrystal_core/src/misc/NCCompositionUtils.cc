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

#include "NCrystal/misc/NCCompositionUtils.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iomanip>
#include <sstream>

namespace NC = NCrystal;
namespace NCCU = NCrystal::CompositionUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace CU = CompositionUtils;

  namespace {
    typedef std::vector< std::tuple<unsigned,unsigned,double> > Flat_ZAfrac;

    template<class T>
    void collect_ZAfrac( Flat_ZAfrac& out,
                         const AtomData& data,
                         double weight,
                         const T& natabprov,
                         CU::ForceIsotopesChoice forceiso )
    {
      if (data.isNaturalElement()) {
        nc_assert(data.Z()!=0);
        if (forceiso==CU::ForceIsotopes) {
          for (auto& Afrac : natabprov(data.Z()))
            out.emplace_back(data.Z(),Afrac.first,Afrac.second*weight);
        } else {
          out.emplace_back(data.Z(),0,weight);
        }
      } else if (data.isSingleIsotope()) {
        nc_assert(data.Z()!=0);
        out.emplace_back(data.Z(),data.A(),weight);
      } else {
        nc_assert(data.isComposite());
        unsigned nc = data.nComponents();
        for ( unsigned i = 0; i < nc; ++i ) {
          const auto& comp = data.getComponent(i);
          collect_ZAfrac(out,*comp.data,weight*comp.fraction,natabprov,forceiso);
        }
      }
    }
  }
}

NC::CU::FullBreakdown NC::CU::createFullBreakdown( const Info::Composition& composition,
                                                   const NC::CU::NaturalAbundanceProvider& natabprov_raw,
                                                   NC::CU::ForceIsotopesChoice forceiso )
{
  Flat_ZAfrac zafrac;
  zafrac.reserve( composition.size()*4 );

  auto natabprov = [&natabprov_raw](unsigned Z)
  {
    if (!natabprov_raw)
      NCRYSTAL_THROW2(CalcError,"Could not determine natural abundances for Z="<<Z<<" (no natural abundance source was provided!)");
    auto natab = natabprov_raw(Z);
    if ( natab.empty() )
      NCRYSTAL_THROW2(CalcError,"Could not determine natural abundances for Z="<<Z);
    StableSum sumnatab;
    for (auto& Afrac : natab)
      sumnatab.add(Afrac.second);
    if ( std::abs(sumnatab.sum()-1.0)>1e-5 )
      NCRYSTAL_THROW2(CalcError,"Invalid (does not add up to 1) natural abundances for Z="<<Z);
    double corrfact = 1.0/sumnatab.sum();
    for (auto& Afrac : natab)
      Afrac.second *= corrfact;
    return natab;
  };

  for ( const auto& ce : composition) {
    nc_assert( ce.fraction > 0.0 && ce.fraction<= 1.0 );
    collect_ZAfrac( zafrac, ce.atom.data(), ce.fraction, natabprov, forceiso );
  }

  //Ok, we now have everything in a flat structure of (Z,A,fraction) values. We
  //now sort the zafrac vector. We don't need stable sort since identical entries
  //are... identical.
  std::sort(zafrac.begin(),zafrac.end());

  //Do one pass where we expand any natural elements that must be expanded (not
  //needed if we already forced iso expansion in collect_ZAfrac):

  auto getZ = [&zafrac] ( unsigned i ) { nc_assert(i<zafrac.size()); return std::get<0>(zafrac[i]); };
  auto getA = [&zafrac] ( unsigned i ) { nc_assert(i<zafrac.size()); return std::get<1>(zafrac[i]); };
  auto getFrac = [&zafrac] ( unsigned i ) ->double& { nc_assert(i<zafrac.size()); return std::get<2>(zafrac[i]); };

  if ( forceiso != ForceIsotopes ) {
    //Bla loop via unsigned idx  and append any expanded isotopes at the end (while setting frac=0 for the A=0 entries that were expanded).
    nc_assert_always(  (uint64_t)zafrac.size() < (uint64_t)(std::numeric_limits<unsigned>::max()) );
    const unsigned nnzafrac = zafrac.size();
    for (unsigned i = 0; i < nnzafrac; ) {
      unsigned zval = getZ(i);
      bool hasNonNat(false);
      StableSum fractNat;
      unsigned inext = i;
      for ( ; (inext < nnzafrac && getZ(inext) == zval ); ++inext ) {
        if ( getA(inext) == 0 ) {
          fractNat.add( getFrac(inext) );
          getFrac(inext) = 0.0;
        } else {
          hasNonNat = true;
        }
      }
      if ( fractNat.sum() > 0.0 ) {
        if ( !hasNonNat ) {
          //oups, we nulled out the fractions of all (zval,A=0) entries when we
          //shouldn't have. Restore.
          assert( getA(i) == 0 && getFrac(i) == 0.0 );
          getFrac(i) = fractNat.sum();
        } else {
          //We must split a natural element into isotopes and append (with
          //weight fractNat):
          for (auto& Afrac : natabprov(zval))
            zafrac.emplace_back(zval,Afrac.first,Afrac.second*fractNat.sum());
        }
      }
      i = inext;
    }
    //Resort:
    std::sort(zafrac.begin(),zafrac.end());
  }

  //Now, loop through zafrac, combining adjacant isotopes with same Z,A values,
  //and filling FullBreakdown result.. Remember to ignore entries with
  //fraction=0. At this point, all Z entries must be either expanded in
  //isotopes, or consist of just natural elements.

  FullBreakdown result;

  nc_assert_always(  (uint64_t)zafrac.size() < (uint64_t)(std::numeric_limits<unsigned>::max()) );
  const unsigned nzafrac = zafrac.size();

  for (unsigned i = 0; i < nzafrac; ) {
    unsigned zval = getZ(i);
    std::vector<std::pair<unsigned,StableSum>> current;
    unsigned inext = i;
    for ( ; (inext < nzafrac && getZ(inext) == zval ); ++inext ) {
      if ( getFrac(inext) == 0.0 )
        continue;
      if ( current.empty() || current.back().first != getA(inext) ) {
        StableSum newsum;
        newsum.add(getFrac(inext));//Todo: StableSum constructor with initial value parameter?
        current.emplace_back( getA(inext), newsum );
      } else {
        current.back().second.add( getFrac(inext) );
      }
    }
    i = inext;

    std::vector<std::pair<unsigned,double>> current_dbls;
    current_dbls.reserve(current.size());
    for (const auto& e : current)
      current_dbls.emplace_back(e.first,e.second.sum());
    result.emplace_back( zval, std::move(current_dbls) );
  }
  return result;
}

std::string NC::CU::ElementBreakdownLW::description(unsigned precision) const
{
  const unsigned zval = Z();
  const unsigned N = nIsotopes();
  auto elemname = elementZToName(zval);
  if ( N == 0 )
    return elemname;//natural element
  std::ostringstream ss;
  ss<<elemname;
  if ( N == 1) {
    ss << firstA();
    return ss.str();
  }
  //Mixture of isotopes:
  ss<<std::setprecision(precision)<<"{";
  for ( unsigned i = 0; i < N; ++i ) {
    ss << fraction(i)<<"*"<<elemname<<A(i);
    if (i+1!=N)
      ss << "+";
  }
  ss<<"}";
  return ss.str();
}

std::string NC::CU::breakdownToStr( const NC::CU::LWBreakdown& bd, unsigned precision )
{
  if ( bd.size()==1 )
    return bd.front().second.description(precision);
  std::ostringstream ss;
  ss<<std::setprecision(precision)<<"Mix{";
  auto Nm1 = bd.size()-1;
  for ( auto&& frac_elem : enumerate(bd) ) {
    ss << frac_elem.val.first<<"*"<<frac_elem.val.second.description(precision);
    if ( frac_elem.idx != Nm1 )
      ss << "+";
  }
  ss<<"}";
  return ss.str();
}


bool NC::CU::ElementBreakdownLW::cmpOthers(const NC::CU::ElementBreakdownLW& o) const
{
  nc_assert( nIsotopes() == o.nIsotopes() );
  nc_assert( !m_other == !o.m_other );
  if (!m_other)
    return false;//neither have composition => equal
  const unsigned ncm1 = nIsotopes()-1;
  for (unsigned i = 0; i < ncm1; ++i) {
    if ( m_other[i] != o.m_other[i] )
      return m_other[i] < o.m_other[i];
  }
  return false;//equal
}

double NC::CU::ElementBreakdownLW::calcFirstFraction() const
{
  const unsigned nc = nIsotopes();
  nc_assert( !!m_other && nc > 1 );
  StableSum tot;
  const unsigned ncm1 = nc - 1;
  for (unsigned i = 0; i < ncm1; ++i)
    tot.add( m_other[i].first );
  return 1.0-tot.sum();
}

NC::CU::ElementBreakdownLW::ElementBreakdownLW(const NC::CU::FullElementBreakdown& eb)
{
  const unsigned Z = eb.first;
  const auto& Afrac = eb.second;
  unsigned N = Afrac.size();
  unsigned Afirst = N ? Afrac.front().first : 0;
  if ( N <= 1 ) {
    //Natural element or single isotope.
    if ( Afrac.empty() || Afrac.front().first == 0 ) {
      Afirst = N = 0;//natural elemnt
    } else {
      //Single isotope.
      N = 1;
    }
  } else {
#if nc_cplusplus >= 201402L
    //Our make_unique for c++11 seems to have problems with arrays
    m_other = std::make_unique<std::pair<double,uint16_t>[]>(N-1);
#else
    m_other = decltype(m_other)(new std::pair<double,uint16_t>[N-1]());
#endif
    StableSum totfrac;
    for (auto af: Afrac)
      totfrac.add( af.second );
    nc_assert( totfrac.sum() > 0.0 );
    double fracnorm = 1.0/totfrac.sum();
    for (unsigned i = 1; i < N; ++i)
      m_other[i-1] = { fracnorm * Afrac.at(i).second, Afrac.at(i).first };
  }
  //encode:
  nc_assert_always( Z < 256 && Afirst < 1024 && N < 16384 && Z > 0);
  m_ZAN = (Z<<24) + (Afirst<<14) + N;
  nc_assert(this->Z()==Z);
  nc_assert(this->firstA()==Afirst);
  nc_assert(N==0||this->firstA()==this->A(0));
  nc_assert(this->nIsotopes()==N);
  nc_assert(bool(m_other)==bool(N>1));
  nc_assert(valid());
}

NC::CU::LWBreakdown NC::CU::createLWBreakdown( const Info::Composition& a,
                                               const NaturalAbundanceProvider& b,
                                               ForceIsotopesChoice c )
{
  //NB: Consider if we (for caching-reasons) want to rounding all final
  //fractions (esp in ElementBreakdownLW constructor)? That could reduce chances
  //that tiny numerical issues in fraction calculations could spoil caching. It
  //is, however, likely to be a rather rare issue...
  FullBreakdown bd = createFullBreakdown(a,b,c);
  LWBreakdown lwbd;
  lwbd.reserve(bd.size());
  for (const auto& e : bd) {
    StableSum totfrac;
    for (auto af: e.second)
      totfrac.add( af.second );
    lwbd.emplace_back(totfrac.sum(),ElementBreakdownLW(e));
  }
  return lwbd;
}

std::string NC::CU::fullBreakdownToJSON( const FullBreakdown& bd )
{
  //encode in json as list of elements of type "[ Z, [ [A_1,fraction_1],...,[A_N,fraction_N] ] ]
  std::ostringstream ss;
  ss << '[';

  for ( auto e : enumerate(bd) ) {
    //e.val type =  std::pair<unsigned,std::vector<std::pair<unsigned,double>>> FullElementBreakdown;//( Z, [(A,fraction),...] )
    ss << '[';
    streamJSON(ss,e.val.first);//Z
    ss <<",[";
    for ( auto eA : enumerate(e.val.second) ) {
      ss<<'[';
      streamJSON(ss,eA.val.first);//A
      ss<<',';
      streamJSON(ss,eA.val.second);//fraction
      ss<<']';
      if ( eA.idx + 1 != e.val.second.size() )
        ss << ',';
    }
    ss << "]]";
    if ( e.idx + 1 != bd.size() )
      ss << ',';
  }

  ss << ']';
  return ss.str();
}
