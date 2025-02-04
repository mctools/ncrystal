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

#include "NCrystal/interfaces/NCAtomData.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include <sstream>
#include <type_traits>

namespace NC = NCrystal;

struct NC::AtomData::Impl {
  //Ensure Components doesn't start to throw during operations below:
  static_assert(std::is_nothrow_destructible<Component>::value,"");
  static_assert(std::is_nothrow_copy_constructible<Component>::value,"");
  static_assert(std::is_nothrow_move_constructible<Component>::value,"");

  static void clearComponents( AtomData* THIS )
  {
    Component* it = THIS->m_components;
    nc_assert( it != nullptr );
    for ( auto i : ncrange( THIS->nComponents() ) ) {
      (void)i;
      (it++)->~Component();
    }
    if ( THIS->m_components ) {
      AlignedAlloc::freeAlignedAlloc<Component>( THIS->m_components );
      THIS->m_components = nullptr;
    }
  }
  static void setComponents(AtomData* THIS, const Component* o_begin, unsigned n )
  {
    if ( THIS->m_components )
      clearComponents(THIS);
    if ( n == 0)
      return;
    nc_assert_always( n < static_cast<unsigned>(-std::numeric_limits<decltype(m_classify)>::lowest()) );
    Component * it = AlignedAlloc::alignedAlloc<Component>(n);
    //RAII, start at m_classify=0 so destructor won't clean up something we didn't construct.
    THIS->m_components = it;
    THIS->m_classify = 0;
    for ( auto i : ncrange( n ) ) {
      (void)i;
      new ( it++ ) Component( *(o_begin++) );//copy construct in place
      --(THIS->m_classify);//one more ready to destruct
    }
    nc_assert_always( THIS->nComponents() == n );
  }

  static void setComponents(AtomData* THIS,const AtomData::ComponentList& cl)
  {
    setComponents( THIS, &cl[0], cl.size() );
  }
  static void setComponents(AtomData* THIS, const AtomData& o)
  {
    setComponents( THIS, o.m_components, o.nComponents() );
  }
};

NC::AtomData::~AtomData()
{
  assert( (this->m_components != nullptr) == ( this->nComponents() > 0 ) );
  if ( m_components )
    Impl::clearComponents(this);
}

std::string NC::AtomData::elementName() const
{
  nc_assert(isElement());
  unsigned z = Z();
  std::string name = elementZToName(z);
  if (name.empty())
    NCRYSTAL_THROW2(BadInput,"Z-value ("<<z<<") of element is out of range");
  return name;
}

NC::AtomData::AtomData( SigmaBound incXS, double cohSL, SigmaAbsorption captureXS, AtomMass avrMassAMU, unsigned Z, unsigned A )
  : m_m(DoValidate,avrMassAMU), m_ixs(incXS.get()), m_csl(cohSL), m_axs(captureXS.get()), m_classify(A), m_z(Z)
{
  nc_assert(m_z>0);
  nc_assert(m_classify>=0);
  nc_assert(isElement());
  nc_assert(A==0||isSingleIsotope());
  nc_assert(!isComposite());
  nc_assert(m_ixs>=0.0);
  nc_assert(m_axs>=0.0);
  elementName();//call to trigger range check on Z
}

NC::AtomData::AtomData( const NC::AtomData::ComponentList& components )
  : m_m(), m_ixs(0.0), m_csl(0.0), m_axs(0.0), m_classify(0), m_z(0)
{
  nc_assert_always(!components.empty());
  nc_assert_always( static_cast<uint64_t>(components.size())
                    < static_cast<uint64_t>(-std::numeric_limits<decltype(m_classify)>::lowest()) );

  constexpr double fractol = 1e-9;

  if ( components.size()==1 ) {
    //Special case, only one component. We simply become a copy of that
    //component.
    nc_assert_always(ncabs(components.front().fraction-1.0)<fractol);
    const AtomData& c = components.front().data;
    m_m = c.m_m;
    m_ixs = c.m_ixs;
    m_csl = c.m_csl;
    m_axs = c.m_axs;
    m_z = c.m_z;
    m_classify = c.m_classify;
    if ( c.m_components )
      Impl::setComponents( this, c );
  } else {
    //Usual case, constructed from a list of components. We must combine these
    //to calculate new values of constants.
    StableSum sum_w, sum_wm, sum_wixs, sum_wcsl, sum_waxs;

    unsigned commonA = 0;
    unsigned commonZ = 0;
    bool allNaturalElements = true;

    bool first = true;
    for ( const auto& cf : components ) {
      const double w = cf.fraction;
      nc_assert_always(w>0.0&&w<=1.0);
      sum_w.add(w);
      const AtomData& c = cf.data;
      sum_wixs.add(w*c.m_ixs);
      sum_wcsl.add(w*c.m_csl);
      sum_waxs.add(w*c.m_axs);
      sum_wm.add(w*c.m_m.dbl());
      if ( allNaturalElements && !c.isNaturalElement() )
        allNaturalElements = false;
      if (first) {
        first = false;
        if (c.isElement())
          commonZ = c.Z();
        if (c.isSingleIsotope())
          commonA = c.A();
      } else {
        if ( commonZ && ( !c.isElement() || c.Z() != commonZ ) )
          commonZ = 0;
        if ( commonA && ( !c.isSingleIsotope() || c.A() != commonA ) )
          commonA = 0;
      }
    }
    const double sumw_val  = sum_w.sum();
    if (ncabs(sumw_val-1.0)>fractol)
      NCRYSTAL_THROW(BadInput,"Inconsistent atom data - component fractions do not add up to 1.0");

    if ( commonZ && ( commonA || allNaturalElements ) ) {
      //Special case: All components are actually the same isotope or natural element!!!
      //Copy over values from the first entry:
      const AtomData& c = components.front().data;
      nc_assert_always(!c.isComposite());//isotope or nat. elem, so not composite.
      m_m = c.m_m;
      m_ixs = c.m_ixs;
      m_csl = c.m_csl;
      m_axs = c.m_axs;
      m_classify = c.m_classify;
      m_z = c.m_z;
      for ( auto i : ncrange ( 1u, (unsigned)components.size() ) ) {
        const AtomData& cc = components.at(i).data;
        if ( cc.isComposite() || m_m != cc.m_m || m_ixs != cc.m_ixs ||
             m_csl != cc.m_csl || m_axs != cc.m_axs ||
             m_classify != cc.m_classify || m_z != cc.m_z ) {
          NCRYSTAL_THROW(BadInput,"Composite atom data constructed from list of supposedly identical parts -- but some values differ!");
        }
      }

    } else {
      m_z = commonZ;
      //Time to calculate new values.
      const double w_correction_factor = 1.0/sumw_val;//"snap" to unity.
      //Some vars are just straight-forward averages:
      m_csl = sum_wcsl.sum()*w_correction_factor;
      m_axs = sum_waxs.sum()*w_correction_factor;
      m_m.set( sum_wm.sum()*w_correction_factor );
      //But the incoherent cross section is more difficult, as it is a
      //variance. We calculate the combined variance based on the combined mean,
      //hence the need for a second loop:
      StableSum sum_w_newixs;
      for ( const auto& cf : components ) {
        const double w = cf.fraction;
        const AtomData& c = cf.data;
        double tmp = c.m_csl-m_csl;
        sum_w_newixs.add( w * (c.m_ixs+k4Pi*tmp*tmp) );
      }
      m_ixs = sum_w_newixs.sum()*w_correction_factor;
      Impl::setComponents(this,components);
    }
  }
  m_m.validate();
}

std::string NC::AtomData::description(bool includeValues) const
{
  std::ostringstream ss;
  descriptionToStream(ss,includeValues);
  return ss.str();
}

void NC::AtomData::descriptionToStream(std::ostream& os, bool includeValues) const
{
  if ( isNaturalElement() ) {
    os<<elementName();//As guaranteed in NCAtomData.hh!!
  } else if ( isSingleIsotope() ) {
    os<<elementName()<<A();//As guaranteed in NCAtomData.hh!!
  } else {
    nc_assert(isComposite());
    if (isElement())
      os<<elementName();
    else
      os<<"Mix";
    os<<"{";
    unsigned nc = -m_classify;
    for ( unsigned i = 0; i < nc; ++i ) {
      os<< m_components[i].fraction*100.0<<"%";
      m_components[i].data->descriptionToStream(os,false);
      if (i+1!=nc)
        os<<"+";
    }
    os<<"}";
  }
  if (!includeValues)
    return;
  os<<"(cohSL="<<coherentScatLenFM()<<"fm"
    <<" cohXS="<<coherentXS()
    <<" incXS="<<incoherentXS()
    <<" absXS="<<captureXS()
    <<" mass="<<averageMassAMU();
  if (isElement())
    os<<" Z="<<Z();
  if (isSingleIsotope())
    os<<" A="<<A();
  os<<")";
}

bool NC::AtomData::operator<(const NC::AtomData & o) const
{
  //Z (multi-Z comes *after* elements)
  unsigned Zval = isElement() ? Z() : 999999;
  unsigned oZval = o.isElement() ? o.Z() : 999999;
  if ( Zval != oZval )
    return Zval < oZval;
  //cheap equality check before wasting more time below:
  if ( getUniqueID() == o.getUniqueID() )
    return false;
  //A (natural elements comes before isotopes):
  unsigned Aval = isSingleIsotope() ? A() : 0;
  unsigned oAval = o.isSingleIsotope() ? o.A() : 0;
  if ( Aval != oAval )
    return Aval < oAval;
  //Only get here rarely, hopefully.
  std::string descr = description();
  std::string odescr = o.description();
  if ( descr != odescr )
    return descr < odescr;
  return getUniqueID() < o.getUniqueID();
}

bool NC::AtomData::sameValuesAs(const AtomData& o, double rtol, double atol) const
{
  if ( m_classify != o.m_classify
       || m_z != o.m_z
       || !floateq(m_m.dbl(),  o.m_m.dbl(),  rtol,atol)
       || !floateq(m_ixs,o.m_ixs,rtol,atol)
       || !floateq(m_csl,o.m_csl,rtol,atol)
       || !floateq(m_axs,o.m_axs,rtol,atol) )
    return false;
  if (!m_components)
    return true;
  //Ok, overall numbers match but must also check components.
  unsigned nc = -m_classify;
  nc_assert_always(nc>0);
  for (unsigned i = 0; i < nc; ++i) {
    if ( !floateq(m_components[i].fraction,o.m_components[i].fraction,rtol,atol)
         || !m_components[i].data->sameValuesAs(o.m_components[i].data,rtol,atol) )
      return false;
  }
  return true;
}

std::size_t NC::AtomData::hash() const
{
  static_assert(std::is_same<NC::HashValue,std::size_t>::value,"HashValue type changed");
  NC::HashValue tmp = calcHash(m_classify);
  hash_combine(tmp,m_z);
  hash_combine(tmp,m_m.dbl());
  hash_combine(tmp,m_ixs);
  hash_combine(tmp,m_csl);
  hash_combine(tmp,m_axs);
  if (isComposite()) {
    unsigned nc = -m_classify;
    nc_assert_always(nc>0);
    for (unsigned i = 0; i < nc; ++i) {
      hash_combine(tmp,m_components[i].fraction);
      hash_combine(tmp,m_components[i].data->hash());
    }
  }
  return tmp;
}

const std::string& NC::AtomData::elementZToName(unsigned z)
{
  return ::NC::elementZToName(z);
}

unsigned NC::AtomData::elementNameToZ(const std::string& name)
{
  return ::NC::elementNameToZ(name);
}
