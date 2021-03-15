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

#include "NCrystal/NCLegacyProcs.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCRandUtils.hh"

namespace NC = NCrystal;
namespace NCL = NCrystal::Legacy;

NCL::RandomBase::~RandomBase() = default;

namespace NCrystal {
  namespace Legacy {
    namespace {

      class RandomLegacyWrapper final : public RandomBase {
        shared_obj<RNG> m_rng;
      public:
        RandomLegacyWrapper( shared_obj<RNG> rng ) : m_rng(std::move(rng)) {}
        double generate() override { return m_rng->generate(); }
        RNG& wrapped() { return *m_rng; }
        const RNG& wrapped() const { return *m_rng; }
        shared_obj<const RNG> wrapped_so() { return m_rng; }
      private:
      };

    }
  }
}

NCL::CalcBase::CalcBase(const char * calculator_type_name)
  : m_name(calculator_type_name)
{
}

NCL::CalcBase::~CalcBase()
{
  for (unsigned i=0;i<m_subcalcs.size();++i)
    m_subcalcs[i]->unref();
}

void NCL::CalcBase::registerSubCalc(CalcBase*sc)
{
  if (sc) {
    sc->ref();
    m_subcalcs.push_back(sc);
  }
}

bool NCL::CalcBase::isSubCalc(const CalcBase* cb) const
{
  for (unsigned i=0;i<m_subcalcs.size();++i)
    if (cb==m_subcalcs[i])
      return true;
  return false;
}

void NCL::CalcBase::initLegacyRng() const
{
  nc_assert(!m_legacyrng);
  if (!m_rng)
    getModernRNG();//trigger setting m_rng
  nc_assert(!!m_rng);
  m_legacyrng = makeRC<RandomLegacyWrapper>(m_rng);
}

NCL::Process::Process(const char * calculator_type_name)
  : CalcBase(calculator_type_name)
{
}

NCL::Process::~Process()
{
}

double NCL::Process::crossSectionNonOriented(double ekin ) const
{
  if (isOriented())
    NCRYSTAL_THROW(BadInput,"Process::crossSectionNonOriented called for oriented object.");
  double indir[3] = { 0., 0., 1. };
  return crossSection(ekin, indir);
}

void NCL::Process::validate()
{
  double test_dir[3] = { 0., 0., 1. };
  double ekin_low;
  double ekin_high;
  double xs_low(0.0),xs_high(0.0);

  domain(ekin_low, ekin_high);
  if ( ekin_low!=ekin_low || ekin_high!=ekin_high || ekin_low<0.0 || !(ekin_low<ekin_high) ) {
    //Invalid domain, unless low=high=inf, which is a special case used to
    //indicate a process with zero cross-section everywhere.
    if ( ! (ncisinf(ekin_low) && ncisinf(ekin_high) ) )
      NCRYSTAL_THROW2(LogicError,getCalcName()<<" returns invalid domain!");
  }
  bool oriented(isOriented());
  double test_ekin((1.0-1.0e-9)*ekin_low);
  if (test_ekin < ekin_low)
    xs_low = oriented ? crossSection(test_ekin, test_dir) : crossSectionNonOriented(test_ekin);
  test_ekin = (1.0+1.0e-9)*ekin_high;
  if (test_ekin > ekin_high)
    xs_high = oriented ? crossSection(test_ekin, test_dir) : crossSectionNonOriented(test_ekin);

  if ( xs_low!=xs_low || xs_high!=xs_high || xs_low!=0.0 || xs_high!=0.0 )
    NCRYSTAL_THROW2(LogicError,getCalcName()<<" returns invalid cross sections outside domain!")
}

bool NCL::Process::isNull() const
{
  double ekin_low, ekin_high;
  domain(ekin_low, ekin_high);
  return ncisinf(ekin_low) || !(ekin_low<ekin_high);
}

NCL::Absorption::Absorption(const char* calc_type_name)
  : Process(calc_type_name)
{
}

NCL::Absorption::~Absorption()
{
}

NCL::Scatter::Scatter(const char * calculator_type_name)
  : Process(calculator_type_name)
{
}

NCL::Scatter::~Scatter()
{
}

void NCL::Scatter::generateScatteringNonOriented( double ekin, double& angle, double& de ) const
{
  if (isOriented())
    NCRYSTAL_THROW(BadInput,"Scatter::generateScatteringNonOriented called for oriented object.");
  Vector indir = { 0., 0., 1. };
  Vector outdir;
  generateScattering( ekin, indir.rawArray(), outdir.rawArray(), de );
  angle = indir.angle( outdir );
}

NCL::ScatterComp::ScatterComp(const char * calculator_type_name)
  : Scatter(calculator_type_name), m_threshold_lower(0.0), m_threshold_upper(kInfinity), m_isOriented(-1)
{
}

NCL::ScatterComp::~ScatterComp()
{
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  for (;it!=itE;++it)
    it->scatter->unref();
}

bool NCL::ScatterComp::Component::operator<(const NCL::ScatterComp::Component& o) const
{
  return o.threshold_lower > threshold_lower;
}

void NCL::ScatterComp::addComponent(Scatter* scat, double thescale )
{
  RCGuard guard(scat);//ensure we always ref/unref scat even in case of exceptions.
  if (!scat)
    NCRYSTAL_THROW(BadInput,"ScatterComp::addComponent Got NULL scatter.");
  if (thescale<0.0)
    NCRYSTAL_THROW(BadInput,"ScatterComp::addComponent Component scale is negative.");
  std::vector<Component>::const_iterator it(m_calcs.begin()), itE(m_calcs.end());
  for (;it!=itE;++it) {
    if (it->scatter == scat)
        NCRYSTAL_THROW(BadInput,"ScatterComp::addComponent got same scatter multiple times.");
  }
  m_calcs.reserve(m_calcs.size()+1);
  scat->validate();
  Component c;
  c.scale = thescale;
  c.scatter = scat;
  scat->domain(c.threshold_lower,c.threshold_upper);
  if (m_calcs.empty() || c.threshold_lower < m_threshold_lower)
    m_threshold_lower = c.threshold_lower;
  if (m_calcs.empty() || c.threshold_upper > m_threshold_upper)
    m_threshold_upper = c.threshold_upper;
  registerSubCalc(scat);
  m_calcs.push_back(c);
  scat->ref();

  //Sort by threshold so lowers comes first, preserving the original order when
  //thresholds are equal:
  std::stable_sort(m_calcs.begin(),m_calcs.end());

  m_isOriented = -1;//invalidate (scatter might be incomplete ScatterComp so we
                    //can't always know already if it is oriented or not).

  validate();
}

double NCL::ScatterComp::crossSection(double ekin, const double (&indir)[3] ) const
{
  double c(0);
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  if (it==itE)
    NCRYSTAL_THROW(BadInput,"ScatterComp::crossSection queried with no components added.");
  for (;it!=itE;++it) {
    if (ekin<it->threshold_lower)
      break;
    if (ekin>it->threshold_upper)
      continue;
    c += it->scatter->crossSection(ekin,indir) * it->scale;
  }
  return c;
}

void NCL::ScatterComp::generateScattering( double ekin, const double (&indir)[3],
                                                double (&outdir)[3], double& de ) const
{
  double rand_choice = getRNG()->generate() * crossSection(ekin,indir);
  double c(0);
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  if (it==itE)
    NCRYSTAL_THROW(BadInput,"ScatterComp::generateScattering queried with no components added.");
  for (;it!=itE;++it) {
    if (ekin<it->threshold_lower)
      break;
    if (ekin>it->threshold_upper)
      continue;
    c += it->scatter->crossSection(ekin,indir) * it->scale;
    if (rand_choice <= c) {
      it->scatter->generateScattering(ekin, indir, outdir, de);
      return;
    }
  }
  //Should get here only in case of rounding errors or if called outside
  //domain(). No cross-section means no action:
  outdir[0] = indir[0];
  outdir[1] = indir[1];
  outdir[2] = indir[2];
  de = 0;
}

bool NCL::ScatterComp::isOriented() const {
  if (m_isOriented==-1)
    checkIsOriented();
  return (bool)m_isOriented;
}

void NCL::ScatterComp::checkIsOriented() const
{
  m_isOriented = 0;
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  for (;it!=itE;++it) {
    if (it->scatter->isOriented()) {
      m_isOriented = 1;
      break;
    }
  }
}

namespace NCrystal {
  namespace Legacy {
    //Wrap new-style classes into the old-style infrastructure:
    class ScatterImplWrapper final : public Scatter {
    public:
      ScatterImplWrapper( ProcImpl::ProcPtr proc )
        : Scatter(proc->name()),
          m_proc(std::move(proc))
      {
        nc_assert_always( m_proc != nullptr );
        nc_assert_always( m_proc->processType() == ProcessType::Scatter );
        auto d = m_proc->domain();
        d.elow.get();
        d.ehigh.get();
        if ( m_proc->isNull() )
          m_domain = PairDD(kInfinity,kInfinity);
        else
          m_domain = PairDD(d.elow.get(),d.ehigh.get());
        validate();
      }

      double crossSection(double ekin, const double (&dir)[3] ) const final
      {
        return m_proc->crossSection(m_cacheptr,NeutronEnergy{ekin},NeutronDirection{dir}).get();
      }
      double crossSectionNonOriented( double ekin ) const final
      {
        return m_proc->crossSectionIsotropic(m_cacheptr,NeutronEnergy{ekin}).get();
      }
      void domain(double& ekin_low, double& ekin_high) const final
      {
        ekin_low = m_domain.first;
        ekin_high = m_domain.second;
      }
      bool isOriented() const final { return m_proc->materialType() != MaterialType::Isotropic; }
      void generateScattering( double ekin, const double (&indir)[3],
                               double (&outdir)[3], double& delta_ekin ) const final
      {
        auto outcome = m_proc->sampleScatter( m_cacheptr, getModernRNG(), NeutronEnergy{ekin}, NeutronDirection{indir} );
        delta_ekin = outcome.ekin.get() - ekin;
        outdir[0] = outcome.direction[0];
        outdir[1] = outcome.direction[1];
        outdir[2] = outcome.direction[2];
      }

      void generateScatteringNonOriented( double ekin, double& angle, double& delta_ekin ) const final
      {
        auto outcome = m_proc->sampleScatterIsotropic( m_cacheptr, getModernRNG(), NeutronEnergy{ekin} );
        delta_ekin = outcome.ekin.get() - ekin;
        nc_assert_always( outcome.mu.get()>=-1.0 && outcome.mu.get()<=1.0 );
        angle = std::acos( outcome.mu.get() );
      }

    protected:
      virtual ~ScatterImplWrapper() = default;
    private:
      ProcImpl::ProcPtr m_proc;
      PairDD m_domain;
      mutable CachePtr m_cacheptr;//NB: not MT safe
    };

    class AbsorptionImplWrapper final : public Absorption {
    public:
      AbsorptionImplWrapper( ProcImpl::ProcPtr proc )
        : Absorption(proc->name()),
          m_proc(std::move(proc))
      {
        nc_assert_always( m_proc != nullptr );
        nc_assert_always( m_proc->materialType() == MaterialType::Isotropic );
        nc_assert_always( m_proc->processType() == ProcessType::Absorption );
        auto d = m_proc->domain();
        d.elow.get();
        d.ehigh.get();
        if ( m_proc->isNull() )
          m_domain = PairDD(kInfinity,kInfinity);
        else
          m_domain = PairDD(d.elow.get(),d.ehigh.get());
        validate();
      }

      double crossSection(double ekin, const double (&)[3] ) const final
      {
        return m_proc->crossSectionIsotropic(m_cacheptr,NeutronEnergy{ekin}).get();
      }
      double crossSectionNonOriented( double ekin ) const final
      {
        return m_proc->crossSectionIsotropic(m_cacheptr,NeutronEnergy{ekin}).get();
      }
      void domain(double& ekin_low, double& ekin_high) const final
      {
        ekin_low = m_domain.first;
        ekin_high = m_domain.second;
      }
      bool isOriented() const final { return false; }

    protected:
      virtual ~AbsorptionImplWrapper() = default;
    private:
      ProcImpl::ProcPtr m_proc;
      PairDD m_domain;
      mutable CachePtr m_cacheptr;//NB: not MT safe
    };
  }
}

NCL::RCHolder<const NCL::Scatter> NCL::wrapModernProcPtrInLegacyScatterClass( NC::ProcImpl::ProcPtr p )
{
  nc_assert_always( p->processType() == ProcessType::Scatter );
  return makeRC<ScatterImplWrapper>( std::move(p) );
}

NCL::RCHolder<const NCL::Absorption> NCL::wrapModernProcPtrInLegacyAbsorptionClass( NC::ProcImpl::ProcPtr p )
{
  nc_assert_always( p->processType() == ProcessType::Absorption );
  return makeRC<AbsorptionImplWrapper>( std::move(p) );
}
