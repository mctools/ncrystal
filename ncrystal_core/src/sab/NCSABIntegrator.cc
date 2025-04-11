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

#include "NCrystal/internal/sab/NCSABIntegrator.hh"
#include "NCrystal/internal/sab/NCSABSamplerModels.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

namespace NC = NCrystal;
namespace NS = NCrystal::SAB;

struct NC::SAB::SABIntegrator::Impl : private NoCopyMove {

  Impl( shared_obj<const SABData>,
        const VectD* egrid,
        std::shared_ptr<const SABExtender> );
  void doit(SABXSProvider *, SABSampler*, Optional<std::string>*);
  double determineEMax( const double ) const;
  double determineEMin( const double ) const;
  void setupEnergyGrid();

  //Input data:
  shared_obj<const SABData> m_data;
  VectD m_egrid;
  std::shared_ptr<const SABExtender> m_extender;

  //Data derived from m_data:
  std::shared_ptr<const SABSamplerAtE_Alg1::CommonCache> m_derivedData;

  //Setting;
  SABSampler::EGridMargin m_egridMargin;

  typedef std::unique_ptr<SABSamplerAtE> SamplerAtE_uptr;
  std::pair<SamplerAtE_uptr,double> analyseEnergyPoint(double ekin, bool doSampler ) const;

  double analyseEnergyPointForXS(double ekin) const
  {
    return analyseEnergyPoint( ekin, false ).second;
  }

};

NS::SABIntegrator::~SABIntegrator() = default;

NS::SABIntegrator::SABIntegrator( shared_obj<const SABData> data,
                                  const VectD* egrid,
                                  std::shared_ptr<const SABExtender> sabextender )
  : m_impl(std::move(data),egrid,std::move(sabextender))
{
}

void NS::SABIntegrator::doit(SABXSProvider * out_xs, SABSampler* out_sampler, Optional<std::string>* json)
{
  m_impl->doit(out_xs,out_sampler,json);
}

NS::SABIntegrator::Impl::Impl( shared_obj<const SABData> data,
                               const VectD* egrid,
                               std::shared_ptr<const SABExtender> sabextender )
  : m_data(std::move(data)),
    m_egrid((egrid&&!egrid->empty())?*egrid:VectD()),
    m_extender(!sabextender?std::make_unique<SABFGExtender>(m_data->temperature(),m_data->elementMassAMU(),m_data->boundXS()):std::move(sabextender)),
    m_egridMargin{ 1.05 }
{
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Derived data factory:
    typedef std::pair<UniqueIDValue, shared_obj<const SABData>* > D2DDKey;
    typedef SAB::SABSamplerAtE_Alg1::CommonCache DerivedData;
    class SABData2DerivedDataFactory : public NC::CachedFactoryBase<D2DDKey,DerivedData> {
    public:
      const char* factoryName() const final { return "SABData2DerivedDataFactory"; }
      std::string keyToString( const D2DDKey& key ) const final
      {
        std::ostringstream ss;
        ss<<"(SABData id="<<key.first.value<<")";
        return ss.str();
      }
    protected:
      virtual ShPtr actualCreate( const D2DDKey& key ) const final
      {
        shared_obj<const SABData> data = *key.second;
        nc_assert( !!data && data->getUniqueID()==key.first );

        const auto& alphaGrid = data->alphaGrid();
        const auto& sab = data->sab();
        const std::size_t nalpha = alphaGrid.size();
        const std::size_t nbeta = data->betaGrid().size();
        const std::size_t nalpham1 = nalpha-1;

        //Calculate log(S) values:
        VectD logsab;
        logsab.reserve(sab.size());
        for (auto e: sab)
          logsab.push_back( e>0.0 ? std::log(e) : -kInfinity);

        //For each beta-idx, integrate each grid cell along alpha:
        VectD alphaintegrals_cumul;
        alphaintegrals_cumul.resize(sab.size(),0.);
        std::size_t global_idx(0);
        for (std::size_t ibeta = 0; ibeta<nbeta; ++ibeta, ++global_idx) {
          double cumul = 0.0;
          for (std::size_t ai = 0; ai < nalpham1; ++ai, ++global_idx) {
            auto next_idx = global_idx + 1;
            double integ = SABUtils::integrateAlphaInterval_fast( vectAt( alphaGrid, ai ), vectAt( sab, global_idx ),
                                                                  vectAt( alphaGrid, ai+1 ), vectAt( sab,  next_idx ),
                                                                  vectAt( logsab, global_idx ), vectAt( logsab,  next_idx ) );
            vectAt( alphaintegrals_cumul, next_idx ) = (cumul += integ);
          }
        }
        nc_assert(global_idx==sab.size());

        //Wrap up and return:
        return std::make_shared<const DerivedData>(DerivedData{data,std::move(logsab),std::move(alphaintegrals_cumul)});
      }
    };
    static SABData2DerivedDataFactory s_SABData2DerivedDataFactory;

  }
}

double NS::SABIntegrator::Impl::determineEMin(const double e_upp) const
{
  // We want to find the energy point below which f(E)=sqrt(E)*sigma(E) essentially
  // becomes a constant (within the given tolerance).
  const double tolerance = 1e-3;
  assert(e_upp>1e-99);
  auto f = [this](double e) { return std::sqrt(e)*analyseEnergyPointForXS(e); };
  double e_old( e_upp*0.9 );
  double f_old( f(e_old) );
  while(true) {
    const double e_new = 0.5*e_old;
    if (e_new <= 1e-30*e_upp)
      return ncmin(e_old,e_upp*0.01);
    const double f_new = f(e_new);
    if (f_new==0.0) {
      NCRYSTAL_WARN("Encountered sqrt(E)*sigma(E)=0 at E="<<e_new
                    <<" while searching for suitable Emin value at which"
                    " to start SAB energy grid. Will revert to using"
                    " Emin=0.001*Emax.");
      return 0.001*e_upp;
    }
    if ( ncabs(f_old/f_new-1.0) < tolerance )
      return e_old;
    e_old = e_new;
    f_old = f_new;
  }
}

double NS::SABIntegrator::Impl::determineEMax(const double eUpperLimit) const
{
  FreeGasXSProvider freeXS( m_data->temperature(),
                            m_data->elementMassAMU(),
                            m_data->boundXS() );
  const double stepsize = 0.95;//NB: changing this to e.g. 0.999 might actually give less precise results! Not so simple to just change it.
  auto xs_dist_to_freeXS = [this,&freeXS](double ekin)
                           {
                             return ncabs( this->analyseEnergyPointForXS(ekin) - freeXS.crossSection(NeutronEnergy{ekin}).dbl() );
                           };
  double prevDist = kInfinity;
  double ekin = eUpperLimit;
  const double elow = eUpperLimit * 1e-4;//1e-4 and stepsize 0.95 ensures this
                                         //fct will have at worst 180 calls to
                                         //analyseEnergyPoint (but usually way
                                         //less).
  while (ekin > elow) {
    double dist = xs_dist_to_freeXS(ekin);
    if ( dist > prevDist ) {
      return stepsize*ekin;//extra stepsize here is safety.
    }
    prevDist = dist;
    ekin *= stepsize;
  }
  return 0.0;
}

void NS::SABIntegrator::Impl::setupEnergyGrid()
{
  //Three numbers: emin (-1 means auto), emax (-1 means auto), npts (-1 means auto).

  if (m_egrid.size()>3) {
    //assume user already passed in a correct grid
  } else {
    if ( !m_egrid.empty() && m_egrid.size()!=3 )
      NCRYSTAL_THROW(BadInput,"SABIntegrator invalid energy grid. It must either be a complete array, empty, or consist of three numbers: {emin, emax, npts}");
    double emin, emax;
    unsigned npts;
    if (m_egrid.size()==3) {
      emin = m_egrid.at(0);
      emax = m_egrid.at(1);
      double npts_fp = m_egrid.at(2);
      npts = static_cast<unsigned>(npts_fp);
      if ( npts != npts_fp )
        NCRYSTAL_THROW(BadInput,"SABIntegrator invalid energy grid. When the array has 3 elements, the third must be an integral number representing number of points.");
    } else {
      emin = 0.0;
      emax = 0.0;
      npts = 0;
    }

    if ( !(emin>=0.0) || !(emax>=0.0) || !( emax==0.0 || emin==0.0 || emax>emin ) )
      NCRYSTAL_THROW(BadInput,"SABIntegrator invalid energy grid. Values for emin/emax must fullfil 0<emin<emax or be 0 indicating automatic determination.");

    if ( npts==0 )
      npts = 300;

    const double kT = m_data->temperature().kT();

    if ( emax == 0.0 && m_data->suggestedEmax() > 0.0 ) {
      if ( emin != 0.0 && m_data->suggestedEmax() <= emin )
        NCRYSTAL_THROW(BadInput,"SABIntegrator invalid energy grid: When emax=0 and table has suggested Emax, the emin value specified must be less than"
                       " this (set emin=0 for automatic emin determination).");
      emax = m_data->suggestedEmax();
    }

    if ( emax == 0.0 ) {

      //The energy associated with the kinematic curve touching (alphamax,betamin)
      //gives us a bound on the upper energy value of the table:
      const double bmin = m_data->betaGrid().front();
      const double amax = m_data->alphaGrid().back();
      const double emax_upper_limit = kT*(bmin-amax)*(bmin-amax)/(4*amax);
      //Try to determine suitable EMax, using very simple algorithm:
      emax = determineEMax(emax_upper_limit);
      if ( !(emax>0.0) ) {
        emax = 0.5 * emax_upper_limit;
        NCRYSTAL_WARN("Algorithm searching for suitable Emax value at"
                      " which to end SAB energy grid failed to provide"
                      " reasonable result. Using crude guess of "<< emax
                      <<"eV. It might be necessary to specify a more"
                      " suitable value directly (e.g. using the \"egrid\""
                      " keyword in .ncmat files). Consider sharing your"
                      " input data with NCrystal developers"
                      " for further debugging.");
      }
    }

    if ( emin == 0.0 ) {
      emin = determineEMin( ncmin(emax*0.01,0.01*kT) );
    } else {
      if ( !(emax>emin) )
        NCRYSTAL_THROW(BadInput,"energy grid does not have emax>emin. Please correct input (possibly by removing hardcoded value of emin).");
    }
    nc_assert_always(emin>0.0);
    nc_assert_always(emax>emin);
    nc_assert_always(npts>=2);
    m_egrid = NC::geomspace(emin,emax,npts);
  }

  if ( m_egrid.size() < 10 )
    NCRYSTAL_THROW(BadInput,"SABIntegrator invalid energy grid - must have at least 10 points.");

  if ( !(m_egrid.front()>0.0) || !nc_is_grid(m_egrid) )
    NCRYSTAL_THROW(BadInput,"SABIntegrator invalid energy grid - must be sorted with non-repeated and positive values.");

}

void NS::SABIntegrator::Impl::doit(SABXSProvider * out_xs, SABSampler* out_sampler, Optional<std::string>* json)
{
  nc_assert_always( out_xs || out_sampler );
  if ( !m_derivedData )
    m_derivedData = s_SABData2DerivedDataFactory.create(D2DDKey(m_data->getUniqueID(),&m_data));

  const bool doSampler = out_sampler!=nullptr;

  //Prepare and validate energy grid:
  setupEnergyGrid();

  SABSampler::SABSamplerAtEList energyPointSamplers;
  if ( doSampler )
    energyPointSamplers.reserve_hint(m_egrid.size());
  VectD xsvals;
  xsvals.reserve(m_egrid.size());

  for (const auto& energy : m_egrid ) {
    nc_assert(energy>0.0);
    auto sampleruptr_and_xs =  analyseEnergyPoint(energy, doSampler );
    if ( doSampler )
      energyPointSamplers.emplace_back(std::move(sampleruptr_and_xs.first));
    xsvals.emplace_back( sampleruptr_and_xs.second );
  }

  energyPointSamplers.shrink_to_fit();

  if ( doSampler )
    out_sampler->setData( m_data->temperature(),
                          VectD(m_egrid.begin(),m_egrid.end()),
                          std::move(energyPointSamplers),
                          m_extender, xsvals.back(), m_egridMargin );
  if ( out_xs )
    out_xs->setData( VectD(m_egrid.begin(),m_egrid.end()),
                     std::move(xsvals),
                     m_extender );

  if (json) {
    std::ostringstream ss;
    {
      std::ostringstream tmp;
      tmp << "nalpha="<<m_data->alphaGrid().size()<<";nbeta="<<m_data->betaGrid().size();
      tmp << ";Emax="<<NeutronEnergy{m_egrid.back()};
      tmp << ";T="<<m_data->temperature();
      tmp << ";M="<<m_data->elementMassAMU();
      tmp << ";sigma_free="<<m_data->boundXS().free(m_data->elementMassAMU());
      streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
    }
    streamJSONDictEntry( ss, "Emax", m_egrid.back()  );
    streamJSONDictEntry( ss, "Emin", m_egrid.front()  );
    streamJSONDictEntry( ss, "negrid", m_egrid.size()  );
    streamJSONDictEntry( ss, "T", m_data->temperature().dbl()  );
    streamJSONDictEntry( ss, "M", m_data->elementMassAMU().dbl()  );
    streamJSONDictEntry( ss, "sigma_bound", m_data->boundXS().dbl()  );
    streamJSONDictEntry( ss, "sigma_free", m_data->boundXS().free(m_data->elementMassAMU()).dbl()  );
    streamJSONDictEntry( ss, "nbeta", m_data->betaGrid().size()  );
    streamJSONDictEntry( ss, "nalpha", m_data->alphaGrid().size(), JSONDictPos::LAST  );
    *json = ss.str();
  }

}

std::pair<NS::SABIntegrator::Impl::SamplerAtE_uptr,double> NS::SABIntegrator::Impl::analyseEnergyPoint(double ekin, bool doSampler ) const
{
  nc_assert_always(ekin>0.0);

  //At each accessible beta-value, figure out the cross-section by integrating
  //over all allowed alpha values at that given (Ekin,beta). The allowed
  //alpha-range will usually encompass a number of whole bins (the "middle"),
  //and have "front" and "back" tail parts into bins just before and after the
  //middle bins.
  //
  //If also setting up sampling, an SamplerAtE instance will also be prepared
  //and returned. It is likely to change in the future, but for now we only
  //support the SABSamplerAtE_Alg1 type, which also needs the same information
  //and breakdown into front/middle/back parts as was used to calculate the
  //cross section.

  const auto& betaGrid = m_data->betaGrid();
  auto alphaGrid_span = Span<const double>(m_data->alphaGrid());
  nc_assert(!!m_derivedData);
  const VectD& logsab = m_derivedData->logsab;
  const VectD& alphaintegrals_cumul = m_derivedData->alphaintegrals_cumul;

  nc_assert(ekin>=0.);
  const double kT = m_data->temperature().kT();
  const double ekin_div_kT = ekin / kT;
  double beta_lower_limit = -ekin_div_kT;
  bool starts_at_kinematic_endpoint = true;
  if (beta_lower_limit<betaGrid.front()) {
    //This can happen at high energies. Push up beta_lower_limit to juuuust
    //before the first grid point, thus keeping code simple without introducing
    //too big artifacts (calculate in three different ways and pick smallest for
    //numerical robustness).
    beta_lower_limit = ncmin( betaGrid.front() - (betaGrid.at(1)-betaGrid.front()) * 1e-6,
                              betaGrid.front() - ncabs(betaGrid.front())*1e-13,
                              std::nexttoward(betaGrid.front(),beta_lower_limit) );
    nc_assert(beta_lower_limit<betaGrid.front());
    starts_at_kinematic_endpoint = false;
  }

  //Figure out which sab alpha-ranges are involved, i.e. intersects the
  //kinematically accessible region (see description of activeGridRanges in
  //NCSABUtils.hh for exactly what is extracted.:
  std::vector<std::pair<uint16_t,uint16_t>> alpharanges;
  std::size_t ibeta_low;
  SABUtils::activeGridRanges( *m_data, ekin_div_kT, alpharanges, ibeta_low );

  if ( ibeta_low >= betaGrid.size() ) {
    //No alpha ranges at all -> cross section is 0 here.
    if (!doSampler)
      return { nullptr, 0.0 };
    return { std::make_unique<SABSamplerAtE_NoScatter>(), 0.0 };
  }

  nc_assert( ibeta_low + alpharanges.size() == betaGrid.size() ) ;

  if ( ibeta_low > 0 && beta_lower_limit < vectAt(betaGrid,ibeta_low-1) ) {
    //move up lower limit, there was apparently no kinematically accessible content further down.
    beta_lower_limit = vectAt(betaGrid,ibeta_low-1);
    starts_at_kinematic_endpoint = false;
  }

  nc_assert( ibeta_low==0 || beta_lower_limit >= betaGrid.front() );
  nc_assert( ibeta_low==0 || beta_lower_limit >= betaGrid.at(ibeta_low-1) );

  //Zoom in on the kinematically allowed part of the beta grid:
  auto relevant_betaGrid = Span<const double>(&betaGrid[0]+ibeta_low,&betaGrid[0] + betaGrid.size());
  nc_assert( alpharanges.size() == (std::size_t)relevant_betaGrid.size() );

  //We know (when beta_lower_limit=-E/kT, otherwise it is just something we
  //assume) that the cross sections take off from 0.0 at beta_lower_limit:
  PairDD prev_betaxs(beta_lower_limit,0.);

  std::vector<PairDD> betasampler_data;
  VectD betasampler_vals,betasampler_weights;
  std::vector<SABSamplerAtE_Alg1::AlphaSampleInfo> sampler_infos;
  const std::size_t nsamplervals = ( doSampler ? relevant_betaGrid.size()+1 : 0 );
  if (doSampler) {
    betasampler_vals.reserve( nsamplervals );
    betasampler_weights.reserve( nsamplervals );
    sampler_infos.reserve(relevant_betaGrid.size());//NB: one less than betasampler_xxx vectors.
    //first point at (beta_lower_limit,xs=0.0)
    betasampler_vals.emplace_back( prev_betaxs.first );
    betasampler_weights.emplace_back( prev_betaxs.second );
  }

  StableSum xs_total_stable;

  bool next_bin_has_kinematic_endpoint = starts_at_kinematic_endpoint;

  for ( auto&& beta : enumerate(relevant_betaGrid) ) {
    if ( beta.val == beta_lower_limit )
      continue;//guard against lower limit falling exactly on a grid point

    nc_assert( beta.val >= beta_lower_limit );
    nc_assert( beta.val >= -ekin_div_kT );
    nc_assert( beta.idx < alpharanges.size() );
    double alow,aupp;
    {
      auto alims = getAlphaLimits(ekin_div_kT, beta.val);
      alow = alims.first;
      aupp = alims.second;
    }

    nc_assert( aupp >= alow );//since beta>=-E/kT this must be the case

    double xs_at_this_beta = 0.0;
    const auto rangeidxs = vectAt( alpharanges, beta.idx );
    unsigned aidx_low(rangeidxs.first), aidx_upp(rangeidxs.second);//unpack shorts to more efficient types


    SABUtils::TailedBreakdown tb;
    if ( aidx_upp > aidx_low && aupp > alow ) {
      //There is non-empty range here and aupp>alow (which can fail for various
      //reasons, including numerical ones).
      const auto nalpha = m_data->alphaGrid().size();
      auto slice_idx = nalpha*(beta.idx + ibeta_low);
      auto sab_slice = Span<const double>(&m_data->sab()[0]+slice_idx,&m_data->sab()[0]+slice_idx+nalpha);
      auto logsab_slice = Span<const double>(&logsab[0]+slice_idx,&logsab[0]+slice_idx+nalpha);
      auto alphaIntegrals_cumul_slice = Span<const double>(&alphaintegrals_cumul[0]+slice_idx,&alphaintegrals_cumul[0]+slice_idx+nalpha);
      tb = SABUtils::createTailedBreakdown( alphaGrid_span, sab_slice, logsab_slice, alphaIntegrals_cumul_slice,
                                            alow, aupp, aidx_low, aidx_upp );
      xs_at_this_beta = tb.xs_front + tb.xs_back + tb.xs_middle;
    } else {
      //No cross-section here!
    }

    nc_assert(xs_at_this_beta>=0.0);
    if (doSampler) {
      sampler_infos.emplace_back();
      auto& info = sampler_infos.back();
      if (xs_at_this_beta > 0.0) {
        info.pt_front.alpha =  tb.front.alpha;
        info.pt_front.sval =  tb.front.sval;
        info.pt_front.logsval =  tb.front.logsval;
        info.pt_back.alpha =  tb.back.alpha;
        info.pt_back.sval =  tb.back.sval;
        info.pt_back.logsval =  tb.back.logsval;
        if (tb.narrow) {
          info.prob_front = 1.0;
        } else {
          info.prob_front = tb.xs_front/xs_at_this_beta;
          info.prob_notback = 1.0 - tb.xs_back/xs_at_this_beta;
          info.pt_front.alpha_idx = tb.imiddle_low;
          info.pt_back.alpha_idx = tb.imiddle_upp;
        }
      } else {
        nc_assert(xs_at_this_beta == 0.0);
        //0 cross-section. For interpolation purposes we fall-back to sampling
        //linearly in allowed alpha range (indicated by special value info.prob_front=2.0)
        info.prob_front = 2.0;
        info.pt_front.alpha =  alow;
        info.pt_back.alpha =  aupp;
      }
    }

    if ( next_bin_has_kinematic_endpoint ) {
      //The piece-wise-linear assumption is too crude in the first bin, we want to
      //use a sqrt(x) shape instead. It just so happen that this increases the
      //integral of the first bin by a factor of 4/3, which can be accommodated
      //(faked!) in the piece-wise-linear machinary by moving the first beta value
      //downwards by 1/3 of the first bin width. The only issue is that when
      //sampling, we should resample values in the first bin according to the
      //sqrt(x) shape, which is a piece of information we need to pass on to the
      //SABSamplerAtE_Alg1 instance.
      next_bin_has_kinematic_endpoint = false;
      double db_real = (beta.val - prev_betaxs.first);
      prev_betaxs.first -= db_real * (1.0 / 3.0);
      if ( doSampler ) {
        nc_assert( betasampler_vals.size() == 1 );
        betasampler_vals.front() = prev_betaxs.first;
      }
    }

    //integrate trapezoidally from (beta,xs) = prev_betaxs
    xs_total_stable.add(0.5 * ( beta.val - prev_betaxs.first ) * ( xs_at_this_beta + prev_betaxs.second ));

    prev_betaxs = { beta.val, xs_at_this_beta };
    if (doSampler) {
      betasampler_vals.emplace_back( prev_betaxs.first );
      betasampler_weights.emplace_back( prev_betaxs.second );
    }
  }

  nc_assert( betasampler_vals.size() == nsamplervals );
  nc_assert( betasampler_weights.size() == nsamplervals );
  nc_assert( !doSampler || sampler_infos.size()+1 ==  nsamplervals );

  //Apply factor C/E, with C=boundXS*kT/4 (cf. eq. 4 in sampling paper):
  double xs_total = xs_total_stable.sum() * m_data->boundXS().get() / (4*ekin_div_kT);

  if (!(xs_total>=0.0))
    xs_total = 0.0;

  if (!doSampler)
    return { nullptr, xs_total };

  if ( xs_total == 0.0 )
    return { std::make_unique<SABSamplerAtE_NoScatter>(), xs_total };

  nc_assert(!!m_derivedData);
  SamplerAtE_uptr up = std::make_unique<SABSamplerAtE_Alg1>( m_derivedData,
                                                             std::move(betasampler_vals),
                                                             std::move(betasampler_weights),
                                                             std::move(sampler_infos),
                                                             ibeta_low,
                                                             ( starts_at_kinematic_endpoint ? beta_lower_limit : 1.0 ) );
  return { std::move(up), xs_total };
}
