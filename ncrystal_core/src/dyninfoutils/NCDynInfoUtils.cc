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

#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/vdos/NCVDOSCache.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace DICache {

    namespace {
      //VDOSDebye cache key:
      //         (reduced vdoslux 0..2 + rounded: elementMass,boundXS,T,TDebye)
      using VDOSDebyeKey
      = std::tuple<unsigned,uint64_t,uint64_t,uint64_t,uint64_t>;

      //For VDOS Debye we can easily share work between different Info objects
      //or phases, since the number of dependent parameters is very low. Thus,
      //we base the key on the rounded values of those parameters, and make sure
      //we only base calculations on values derived from those rounded values:
      struct VDOSDebyePars {
        unsigned reduced_vdoslux;
        AtomMass elementMass;
        Temperature temperature;
        DebyeTemperature debyeTemperature;
        SigmaBound boundXS;
      };

      VDOSDebyeKey getKey(unsigned nonreduced_vdoslux, Temperature t,
                          DebyeTemperature dt, SigmaBound sb, AtomMass mass ) {
        dt.validate();
        t.validate();
        sb.validate();
        mass.validate();
        nc_assert( sb.get() > 0.0 && sb.get() < 1.0e6 );
        nc_assert( dt.get() > 0.0 && dt.get() < 1.0e5 );
        nc_assert( mass.get() > 0.0 && mass.get() < 1.0e5 );
        nc_assert(nonreduced_vdoslux<=5);
        unsigned reduced_vdoslux = static_cast<unsigned>
          (std::max<int>(0,static_cast<int>(nonreduced_vdoslux)-3));
        nc_assert(reduced_vdoslux<=2);
        auto roundFct = [](double x)
        {
          nc_assert_always(x>0.0&&x<1.0e11);
          return static_cast<uint64_t>(1e7*x+0.5);
        };
        return VDOSDebyeKey( reduced_vdoslux,
                             roundFct(mass.get()),
                             //sb might be zero in rare corner cases:
                             (sb.get()?roundFct(sb.get()):uint64_t(0)),
                             roundFct(t.get()),
                             roundFct(dt.get()) );
      }

      VDOSDebyeKey getKey(unsigned nonreduced_vdoslux, const DI_VDOSDebye& di) {
        return getKey(nonreduced_vdoslux,
                      di.temperature(),
                      di.debyeTemperature(),
                      di.atomData().scatteringXS(),
                      di.atomData().averageMassAMU());
      }

      VDOSDebyePars debyekey2params( const VDOSDebyeKey& key ) {
        //always base calculations only on what can be extracted using the key
        //(this is important due to rounding):
        return { std::get<0>(key),
                 AtomMass{std::get<1>(key)*1e-7},
                 Temperature{std::get<3>(key)*1e-7},
                 DebyeTemperature{std::get<4>(key)*1e-7},
                 SigmaBound{std::get<2>(key)*1e-7} };
      }

      //For regular VDOS expansion, we need a bit more care.
      struct VDOSExpandInput {
        const VDOSData* vdosdata = nullptr;
        double requestedEMax = 0.0;
        double vdos2sabExcludeScaleFactor = 0.0;
        std::uint32_t vdoslux;//0..5
        std::uint32_t vdos2sabExcludeFlag = 0;
      };

      struct VDOSExpandCacheKey {
        UniqueIDValue uid_vdosdata = {0};
        double requestedEMax = 0.0;
        double vdos2sabExcludeScaleFactor = 0.0;
        std::uint32_t vdos2sabExcludeFlag = 0;
        std::uint32_t vdoslux = 3;//0..5
        bool operator<(const VDOSExpandCacheKey& o ) const
        {
          if ( uid_vdosdata != o.uid_vdosdata )
            return uid_vdosdata < o.uid_vdosdata;
          if ( requestedEMax != o.requestedEMax )
            return requestedEMax < o.requestedEMax;
          if ( vdos2sabExcludeFlag != o.vdos2sabExcludeFlag )
            return vdos2sabExcludeFlag < o.vdos2sabExcludeFlag;
          if ( vdos2sabExcludeScaleFactor != o.vdos2sabExcludeScaleFactor )
            return vdos2sabExcludeScaleFactor < o.vdos2sabExcludeScaleFactor;
          return vdoslux < o.vdoslux;
        }

        void adopt( const VDOSExpandInput& i )
        {
          vdos2sabExcludeScaleFactor = i.vdos2sabExcludeScaleFactor;
          vdos2sabExcludeFlag = i.vdos2sabExcludeFlag;
          vdoslux = i.vdoslux;
          requestedEMax = i.requestedEMax;
        }
      };

      double requestedEMaxFromEgrid( const std::shared_ptr<const VectD>& egrid )
      {
        //egrid is either the grid pts directly (when size>3), of the form
        //[emin, emax, npts], where 0 entries (or null ptr) indicates no value.
        if ( !egrid || egrid->empty() )
          return 0.0;
        nc_assert_always(egrid->size()>=3);
        return egrid->size()==3 ? egrid->at(1) : egrid->back();
      }

      shared_obj<const SABData>
      extractFrom_VDOSNoCache( const VDOSExpandInput& input )
      {
        //If user specified an energy-grid with a specific upper energy, Emax,
        //this is essentially a request to expand the vdos out to that energy:
        nc_assert(input.vdosdata);
        const VDOSData& vd = *input.vdosdata;

        ScaleGnContributionFct scaleGnFct = nullptr;
        if ( input.vdos2sabExcludeFlag > 0 ) {
          //vdos2sabExcludeFlag = MODE + 4*LOW + 40000*HIGH
          unsigned high = input.vdos2sabExcludeFlag / 40000;
          unsigned low = (input.vdos2sabExcludeFlag / 4)%10000;
          unsigned mode = input.vdos2sabExcludeFlag % 4;
          if (low>=9999)
            low = std::numeric_limits<unsigned>::max();
          if (high>=9999)
            high = std::numeric_limits<unsigned>::max();
          nc_assert_always(high>=low);
          nc_assert_always(low>=1);
          nc_assert_always(mode>0);
          //Must reduce contribution of selected Gn functions (Gnlow,...,Gnhigh)
          //with scale factor dependening on mode.
          const double scalefact = input.vdos2sabExcludeScaleFactor;
          nc_assert_always( scalefact>=0.0 && scalefact<=1.0 );
          scaleGnFct = [scalefact,low,high](unsigned n)
          {
            return ( n >= low && n<= high ) ? scalefact : 1.0;
          };
        }
        SABData sabdata
          = SABUtils::transformKernelToStdFormat
          ( createScatteringKernel( vd,
                                    static_cast<unsigned>(input.vdoslux),
                                    input.requestedEMax,
                                    VDOSGn::TruncAndThinningChoices::Default,
                                    scaleGnFct ) );
        return std::make_shared<const SABData>(std::move(sabdata));
      }

      using VDOSKeyThin = VDOSExpandCacheKey;
      struct VDOSKeyThick {
        VDOSExpandCacheKey key_thin;
        VDOSExpandInput key_thick;
      };

      struct VDOSKeyThinner {
        using key_type = VDOSKeyThick;
        using thinned_key_type = VDOSKeyThin;
        template <class TMap>
        static typename TMap::mapped_type&
        cacheMapLookup( TMap& map, const key_type& key,
                        Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.key_thin;
          return map[tkey.value()];
        }
      };

      class VDOS2SABFactory final
        : public CachedFactoryBase<VDOSKeyThick, SABData,
                                   30/*NStrongRefsKept*/,VDOSKeyThinner> {
      public:
        const char* factoryName() const override { return "VDOS2SABFactory"; }
        std::string keyToString( const VDOSKeyThick& key_thick ) const override
        {
          const auto& k = key_thick.key_thin;
          std::ostringstream ss;
          ss<<"(VDOSData id="<<k.uid_vdosdata.value
            <<";vdoslux="<<k.vdoslux;
          if ( k.requestedEMax != 0.0 )
            ss<<";requestedEMax="<<k.requestedEMax;
          if ( k.vdos2sabExcludeFlag ) {
            ss<<";vdos2sabExcludeFlag="<<k.vdos2sabExcludeFlag;
            ss<<";vdos2sabScaleFactor="<<k.vdos2sabExcludeScaleFactor;
          } else {
            nc_assert_always(k.vdos2sabExcludeScaleFactor==0.0);
          }
          ss<<")";
          return ss.str();
        }
      protected:
        ShPtr actualCreate( const VDOSKeyThick& key ) const override
        {
          return extractFrom_VDOSNoCache( key.key_thick );
        }
      };

      shared_obj<const SABData>
      extractFrom_DI_VDOS( const DI_VDOS& di_vdos, bool try_use_cache,
                           unsigned vdoslux, uint32_t vdos2sabExcludeFlag )
      {
        nc_assert_always( vdoslux <= 5 );

        VDOSExpandInput input;
        const auto& vd = di_vdos.vdosData();
        input.vdosdata = &vd;
        input.requestedEMax = requestedEMaxFromEgrid( di_vdos.energyGrid() );
        input.vdoslux = static_cast<std::uint32_t>(vdoslux);
        if ( vdos2sabExcludeFlag%4 ) {
          const auto& ad = di_vdos.atomData();
          if ( ad.scatteringXS() != vd.boundXS() )
            NCRYSTAL_THROW(LogicError,"VDOSData from DI_VDOS has boundXS which"
                           " is not consistent with total scatteringXS of"
                           " associated atom");
          input.vdos2sabExcludeFlag = vdos2sabExcludeFlag;
          double scalefact = 0.0;
          //vdos2sabExcludeFlag = MODE + 4*LOW + 40000*HIGH
          unsigned mode = input.vdos2sabExcludeFlag % 4;
          if ( mode == 1 ) {
            //exclude sigma_coh
            scalefact = ad.incoherentXS().dbl() / ad.scatteringXS().dbl();
          } else if ( mode == 2 ) {
            //exclude sigma_incoh
            scalefact = ad.coherentXS().dbl() / ad.scatteringXS().dbl();
          } else {
            nc_assert_always( mode == 3 );
            scalefact = 0.0;//exclude both sigma_coh and sigma_incoh
          }
          nc_assert_always( scalefact>=0.0 && scalefact<=1.0 );
          input.vdos2sabExcludeScaleFactor = scalefact;
        }
        const DI_VDOSShPtr* di_vdosshptr
          = ( try_use_cache
              ? dynamic_cast<const DI_VDOSShPtr*>(&di_vdos) : nullptr );
        if (!di_vdosshptr)
          return extractFrom_VDOSNoCache(input);

        VDOSExpandCacheKey key;
        key.adopt(input);
        key.uid_vdosdata = di_vdosshptr->vdosDataHashPtr().data().getUniqueID();
        VDOSKeyThick cachekey;
        cachekey.key_thin = key;
        cachekey.key_thick = input;
        static VDOS2SABFactory s_vdos2sabfactory;
        return s_vdos2sabfactory.create(cachekey);
      }

      shared_obj<const SABData>
      extractFromDIVDOSDebyeNoCache( const VDOSDebyeKey& key )
      {
        auto param = debyekey2params( key );
        //Setup VDOS data from Debye Model. We only specify points in the upper
        //50% of [0,debye_energy], to benefit from the quadratic scaling below
        //the first grid point implemented in VDOSEval (i.e. we get a more
        //precise G1 function constructed):
        auto vdosdata = createVDOSDebye( param.debyeTemperature,
                                         param.temperature,
                                         param.boundXS, param.elementMass );
        SABData sabdata = SABUtils::transformKernelToStdFormat
          ( createScatteringKernel( vdosdata, param.reduced_vdoslux ) );
        return std::make_shared<const SABData>(std::move(sabdata));
      }

      class VDOSDebye2SABFactory
        : public CachedFactoryBase<VDOSDebyeKey,SABData,10> {
      public:
        const char* factoryName() const final { return "VDOSDebye2SABFactory"; }
        std::string keyToString( const VDOSDebyeKey& key ) const final
        {
          auto p = debyekey2params( key );
          std::ostringstream ss;
          ss<<"(reduced_vdoslux="<<p.reduced_vdoslux
            <<";M="<<p.elementMass
            <<";T="<<p.temperature
            <<";TDebye="<<p.debyeTemperature
            <<";boundXS="<<p.boundXS<<")";
          return ss.str();
        }
      protected:
        virtual ShPtr actualCreate( const VDOSDebyeKey& key ) const final
        {
          return extractFromDIVDOSDebyeNoCache(key);
        }
      };

      shared_obj<const SABData>
      extractFromDIVDOSDebye( const VDOSDebyeKey& key )
      {
        static VDOSDebye2SABFactory s_vdosdebye2sabfactory;
        return s_vdosdebye2sabfactory.create(key);
      }
    }
  }
}

NC::shared_obj<const NC::SABData>
NC::extractSABDataFromVDOSDebyeModel( DebyeTemperature debT,
                                      Temperature temperature,
                                      SigmaBound boundXS,
                                      AtomMass elementMassAMU,
                                      unsigned vdoslux, bool useCache )
{
  auto key = DICache::getKey(vdoslux,temperature,debT,boundXS,elementMassAMU);
  if (!useCache)
    return DICache::extractFromDIVDOSDebyeNoCache(key);
  return DICache::extractFromDIVDOSDebye(key);
}

NC::shared_obj<const NC::SABData>
NC::extractSABDataFromDynInfo( const NC::DI_ScatKnl* di,
                               unsigned vdoslux,
                               bool useCache,
                               std::uint32_t vdos2sabExcludeFlag )
{
  nc_assert( di );
  nc_assert( vdoslux <= 5 );

  //==> VDOSDebye
  auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di);
  if (di_vdosdebye) {
    auto key = DICache::getKey(vdoslux,*di_vdosdebye);
    if (!useCache)
      return DICache::extractFromDIVDOSDebyeNoCache(key);
    return DICache::extractFromDIVDOSDebye(key);
  }

  //==> Directly specified kernels:
  auto di_direct = dynamic_cast<const DI_ScatKnlDirect*>(di);
  if (di_direct)
    return di_direct->ensureBuildThenReturnSAB();

  //==> VDOS:
  auto di_vdos = dynamic_cast<const DI_VDOS*>(di);
  if (di_vdos)
    return DICache::extractFrom_DI_VDOS( *di_vdos, useCache,
                                         vdoslux, vdos2sabExcludeFlag );

  //==> Unknown:
  NCRYSTAL_THROW(LogicError,"Unknown DI_ScatKnl sub class");
  return std::shared_ptr<const SABData>{nullptr};
}

//Idealised VDOS based only on Debye temperature:
NC::VDOSData NC::createVDOSDebye( DebyeTemperature debyeTemperature,
                                  Temperature temperature,
                                  SigmaBound boundXS,
                                  AtomMass elementMassAMU )
{
  //NB: Must keep function exactly synchronised with createVDOSDebye function in
  //Python interface:
  const double debye_energy = constant_boltzmann*debyeTemperature.get();
  auto vdos_egrid = linspace(0.5*debye_energy,debye_energy,20);
  double scale = 1.0 / (debye_energy*debye_energy);
  auto vdos_density = vectorTrf( vdos_egrid,
                                 [scale](double e){ return e*e*scale; } );
  //Actual returned egrid should contain only first and last value:
  return VDOSData( PairDD(vdos_egrid.front(),vdos_egrid.back()),
                   std::move(vdos_density),
                   temperature,
                   boundXS,
                   elementMassAMU );
}
