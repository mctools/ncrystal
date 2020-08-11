////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCDynInfoUtils.hh"
#include "NCFactoryUtils.hh"
#include "NCNeutronSCL.hh"
#include "NCMath.hh"
#include "NCVDOSToScatKnl.hh"
#include "NCSABUtils.hh"
namespace NC = NCrystal;

namespace NCrystal {
  namespace DICache {
    //Cache keys:
    using VDOSKey = std::tuple<uint64_t,unsigned,const DI_VDOS*>;//(DI unique id, vdoslux 0..5, DI object)
    using VDOSDebyeKey = std::tuple<unsigned,uint64_t,uint64_t,uint64_t,uint64_t>;//(reduced vdoslux 0..2 + rounded: elementMass, boundXS, T, TDebye)

    //For VDOS Debye we can potentially share work between different Info
    //objects, since the number of dependent parameters is very low. Thus, we
    //base the key on the rounded values of those parameters, and make sure we
    //only base calculations on values derived from those rounded values:
    struct VDOSDebyePars {
      unsigned reduced_vdoslux;
      double elementMass, temperature, debyeTemperature;
      SigmaBound boundXS;
    };

    VDOSDebyeKey getKey(unsigned reduced_vdoslux, double temperature, double debyeTemperature, SigmaBound boundXS, double elementMassAMU ) {
      nc_assert(reduced_vdoslux<=2);
      auto roundFct = [](double x) { nc_assert(x>0.0&&x<1.0e15); return static_cast<uint64_t>(1000.0*x+0.5); };
      return VDOSDebyeKey( reduced_vdoslux,
                           roundFct(elementMassAMU),
                           roundFct(boundXS.val),
                           roundFct(temperature),
                           roundFct(debyeTemperature) );
    }
    VDOSDebyeKey getKey(unsigned reduced_vdoslux, const DI_VDOSDebye& di) {
      return getKey(reduced_vdoslux,
                    di.temperature(),
                    di.debyeTemperature(),
                    SigmaBound{NeutronSCL::instance()->getBoundXS(di.elementName())},
                    NeutronSCL::instance()->getAtomicMass(di.elementName()));
    }

    VDOSDebyePars debyekey2params( const VDOSDebyeKey& key ) {
      //always base calculations only on what can be extracted using the key
      //(this is important due to rounding):
      return { std::get<0>(key),
               std::get<1>(key)*0.001,
               std::get<3>(key)*0.001,
               std::get<4>(key)*0.001,
               SigmaBound{std::get<2>(key)*0.001} };
    }

    //Actual worker functions producing results:
    std::shared_ptr<const SABData> extractFromDIVDOSNoCache( unsigned vdoslux, const DI_VDOS& );
    std::shared_ptr<const SABData> extractFromDIVDOSDebyeNoCache( const VDOSDebyeKey& );

    //Factories:
    class VDOS2SABFactory : public NC::CachedFactoryBase<VDOSKey,SABData> {
    public:
      const char* factoryName() const final { return "VDOS2SABFactory"; }
      std::string keyToString( const VDOSKey& key ) const final
      {
        std::ostringstream ss;
        ss<<"(DI_VDOS id="<<std::get<0>(key)<<";vdoslux="<<std::get<1>(key)<<")";
        return ss.str();
      }
    protected:
      virtual ShPtr actualCreate( const VDOSKey& key ) final
      {
        unsigned vdoslux = std::get<1>(key);
        const DI_VDOS* di_vdos = std::get<2>(key);
        nc_assert_always( di_vdos && di_vdos->getUniqueID().value == std::get<0>(key) );
        return extractFromDIVDOSNoCache( vdoslux, *di_vdos );
      }
    };

    class VDOSDebye2SABFactory : public NC::CachedFactoryBase<VDOSDebyeKey,SABData> {
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
          <<";boundXS="<<p.boundXS.val<<")";
        return ss.str();
      }
    protected:
      virtual ShPtr actualCreate( const VDOSDebyeKey& key ) final
      {
        return extractFromDIVDOSDebyeNoCache(key);
      }
    };

    static VDOS2SABFactory s_vdos2sabfactory;
    static VDOSDebye2SABFactory s_vdosdebye2sabfactory;

    std::shared_ptr<const SABData> extractFromDIVDOS( unsigned vdoslux, const DI_VDOS& di )
    {
      VDOSKey key( di.getUniqueID().value, vdoslux, &di );
      return s_vdos2sabfactory.create(key);
    }

    std::shared_ptr<const SABData> extractFromDIVDOSDebye( const VDOSDebyeKey& key )
    {
      return s_vdosdebye2sabfactory.create(key);
    }

  }
}

std::shared_ptr<const NC::SABData> NC::extractSABDataFromVDOSDebyeModel( double debyeTemperature,
                                                                         double temperature,
                                                                         SigmaBound boundXS,
                                                                         double elementMassAMU,
                                                                         unsigned vdoslux, bool useCache )
{
  nc_assert( vdoslux <= 5 );
  nc_assert( temperature > 0.0 && temperature < 1.0e5 );
  nc_assert( boundXS.val > 0.0 && boundXS.val < 1.0e6 );
  nc_assert( debyeTemperature > 0.0 && debyeTemperature < 1.0e5 );
  nc_assert( elementMassAMU > 0.0 && elementMassAMU < 1.0e5 );
  unsigned reduced_vdoslux = static_cast<unsigned>(std::max<int>(0,static_cast<int>(vdoslux)-3));//nb: replicated below
  auto key = DICache::getKey(reduced_vdoslux,temperature,debyeTemperature,boundXS,elementMassAMU);
  if (!useCache)
    return DICache::extractFromDIVDOSDebyeNoCache(key);
  return DICache::extractFromDIVDOSDebye(key);
}

std::shared_ptr<const NC::SABData> NC::extractSABDataFromDynInfo( const NC::DI_ScatKnl* di, unsigned vdoslux, bool useCache )
{
  nc_assert( di );
  nc_assert( vdoslux <= 5 );

  //==> VDOSDebye
  auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di);
  if (di_vdosdebye) {
    unsigned reduced_vdoslux = static_cast<unsigned>(std::max<int>(0,static_cast<int>(vdoslux)-3));//nb: replicated above
    auto key = DICache::getKey(reduced_vdoslux,*di_vdosdebye);
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
  if (di_vdos) {
    if (!useCache)
      return DICache::extractFromDIVDOSNoCache(vdoslux,*di_vdos);
    return DICache::extractFromDIVDOS(vdoslux,*di_vdos);
  }

  //==> Unknown:
  NCRYSTAL_THROW(LogicError,"Unknown DI_ScatKnl sub class");
  return nullptr;
}

void NC::clearSABDataFromDynInfoCaches()
{
  DICache::s_vdos2sabfactory.cleanup();
  DICache::s_vdosdebye2sabfactory.cleanup();
}

std::shared_ptr<const NC::SABData> NC::DICache::extractFromDIVDOSNoCache( unsigned vdoslux, const DI_VDOS& di )
{
  //If user specified an energy-grid with a specific upper energy, Emax,
  //this is essentially a request to expand the vdos out to that energy:
  double requested_Emax = 0.0;
  auto egrid = di.energyGrid();
  if ( !!egrid && ! egrid->empty() ) {
    //egrid is either the grid pts directly (when size>3), of the form
    //[emin, emax, npts], where 0 entries indicates no value.
    nc_assert_always(egrid->size()>=3);
    requested_Emax = egrid->size()==3 ? egrid->at(1) : egrid->back();
  }
  const auto& vd = di.vdosData();
  SABData sabdata = SABUtils::transformKernelToStdFormat( createScatteringKernel( vd, vdoslux,requested_Emax ) );
  return std::make_shared<const SABData>(std::move(sabdata));

}

std::shared_ptr<const NC::SABData> NC::DICache::extractFromDIVDOSDebyeNoCache( const VDOSDebyeKey& key )
{
  auto param = debyekey2params( key );

  //Setup VDOS data from Debye Model. We only specify points in the upper 50% of
  //[0,debye_energy], to benefit from the quadratic scaling below the first grid
  //point implemented in VDOSEval (i.e. we get a more precise G1 function
  //constructed):
  auto vdosdata = createVDOSDebye( param.debyeTemperature, param.temperature, param.boundXS, param.elementMass );
  SABData sabdata = SABUtils::transformKernelToStdFormat( createScatteringKernel( vdosdata, param.reduced_vdoslux ) );
  return std::make_shared<const SABData>(std::move(sabdata));
}


//Idealised VDOS based only on Debye temperature:
NC::VDOSData NC::createVDOSDebye(double debyeTemperature, double temperature, SigmaBound boundXS, double elementMassAMU)
{
  //NB: Must keep function exactly synchronised with createVDOSDebye function in Python interface:
  const double debye_energy = constant_boltzmann*debyeTemperature;
  auto vdos_egrid = linspace(0.5*debye_energy,debye_energy,20);
  double scale = 1.0 / (debye_energy*debye_energy);
  auto vdos_density = vectorTrf(vdos_egrid,[scale](double e){ return e*e*scale; });
  //Actual returned egrid should contain only first and last value:
  return VDOSData( PairDD(vdos_egrid.front(),vdos_egrid.back()),
                   std::move(vdos_density),
                   temperature,
                   boundXS,
                   elementMassAMU );
}
