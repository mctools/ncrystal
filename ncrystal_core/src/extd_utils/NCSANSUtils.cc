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

#include "NCrystal/internal/extd_utils/NCSANSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    const Info::CustomSectionData* getATypePhaseData( const Info& info, const std::string& customsectionname )
    {
      //NB: We don't currently allow A-type phases (the ones with the custom data) to be themselves multiphase.
      nc_assert_always(!info.isMultiPhase());
      auto n = info.countCustomSections(customsectionname);
      if ( n>1 )
        NCRYSTAL_THROW2(BadInput,"Info object (from "<<info.getDataSourceName()
                        <<") contains more than one @CUSTOM_"<<customsectionname<<" section");
      //Respect if sub-phases have ";sans=0" (to e.g. allow usage of a nanodiamond subphase where sans is disabled).
      return ( ( n && Cfg::CfgManip::get_sans(info.getCfgData()) )
               ? &info.getCustomSection(customsectionname)
               : nullptr );
    }

    using CSPDResultsVect = SmallVector<CustomSansPluginData,2>;
    enum class SearchMode { Normal, JustCheckIfPresent };

    CSPDResultsVect extractCDSANSHelper( const Info::PhaseList& phases,
                                         const std::string& customsectionname,
                                         double overallScaleFactor = 1.0,
                                         SearchMode searchMode = SearchMode::Normal )
    {
      CSPDResultsVect res;

      //Total number density normalisation:
      StableSum sum_nd;
      for ( auto& ph : phases )
        sum_nd.add( ph.first * ph.second->getNumberDensity().dbl() );
      nc_assert( sum_nd.sum() > 0.0 );
      const double inv_total_nd = 1.0 / sum_nd.sum();


      //Loop through phases, looking for A-type (the ones with the custom data)
      //and B-type (the ones right after) pairs. Additionally, dive down
      //recursively into any phases which are themselves multiphase:

      auto it = phases.begin();
      const auto itE = phases.end();
      struct AType {
        ScatLenDensity sld;
        NumberDensity nd;
        double volfrac;
        const Info::CustomSectionData* customData = nullptr;
      };
      Optional<AType> currentAType;
      decltype(it) itATypeData = itE;
      for ( ; it!=itE; ++it ) {
        double info_volfrac = it->first;
        const Info& info_ph = *it->second;

        auto aTypeCData = info_ph.isMultiPhase() ? nullptr : getATypePhaseData(it->second,customsectionname);
        if ( currentAType.has_value() ) {
          if ( aTypeCData )
            NCRYSTAL_THROW2(BadInput,"Two consecutive entries in the phase list both contain @CUSTOM_"<<customsectionname<<" sections");
          nc_assert(currentAType.value().customData!=nullptr);
          const double fracA = currentAType.value().volfrac;
          const double fracB = info_volfrac;
          nc_assert_always( fracA > 0.0 && fracB > 0.0 );
          const auto sldA = currentAType.value().sld;
          const auto sldB = info_ph.getSLD();
          const auto drho = sldA.contrast( sldB );
          const auto numberDensity = NumberDensity{ fracA * currentAType.value().nd.dbl() + fracB * info_ph.getNumberDensity().dbl() };

          //For phi, we simply assume that the scattering objects are the
          //smaller (in terms of volume) of the two involved phases. After all,
          //the theoretical SANS models are anyway unlikely to work for samples
          //where the scattering objects occupy more volumes than the solution
          //(e.g. if >50% of the volume is occupied by spheres it is unlikely
          //that approximations regarding diluteness etc. would still hold).
          double phi = std::min<double>( fracA, fracB ) / ( fracA+fracB );

          //Absorp any extra scale factors (<1.0 if the two phases in question
          //are not the only phases) into phi:
          phi *= overallScaleFactor;
          if ( phases.size() != 2 )
            phi *= numberDensity.dbl() * inv_total_nd;

          //Record the entry (if not obviously vanishing):
          auto sansScaleFactor = SANSScaleFactor( drho, numberDensity, phi );
          if ( sansScaleFactor.dbl() > 0 ) {
            res.push_back( { sansScaleFactor, {} } );
            if ( searchMode != SearchMode::JustCheckIfPresent )
              res.back().customData = *(currentAType.value().customData);//copy the whole thing, no life-time hassle for plugin developers
          }
          currentAType = NullOpt;
        } else {
          if ( aTypeCData ) {
            itATypeData = it;
            currentAType = AType{};
            currentAType.value().sld = info_ph.getSLD();
            currentAType.value().nd = info_ph.getNumberDensity();
            currentAType.value().volfrac = info_volfrac;
            currentAType.value().customData = aTypeCData;
          }
        }
        if ( info_ph.isMultiPhase() ) {
          //Search recursively for relevant entries as well, so plugin continues
          //to work when material is used in a composition.
          auto subv = extractCDSANSHelper( info_ph.getPhases(),
                                           customsectionname,
                                           overallScaleFactor * ( info_volfrac * info_ph.getNumberDensity().dbl() * inv_total_nd ),
                                           searchMode );
          for ( auto& e : subv )
            res.push_back( std::move(e) );
        }

      }
      if ( currentAType.has_value())
        NCRYSTAL_THROW2(BadInput,"The phase with a @CUSTOM_"<<customsectionname
                        <<" section must always be followed by another phase (which provides the contrast).");

      return res;
    }
  }
}

bool NC::hasCustomDataForSANSPlugin( const Info& info, const std::string& customsectionname )
{
  nc_assert_always(!customsectionname.empty());
  return ( info.isMultiPhase()
           && !extractCDSANSHelper( info.getPhases(), customsectionname, 1.0, SearchMode::JustCheckIfPresent ).empty() );
}

std::vector<NC::CustomSansPluginData> NC::extractCustomDataForSANSPlugin( const Info& info,
                                                                          const std::string& customsectionname )
{
  nc_assert_always(!customsectionname.empty());
  std::vector<NC::CustomSansPluginData> res;
  if (!info.isMultiPhase())
    return res;
  auto v = extractCDSANSHelper( info.getPhases(),
                                customsectionname );
  if ( !v.empty() ) {
    res.reserve(v.size());
    for ( auto& e : v )
      res.push_back( std::move(e) );
  }

  return res;
}
