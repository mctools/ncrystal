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

#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCMMC_Sphere.hh"
#include "NCMMC_ParseCfg.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      template<class TVolume>
      class GeometryImpl final : public Geometry {
        TVolume m_vol;
      public:

        template<typename ...Args>
        GeometryImpl( Args&& ...args )
          : m_vol(std::forward<Args>(args)... )
        {
        }

        void distToVolumeEntry( const NeutronBasket& nb,
                                Span<double> tgt ) const override
        {
          m_vol.distToVolumeEntry(nb,tgt);
        }

        void distToVolumeExit( const NeutronBasket& nb,
                               Span<double> tgt ) const override
        {
          m_vol.distToVolumeExit(nb,tgt);
        }


        bool pointIsInside( const Vector& v ) const override
        {
          return m_vol.pointIsInside(v);
        }

      };

      GeometryPtr createGeometryImpl( const char * raw_geomstr )
      {
        namespace PMC = parseMMCCfg;

        auto tokens = PMC::tokenize( raw_geomstr );
        auto geom_name =  PMC::mainName( tokens );
        if ( !geom_name.has_value() )
          NCRYSTAL_THROW2(BadInput,"Invalid geom cfg: \""<<raw_geomstr<<"\"");

        if ( geom_name == "sphere" ) {
          PMC::applyDefaults( tokens, "r=0.01" );//0.01m = 1cm
          PMC::checkNoUnknown(tokens,"r","geometry");
          static const int dummy = [](){ Sphere::unit_test(); return 1; }();
          (void)dummy;
          return makeSO<GeometryImpl<Sphere>>( Length{ PMC::getValue_dbl(tokens,"r") } );
        } else {
          NCRYSTAL_THROW2(BadInput,"Unknown geometry type requested: \""<<geom_name<<"\"");
        }
      }
    }
  }
}

NCMMC::GeometryPtr NCMMC::createGeometry( const char * raw_geomstr )
{
  return createGeometryImpl( raw_geomstr );
}
