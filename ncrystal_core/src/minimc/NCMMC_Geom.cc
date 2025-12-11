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

#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCMMC_Sphere.hh"
#include "NCMMC_Box.hh"
#include "NCMMC_Slab.hh"
#include "NCMMC_Cyl.hh"
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

        void toJSON( std::ostream& os ) const override
        {
          std::ostringstream cfgstr;
          m_vol.toCfgString(cfgstr);
          streamJSONDictEntry( os, "cfgstr", cfgstr.str(), JSONDictPos::FIRST );
          os << ",\"decoded\":{";
          m_vol.toJSONDecodedItems(os);
          os << "}}";
        }

        void toString(std::ostream& os) const
        {
          m_vol.toCfgString(os);
        }

        void distToVolumeEntry( const NeutronBasket& nb,
                                Span<double> tgt,
                                std::size_t offset ) const override
        {
          //fixme: maybe just handle the offset logic here and then pass on
          //final x,y,z,ux,uy,uz,n
          nc_assert( tgt.size() >= nb.nused);
          nc_assert( offset < nb.nused);
          nb.validateIfDbg();
          m_vol.distToVolumeEntry(nb,tgt,offset);
          nb.validateIfDbg();
        }

        bool hasUnboundedDistToVolExit() const override
        {
          return m_vol.hasUnboundedDistToVolExit();
        }

        void distToVolumeExit( const NeutronBasket& nb,
                               Span<double> tgt ) const override
        {
          nb.validateIfDbg();
          m_vol.distToVolumeExit(nb,tgt);
          nb.validateIfDbg();
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
          return makeSO<GeometryImpl<Sphere>>
            ( Length{ PMC::getValue_dbl(tokens,"r") } );

        } else if ( geom_name == "cyl" ) {
          PMC::applyDefaults( tokens, "r=0.01;dy=0" );//0.01m = 1cm
          PMC::checkNoUnknown(tokens,"r;dy","geometry");
          return makeSO<GeometryImpl<Cyl>>
            ( Length{ PMC::getValue_dbl(tokens,"r") },
              Length{ PMC::getValue_dbl(tokens,"dy") } );
        } else if ( geom_name == "slab" ) {
          PMC::applyDefaults( tokens, "dz=0.01" );//0.01m = 1cm
          PMC::checkNoUnknown(tokens,"dz","geometry");
          return makeSO<GeometryImpl<Slab>>
            ( Length{ PMC::getValue_dbl(tokens,"dz") } );
        } else if ( geom_name == "box" ) {
          PMC::applyDefaults( tokens, "dx=0.01;dy=0.01;dz=0.01" );//0.01m = 1cm
          PMC::checkNoUnknown(tokens,"dx;dy;dz","geometry");
          return makeSO<GeometryImpl<Box>>
            ( Length{ PMC::getValue_dbl(tokens,"dx") },
              Length{ PMC::getValue_dbl(tokens,"dy") },
              Length{ PMC::getValue_dbl(tokens,"dz") } );
        } else {
          NCRYSTAL_THROW2(BadInput,
                          "Unknown geometry type requested: \""
                          <<geom_name<<"\"");
        }
      }
    }
  }
}

NCMMC::GeometryPtr NCMMC::createGeometry( const char * raw_geomstr )
{
  return createGeometryImpl( raw_geomstr );
}
//FIXME: idea for geom checks: random rays inbound, ensure that if propagating
//the particle forward the distToVolumeEntry value, the next distToVolumeEntry
//should be 0 and distToVolumeExit should be >0 (in non-degenerate cases).
