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
#include "NCMMC_BasketUtils.hh"
#include "NCMMC_Sphere.hh"//NB: Sphere first, since it contains a few utils.
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
          streamJSONDictEntry( os, "short_description", shortDescription() );
          os << "}}";
        }

        std::string shortDescription() const override
        {
          std::ostringstream os;
          m_vol.shortDescr( os );
          return os.str();
        }

        void toString(std::ostream& os) const override
        {
          m_vol.toCfgString(os);
        }

        void distToVolumeEntry( const NeutronBasket& nb,
                                BasketValBufDbl& tgt,
                                std::size_t offset ) const override
        {
          nc_assert( nb.nused <= basket_N );
          nc_assert( offset < nb.nused);
          BasketUtils::basket_validateIfDbg(nb);
          m_vol.distToVolumeEntry(nb,tgt,offset);
          BasketUtils::basket_validateIfDbg(nb);
        }

        bool hasUnboundedDistToVolExit() const override
        {
          return m_vol.hasUnboundedDistToVolExit();
        }

        void distToVolumeExit( const NeutronBasket& nb,
                               BasketValBufDbl& tgt ) const override
        {
          BasketUtils::basket_validateIfDbg(nb);
          m_vol.distToVolumeExit(nb,tgt);
          BasketUtils::basket_validateIfDbg(nb);
        }


        bool pointIsInside( const Vector& v ) const override
        {
          return m_vol.pointIsInside(v);
        }

      };

      GeometryPtr createGeometryImpl( const StrView& raw_geomstr )
      {
        namespace PMC = parseMMCCfg;

        auto tokeninfo = PMC::tokenize( raw_geomstr );
        auto& tokens = tokeninfo.tokens;
        auto geom_name =  tokeninfo.mainName;
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

NCMMC::GeometryPtr NCMMC::createGeometry( const StrView& raw_geomstr )
{
  return createGeometryImpl( raw_geomstr );
}
//TODO: idea for geom checks: random rays inbound, ensure that if propagating
//the particle forward the distToVolumeEntry value, the next distToVolumeEntry
//should be 0 and distToVolumeExit should be >0 (in non-degenerate cases).

void NCMMC::geometryOptsDocsToJSON( std::ostream& os )
{
  os << "{\"intro_text\":";
  streamJSON(os, "Four MiniMC simulation geometries are currently available,"
             " all representing a single convex volume centered around"
             " the origin and axis aligned where appropriate.");
  os << ",\"geom_list\":[";
  using VV3 = SmallVector<std::array<StrView,3>,9>;
  unsigned sort_key = 0;
  auto addGeom = [&os,&sort_key]( StrView name,
                                  const VV3& params,
                                  StrView descr )
  {
    os << "{\"name\":";
    streamJSON(os,name);
    os<< ",\"descr\":";
    streamJSON(os,descr);
    os<< ",\"params\":";
    streamJSON(os,params);
    os<< ",\"sort_key\":";
    streamJSON(os,sort_key++);
    os<<'}';
  };
  {
    VV3 v;
    v.push_back({"r","0.01","sphere radius [m]."});
    addGeom("sphere",v,
            "A spherical volume centered at (0,0,0).");
  }
  os << ',';
  {
    VV3 v;
    v.push_back({"dx","0.01","half-width [m]."});
    v.push_back({"dy","0.01","half-height [m]."});
    v.push_back({"dz","0.01","half-depth [m]."});
    addGeom("box",v,
            "An axis-aligned box with corners at (+-dx,+-dy,+-dz).");
  }
  os << ',';
  {
    VV3 v;
    v.push_back({"dz","0.01","half-depth [m]."});
    addGeom("slab",v,
            "An unbounded slab bound by planes through (0,0,+-dz) and normals (0,0,1). This is essentially a box with dx=dy=infinity.");
  }
  os << ',';
  {
    VV3 v;
    v.push_back({"r","0.01","cylinder radius radius [m]."});
    v.push_back({"dy","0","Half-length of cylinder (0 means infinitely long) [m]."});
    addGeom("cyl",v,
            "A cylinder around the y-axis. If dy!=0 the cylinder has finite length, with ends at y=+-dy.");
  }
  os << ']';

  int iexample = 0;
  auto addExample = [&iexample,&os]( const char * ex, const char * txt )
  {
    if (iexample++)
      os <<',';
    std::array<StrView,2> v = {ex,txt};
    streamJSON(os,v);
  };
  os << ",\"examples\":[";
  addExample("sphere;r=0.05",
             "A sphere centered at (x,y,z)=(0,0,0) with a radius of 5cm.");
  addExample("cyl;r=0.001",
             "An infinite cylinder around the y-axis with a radius of 1mm.");
  addExample("cyl;r=0.007;dy=0.03",
             "A 6cm long cylinder around the y-axis with radius of 7mm."
             " The cylinder endcaps are located at y=+-10cm.");
  addExample("slab;dz=0.1",
             "An infinite slab bounded by planes at z=+-10cm.");
  addExample("box;dx=0.1;dy=0.2;dz=0.3",
             "An axis-aligned box centered at (x,y,z)=(0,0,0) with the 8"
             " corners located at (x,y,z)=(+-10cm,+-20cm,+-30cm).");
  os <<"]}";
}
