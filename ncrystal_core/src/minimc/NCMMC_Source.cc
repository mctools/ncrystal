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

#include "NCrystal/internal/minimc/NCMMC_Source.hh"
#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCMMC_ParseCfg.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    namespace {
      class StatCount {
      public:
        StatCount( std::size_t n )
          : m_nused( 0 ),
            m_norig( n )
        {
        }
        void apply(SourceMetaData& md) const
        {
          if ( m_norig == 0 ) {
            md.isInfinite = true;
          } else {
            md.isInfinite = false;
            md.totalSize = m_norig;
          }
        }
        std::size_t nUsed() const
        {
          std::size_t raw = m_nused.load();
          return ( ( m_norig != 0 && raw > m_norig )
                   ? m_norig
                   : raw );
        }
        struct FillCounts {
          std::size_t i0, N;
          bool empty() const { return i0 == N; }
        };
        FillCounts updateBasketCounts( NeutronBasket& nb )
        {
          FillCounts fc;
          fc.i0 = nb.size();
          const std::uint64_t n_needs = (std::uint64_t)(basket_N - fc.i0);
          nc_assert( n_needs>0 );
          const std::uint64_t nused_after_add = m_nused.fetch_add(std::uint64_t(n_needs), std::memory_order_relaxed) + n_needs;
          const std::uint64_t overshot = ( ( m_norig == 0 || nused_after_add < m_norig ) ? 0 : nused_after_add-m_norig );
          const std::uint64_t n_use = ( overshot >= n_needs ? 0 : n_needs - overshot );
          fc.N = fc.i0 + (std::size_t)(n_use);
          nb.nused = fc.N;
          return fc;
        }
      private:
        std::atomic<std::uint64_t> m_nused;
        std::size_t m_norig;
      };

      class SourceIsotropic final : public Source {
        StatCount m_stat;
        Length m_x;
        Length m_y;
        Length m_z;
        NeutronEnergy m_ekin;
      public:

        SourceIsotropic( std::size_t n,
                         NeutronEnergy ekin,
                         Length x = Length{0},
                         Length y = Length{0},
                         Length z = Length{0} )
          : m_stat(n), m_x(x), m_y(y), m_z(z), m_ekin(ekin)
        {
        }

        SourceMetaData metaData() const override
        {

          SourceMetaData md;
          {
            std::ostringstream ss;
            ss << "SourceIsotropic("<<m_ekin<<", pos=["<< m_x<<", "<< m_y<<", "<< m_z<<"])";
            md.description = ss.str();
          }
          md.concurrent = true;
          m_stat.apply( md );
          return md;
        }

        bool particlesMightBeOutside( const Geometry& geom ) const override
        {
          return !geom.pointIsInside( { m_x.dbl(), m_y.dbl(), m_z.dbl() } );
        }

        void fillBasket( RNG& rng, NeutronBasket& nb ) override
        {
          auto counts = m_stat.updateBasketCounts( nb );
          NCRYSTAL_DEBUGMMCMSG("Source Filling n="<<(counts.N-counts.i0)
                               <<" with xyz: " <<m_x<<", "<<m_y<<", "<<m_z);
          for ( std::size_t i = counts.i0; i < counts.N; ++i ) {
            auto v = randIsotropicDirection( rng );
            nb.ux[i] = v.x();
            nb.uy[i] = v.y();
            nb.uz[i] = v.z();
          }
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.x[i] = m_x.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.y[i] = m_y.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.z[i] = m_z.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.w[i] = 1.0;
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.ekin[i] = m_ekin.dbl();
        }
      };

      class SourceConstant final : public Source {
        StatCount m_stat;
        Length m_x;
        Length m_y;
        Length m_z;
        NeutronDirection m_dir;
        NeutronEnergy m_ekin;
      public:

        SourceConstant( std::size_t n,
                        NeutronEnergy ekin,
                        Length x = Length{0},
                        Length y = Length{0},
                        Length z = Length{0},
                        NeutronDirection direction = NeutronDirection{0.0,0.0,1.0} )
          : m_stat(n), m_x(x), m_y(y), m_z(z),
            m_dir( direction.as<Vector>().unit().as<NeutronDirection>() ),
            m_ekin(ekin)
        {
        }

        SourceMetaData metaData() const override
        {

          SourceMetaData md;
          {
            std::ostringstream ss;
            ss << "SourceConstant("<<m_ekin<<", pos=["<< m_x<<", "<< m_y<<", "<< m_z<<"], dir="<<m_dir<<")";
            md.description = ss.str();
          }
          md.concurrent = true;
          m_stat.apply( md );
          return md;
        }

        bool particlesMightBeOutside( const Geometry& geom ) const override
        {
          return !geom.pointIsInside( { m_x.dbl(), m_y.dbl(), m_z.dbl() } );
        }

        void fillBasket( RNG&, NeutronBasket& nb ) override
        {
          auto counts = m_stat.updateBasketCounts( nb );
          NCRYSTAL_DEBUGMMCMSG("Source Filling n="<<(counts.N-counts.i0)
                               <<" with xyz: " <<m_x<<", "<<m_y<<", "<<m_z
                               <<" and dir: "<<m_dir);
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.ux[i] = m_dir[0];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.uy[i] = m_dir[1];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.uz[i] = m_dir[2];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.x[i] = m_x.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.y[i] = m_y.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.z[i] = m_z.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.w[i] = 1.0;
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            nb.ekin[i] = m_ekin.dbl();
        }
      };

      SourcePtr createSourceImpl( const char * raw_srcstr )
      {
        namespace PMC = parseMMCCfg;

        auto tokens = PMC::tokenize( raw_srcstr );
        auto src_name =  PMC::mainName( tokens );
        if ( !src_name.has_value() )
          NCRYSTAL_THROW2(BadInput,"Invalid src cfg: \""<<raw_srcstr<<"\"");

        const char * common_defaults = "n=1e6;w=1.0";
        auto common_default_energy = NeutronWavelength{ 1.8 }.energy();

        if ( src_name == "constant" ) {
          PMC::applyDefaults( tokens, common_defaults );
          PMC::applyDefaults( tokens, "x=0;y=0;z=0;ux=0;uy=0;uz=1;n=1e6" );
          PMC::checkNoUnknown(tokens,"ekin;wl;n;w;;x;y;z;ux;uy;uz","source");
          return makeSO<SourceConstant>( PMC::getValue_sizet(tokens,"n"),
                                         PMC::getValue_Energy( tokens, common_default_energy ),
                                         Length{ PMC::getValue_dbl(tokens,"x") },
                                         Length{ PMC::getValue_dbl(tokens,"y") },
                                         Length{ PMC::getValue_dbl(tokens,"z") },
                                         NeutronDirection{ PMC::getValue_dbl(tokens,"ux"),
                                                           PMC::getValue_dbl(tokens,"uy"),
                                                           PMC::getValue_dbl(tokens,"uz") } );
        } else if ( src_name == "isotropic" ) {
          PMC::applyDefaults( tokens, common_defaults );
          PMC::applyDefaults( tokens, "x=0;y=0;z=0;n=1e6" );
          PMC::checkNoUnknown(tokens,"ekin;wl;n;w;;x;y;z","source");
          return makeSO<SourceIsotropic>( PMC::getValue_sizet(tokens,"n"),
                                          PMC::getValue_Energy( tokens, common_default_energy ),
                                          Length{ PMC::getValue_dbl(tokens,"x") },
                                          Length{ PMC::getValue_dbl(tokens,"y") },
                                          Length{ PMC::getValue_dbl(tokens,"z") } );
        } else {
          NCRYSTAL_THROW2(BadInput,"Unknown source type requested: \""<<src_name<<"\"");
        }
      }
    }
  }
}


NCMMC::SourcePtr NCMMC::createSource( const char * raw_srcstr )
{
  auto src = createSourceImpl( raw_srcstr );
  {
    //Sanity checks:
    auto md = src->metaData();
    if ( md.totalSize.has_value() && md.isInfinite )
      NCRYSTAL_THROW(LogicError,"Inconsistent source metadata:"
                     "infinite sources can not have a totalSize");
  }
  return src;
}
