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

#include "NCrystal/internal/minimc/NCMMC_Source.hh"
#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/extd_utils/NCABIUtils.hh"
#include "NCMMC_ParseCfg.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    namespace {
      using EParsed = parseMMCCfg::EParsed;
      using FlexRangeValue = parseMMCCfg::FlexRangeValue;

      class StatCount {
      public:
        StatCount( std::size_t n )
          : m_nused( 0 ),
            m_norig( n )
        {
        }
        std::size_t nOrig() const { return m_norig; }
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

      void setEnergy_Maxwell( double halfkT,
                              const StatCount::FillCounts& counts,
                              RNG& rng,
                              BasketValBufDbl& buf_ekin )
      {
        //Essentially sample 3 std gaussian numbers, G1, G2, G3, then
        //
        //  E = ( G1^2+G2^2+G3^2 )*kT/2
        //
        //(the physics here is that we sample three velocity components
        //independently, hence 3 gaussian numbers)

        nc_assert( counts.N > counts.i0 );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] = randNorm( rng );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] *= buf_ekin[i];
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] += ncsquare( randNorm( rng ) );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] += ncsquare( randNorm( rng ) );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] *= halfkT;
      }

      void convertBufWl2E( const StatCount::FillCounts& counts,
                           BasketValBufDbl& buf )
      {
        //Converts lambda -> E, as E=k/lambda^2 in a way which handles
        //infinities and zero divisions appropriately and without triggering
        //zero division FPEs.
        nc_assert( counts.N > counts.i0 );

        //square:
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf.data[i] *= buf[i];

        //Safe reciprocal which does not trigger zero division FPE and correctly
        //maps 0 -> inf:
        constexpr double tiny = std::numeric_limits<double>::denorm_min();
        nc_assert( ncisinf(1.0/tiny) );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf.data[i] = 1.0/ncmax(tiny,buf[i]);

        //Final factor:
        constexpr double conv_factor = wlsq2ekin(1.0);
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf.data[i] *= conv_factor;
      }

      void setMetaData_EnergyInfo( const EParsed& ecfg, SourceMetaData& md )
      {
        md.energyDescription = ecfg.description;
        if ( ecfg.fixed_ekin.has_value() ) {
          NeutronEnergy E = ecfg.fixed_ekin.value();
          md.fixedEnergy = E;
          md.meanEnergy = E;
          md.approxERange.emplace( E, E );
          return;
        }
        //Ok, so energy is not fixed but we can set the mean energy (or
        //wavelength if the wl=... keyword was used).
        md.fixedEnergy = NullOpt;
        md.meanEnergy = NullOpt;
        md.approxERange = NullOpt;
        if ( ecfg.maxwell.has_value() ) {
          //mean energy is (3/2)kT, most probable is (1/2)kT. Although a value
          //of kT might feel natural, we use the mean energy for consistency
          //(also for consistency with the documentation):
          const double kT = ecfg.maxwell.value().kT();
          md.meanEnergy = kT*1.5;
          //FIXME: Verify that the values correspond to a 0.997300203936740
          //confidence interval (i.e. +-3sigma std dev).
          md.approxERange.emplace( 0.073*kT, 13.688*kT );
          return;
        }
        auto& fr = ecfg.flexrange.value();
        nc_assert( fr.fr.mode == FlexRangeValue::Mode::UniformRange
                   || fr.fr.mode == FlexRangeValue::Mode::LogNormal );
        //Take the middle of UniformRange or mean of LogNormal, in either
        //wavelength or energy depending on the mode:
        nc_assert( fr.fr.mode == FlexRangeValue::Mode::LogNormal
                   || fr.fr.mode == FlexRangeValue::Mode::UniformRange );
        const double meanval = ( fr.fr.mode == FlexRangeValue::Mode::LogNormal
                                 ? fr.fr.value
                                 : 0.5*( fr.fr.value
                                         + fr.fr.secondary_value.value() ) );
        if ( fr.mode == EParsed::Mode::Energy )
          md.meanEnergy = NeutronEnergy( meanval );
        else
          md.meanEnergy = NeutronWavelength( meanval );
        double a,b;
        if ( fr.fr.mode == FlexRangeValue::Mode::UniformRange ) {
          //The range is trivial in case of UniformRange:
          a = fr.fr.value;
          b = fr.fr.secondary_value.value();
        } else {
          //For log-normal we take the +-3 sigma confidence interval of the
          //underlying gaussian:
          const double mu_n = ecfg.cachevals.first;
          const double sigma_n = ecfg.cachevals.second;
          // const double mu    = fr.fr.value;
          // const double sigma = fr.fr.secondary_value.value();
          a = std::exp( mu_n - 3*sigma_n );
          b = std::exp( mu_n + 3*sigma_n );
        }
        if ( fr.mode == EParsed::Mode::Wavelength )
          md.approxERange.emplace( NeutronWavelength(a),
                                   NeutronWavelength(b) );
        else
          md.approxERange.emplace( NeutronEnergy(a),
                                   NeutronEnergy(b) );
      }

      void setEnergy_UniformRange( const EParsed& ecfg,
                                   const StatCount::FillCounts& counts,
                                   RNG& rng,
                                   BasketValBufDbl& buf_ekin )
      {
        auto& fr = ecfg.flexrange.value();
        nc_assert( fr.fr.mode == FlexRangeValue::Mode::UniformRange );

        const double genval1 = fr.fr.value;
        const double genval2 = fr.fr.secondary_value.value();
        const double delta_genval = genval2 - genval1;
        nc_assert( delta_genval > 0 );
        nc_assert( genval2 > genval1 );

        //uniform from genval1 to genval2:
        nc_assert( counts.N > counts.i0 );
        std::size_t nvals = counts.N - counts.i0;
        nc_assert( nvals <= basket_N );
        NewABI::generateMany( rng, nvals, buf_ekin.data + counts.i0 );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] *= delta_genval;
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] += genval1;

        //clamp upper edge as added safety:
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] = ncmin( genval2, buf_ekin[i] );

        if ( fr.mode == EParsed::Mode::Wavelength )
          convertBufWl2E( counts, buf_ekin );

#ifndef NDEBUG
        static int validate = []()
        {
          BasketValBufDbl vbuf;
          vbuf.data[0] = 0.0;
          vbuf.data[1] = std::numeric_limits<double>::infinity();
          vbuf.data[2] = 1.8;
          StatCount::FillCounts vcounts;
          vcounts.i0 = 0;
          vcounts.N = 3;
          convertBufWl2E( vcounts, vbuf );
          nc_assert( !(ncisnan(vbuf[0])||ncisnan(vbuf[1])||ncisnan(vbuf[2])) );
          nc_assert( ncisinf(vbuf[0]) );
          nc_assert( vbuf[1] == 0 );
          nc_assert( floateq(ekin2wl(vbuf[2]),1.8) );
          return 1;
        }();
        (void)validate;
#endif
      }

      void setEnergy_LogNormal( const EParsed& ecfg,
                                const StatCount::FillCounts& counts,
                                RNG& rng,
                                BasketValBufDbl& buf_ekin )
      {
        auto& fr = ecfg.flexrange.value();
        const double& mu_n = ecfg.cachevals.first;
        const double& sigma_n = ecfg.cachevals.second;
        nc_assert( counts.N > counts.i0 );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] = randNorm( rng );
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] *= sigma_n;
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] += mu_n;
        for ( std::size_t i = counts.i0; i < counts.N; ++i )
          buf_ekin.data[i] = std::exp(buf_ekin[i]);
        if ( fr.mode == EParsed::Mode::Wavelength )
          convertBufWl2E( counts, buf_ekin );
      }

      void setEnergy( const EParsed& ecfg,
                      const StatCount::FillCounts& counts,
                      RNG& rng,
                      BasketValBufDbl& buf_ekin ) {
        if ( ecfg.fixed_ekin.has_value() ) {
          //Fast path for fixed energy:
          nc_assert( ecfg.flexrange.value().fr.mode
                     == FlexRangeValue::Mode::Fixed );
          const double eval = ecfg.fixed_ekin.value().get();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            buf_ekin.data[i] = eval;
          return;
        }
        if ( ecfg.maxwell.has_value() ) {
          nc_assert( floateq( ecfg.maxwell.value().kT()*0.5,
                              ecfg.cachevals.first ) );
          setEnergy_Maxwell( ecfg.cachevals.first, counts, rng, buf_ekin);
        } else {
          auto& fr = ecfg.flexrange.value();
          if (fr.fr.mode == FlexRangeValue::Mode::UniformRange ) {
            setEnergy_UniformRange(ecfg,counts,rng,buf_ekin);
          } else {
            nc_assert( fr.fr.mode == FlexRangeValue::Mode::LogNormal );
            setEnergy_LogNormal(ecfg,counts,rng,buf_ekin);
          }
        }
      }
      void sourceJSONHelper_energyItems( std::ostream& os, const EParsed& ep )
      {
        if ( ep.maxwell.has_value() ) {
          os << ",\"wl\":null,\"ekin\":";
          os << "{\"mode\":\"maxwell\",\"temperature\":";
          streamJSON(os,ep.maxwell.value().get());
          os << '}';
        } else {
          auto& fr = ep.flexrange.value();
          if ( fr.mode == EParsed::Mode::Energy ) {
            os << ",\"wl\":null,\"ekin\":";
            fr.fr.toJSON(os);
          } else {
            os << ",\"wl\":";
            fr.fr.toJSON(os);
            os << ",\"ekin\":null";
          }
        }
      }

      void sourceJSONHelper_direction( std::ostream& os, NeutronDirection dir )
      {
        os << ",\"ux\":";
        streamJSON(os,dir[0]);
        os << ",\"uy\":";
        streamJSON(os,dir[1]);
        os << ",\"uz\":";
        streamJSON(os,dir[2]);
      }

      void sourceJSONHelper_namexyzwn( std::ostream& os,
                                       const char * name,
                                       Length x,
                                       Length y,
                                       Length z,
                                       double w,
                                       const StatCount& stat )
      {
        os << "\"name\":";
        streamJSON(os,name);
        os << ",\"n\":";
        streamJSON(os,stat.nOrig());
        os << ",\"x\":";
        streamJSON(os,x.get());
        os << ",\"y\":";
        streamJSON(os,y.get());
        os << ",\"z\":";
        streamJSON(os,z.get());
        os << ",\"w\":";
        streamJSON(os,w);
      }

      void sourceJSONHelper_md( std::ostream& os, const SourceMetaData& md )
      {
        //NB: skipping md.concurrent for now
        streamJSONDictEntry( os, "fixed_dir",md.fixedDirection,
                             JSONDictPos::FIRST );
        streamJSONDictEntry( os, "fixed_energy",md.fixedEnergy );
        streamJSONDictEntry( os, "mean_dir",md.meanDirection );
        streamJSONDictEntry( os, "mean_energy",md.meanEnergy );
        streamJSONDictEntry( os, "energy_description",md.energyDescription,
                             JSONDictPos::LAST );
      }

      template <class TSource>
      void sourceJSONHelper(std::ostream& os, const TSource& src )
      {
        std::ostringstream cfgstr;
        src.toString(cfgstr);
        streamJSONDictEntry( os, "cfgstr", cfgstr.str(), JSONDictPos::FIRST );
        os << ",\"decoded\":{";
        src.toJSONDecodedCfgItems(os);
        os << "},\"metadata\":";
        sourceJSONHelper_md(os,src.metaData());
        os << '}';
      }

      void sourceCfgStrHelper_direction( std::ostream& os,
                                         NeutronDirection dir )
      {
        if ( dir == NeutronDirection{0,0,1} )
          return;
        os << ";ux="<<fmt(dir[0]);
        os << ";uy="<<fmt(dir[1]);
        os << ";uz="<<fmt(dir[2]);
      }

      void sourceCfgStrHelper_energyItems( std::ostream& os,
                                           const EParsed& ep )
      {
        if ( ep.maxwell.has_value() ) {
          os << ";ekin=thermal:"<<fmt(ep.maxwell.value().get())<<'K';
        } else {
          auto& fr = ep.flexrange.value();
          os << ( fr.mode == EParsed::Mode::Energy ? ";ekin=" : ";wl=" );
          fr.fr.toString(os);
        }
      }

      void sourceCfgStrHelper_namexyzwn( std::ostream& os,
                                         const char * name,
                                         Length x,
                                         Length y,
                                         Length z,
                                         double w,
                                         const StatCount& stat)
      {
        os << name;
        os << ";n="<<stat.nOrig();
        if ( x.dbl() )
          os << ";x="<<fmt(x.dbl());
        if ( y.dbl() )
          os << ";y="<<fmt(x.dbl());
        if ( z.dbl() )
          os << ";z="<<fmt(x.dbl());
        nc_assert_always(w==1.0);
        // if ( w != 1.0 )
        //   os << ";w="<<fmt(w);
      }

      class SourceIsotropic final : public Source {
        StatCount m_stat;
        Length m_x;
        Length m_y;
        Length m_z;
        double m_w = 1.0;
        EParsed m_einfo;
        double m_minusr = 0.0;
        SourceMetaData m_md;
     public:

        //Note: If radius is set to a positive value, we move particles
        //backwards that amount to effectively get a sphere of particles
        //radiating inwards.  If set to a negative value, move particles forward
        //that amount (i.e. a sphere radiating outwards).

        SourceIsotropic( std::size_t n,
                         EParsed ekin_parsed,
                         double weight,
                         Length x = Length{0},
                         Length y = Length{0},
                         Length z = Length{0},
                         double radius = 0.0 )
          : m_stat(n), m_x(x), m_y(y), m_z(z),
            m_w(weight),
            m_einfo(std::move(ekin_parsed)),
            m_minusr(radius?-radius:0.0)
        {
          nc_assert_always(m_w>0.0&&std::isfinite(m_w));
          nc_assert_always(std::isfinite(radius)&&ncabs(radius)<1e99);
          m_md.concurrent = true;
          m_md.meanDirection = NullOpt;//not well defined
          m_md.fixedDirection = NullOpt;//not well defined
          setMetaData_EnergyInfo( m_einfo, m_md );
        }

        void toJSONDecodedCfgItems(std::ostream& os) const
        {
          sourceJSONHelper_namexyzwn(os,"isotropic",m_x,m_y,m_z,m_w,m_stat);
          sourceJSONHelper_energyItems(os,m_einfo);
          os << ",\"r\":";
          streamJSON(os,(m_minusr?-m_minusr:0.0));
        }

        void toJSON(std::ostream& os) const override
        {
          return sourceJSONHelper(os,*this);
        }

        void toString(std::ostream& os) const override
        {
          sourceCfgStrHelper_namexyzwn(os, "isotropic",
                                       m_x, m_y, m_z, m_w, m_stat );
          sourceCfgStrHelper_energyItems( os, m_einfo );
          if ( m_minusr != 0.0 )
            os << ";r="<<fmt((m_minusr?-m_minusr:0.0));
        }

        ParticleCountSum particlesProvided() const override {
          ParticleCountSum s;
          s.count = m_stat.nUsed();
          s.weight = m_stat.nUsed() * m_w;
          return s;
        }

        const SourceMetaData& metaData() const override
        {
          return m_md;
        }

        bool particlesMightBeOutside( const Geometry& geom ) const override
        {
          return !geom.pointIsInside( { m_x.dbl(), m_y.dbl(), m_z.dbl() } );
        }

        void fillBasket( RNG& rng, NeutronBasket& nb ) override
        {
          nc_assert( !nb.full() );
          auto counts = m_stat.updateBasketCounts( nb );
          NCRYSTAL_DEBUGMMCMSG("Source Filling n="<<(counts.N-counts.i0)
                               <<" with xyz: " <<m_x<<", "<<m_y<<", "<<m_z);
          if ( counts.N == counts.i0 )
            return;
          nc_assert( counts.N > counts.i0 );
          auto& f = nb.fields;
          for ( std::size_t i = counts.i0; i < counts.N; ++i ) {
            auto v = randIsotropicDirection( rng );
            f.ux.data[i] = v.x();
            f.uy.data[i] = v.y();
            f.uz.data[i] = v.z();
          }
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.x.data[i] = m_x.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.y.data[i] = m_y.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.z.data[i] = m_z.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.w.data[i] = m_w;
          setEnergy( m_einfo, counts, rng, f.ekin );
          if ( m_minusr ) {
            //Move particles u*(-r):
            for ( std::size_t i = counts.i0; i < counts.N; ++i )
              f.x.data[i] += f.ux[i] * m_minusr;
            for ( std::size_t i = counts.i0; i < counts.N; ++i )
              f.y.data[i] += f.uy[i] * m_minusr;
            for ( std::size_t i = counts.i0; i < counts.N; ++i )
              f.z.data[i] += f.uz[i] * m_minusr;
          }
        }
      };

      class SourceConstant final : public Source {
        //A source which always fires the same neutron state.
        StatCount m_stat;
        Length m_x;
        Length m_y;
        Length m_z;
        NeutronDirection m_dir;
        double m_w = 1.0;
        EParsed m_einfo;
        SourceMetaData m_md;
      public:

        SourceConstant( std::size_t n,
                        EParsed ekin_parsed,
                        double weight,
                        Length x,
                        Length y,
                        Length z,
                        NeutronDirection direction )
          : m_stat(n), m_x(x), m_y(y), m_z(z),
            m_dir( direction.as<Vector>().unit().as<NeutronDirection>() ),
            m_w(weight),
            m_einfo(std::move(ekin_parsed))
        {
          nc_assert_always(m_w>0.0&&std::isfinite(m_w));
          m_md.concurrent = true;
          m_md.fixedDirection = m_dir;
          m_md.meanDirection = m_dir;
          setMetaData_EnergyInfo( m_einfo, m_md );
        }

        void toJSONDecodedCfgItems(std::ostream& os) const
        {
          sourceJSONHelper_namexyzwn(os,"constant",m_x,m_y,m_z,m_w,m_stat);
          sourceJSONHelper_energyItems(os,m_einfo);
          sourceJSONHelper_direction(os,m_dir);
        }

        void toJSON(std::ostream& os) const override
        {
          return sourceJSONHelper(os,*this);
        }

        void toString(std::ostream& os) const override
        {
          sourceCfgStrHelper_namexyzwn( os, "constant",
                                        m_x, m_y, m_z, m_w, m_stat );
          sourceCfgStrHelper_energyItems( os, m_einfo );
          sourceCfgStrHelper_direction( os, m_dir );
        }

        ParticleCountSum particlesProvided() const override {
          ParticleCountSum s;
          s.count = m_stat.nUsed();
          s.weight = m_stat.nUsed() * m_w;
          return s;
        }

        const SourceMetaData& metaData() const override
        {
          return m_md;
        }

        bool particlesMightBeOutside( const Geometry& geom ) const override
        {
          return !geom.pointIsInside( { m_x.dbl(), m_y.dbl(), m_z.dbl() } );
        }

        void fillBasket( RNG& rng, NeutronBasket& nb ) override
        {
          nc_assert( !nb.full() );
          auto counts = m_stat.updateBasketCounts( nb );
          NCRYSTAL_DEBUGMMCMSG("Source Filling n="<<(counts.N-counts.i0)
                               <<" with xyz: " <<m_x<<", "<<m_y<<", "<<m_z
                               <<" and dir: "<<m_dir);
          if ( counts.N == counts.i0 )
            return;
          nc_assert( counts.N > counts.i0 );
          auto& f = nb.fields;
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.ux.data[i] = m_dir[0];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.uy.data[i] = m_dir[1];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.uz.data[i] = m_dir[2];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.x.data[i] = m_x.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.y.data[i] = m_y.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.z.data[i] = m_z.dbl();
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.w.data[i] = m_w;
          setEnergy( m_einfo, counts, rng, f.ekin );
        }
      };

      class SourceUniformCircularBeam final : public Source {
        //A source which fires a monochromatic beam of neutrons, with a circular
        //and uniform beam profile. Energy and radius are both required
        //parameters.
        StatCount m_stat;
        //Length m_x, m_y, m_z;//beam center
        Vector m_center;//beam center
        NeutronDirection m_dir;//beam direction
        Vector m_a, m_b;//basis vectors orthogonal to m_dir;
        double m_w = 1.0;
        Length m_radius;//for metadata
        EParsed m_einfo;
        SourceMetaData m_md;
      public:

        SourceUniformCircularBeam( std::size_t n,
                                   EParsed ekin_parsed,
                                   double weight,
                                   Length radius,
                                   Length x,
                                   Length y,
                                   Length z,
                                   NeutronDirection direction )
          : m_stat(n), m_center{ x.get(), y.get(), z.get() },
            m_dir( direction.as<Vector>().unit().as<NeutronDirection>() ),
            m_w(weight),
            m_radius(radius),
            m_einfo(std::move(ekin_parsed))
        {
          nc_assert_always(m_w>0.0&&std::isfinite(m_w));
          const Vector& vdir = m_dir.as<Vector>();
          if ( m_radius.get() > 0.0 ) {
            //We must find m_a and m_b as basis vectors of the disk on which to
            //generate neutrons:
            //-> First find a unit vector not too colinear with dir:
            m_a.set( 1.0, 0.0, 0.0 );
            if ( m_a.dot(vdir) > 0.8 )
              m_a.set( 0.0, 1.0, 0.0 );
            if ( m_a.dot(vdir) > 0.8 )
              m_a.set( 0.0, 0.0, 1.0 );
            nc_assert_always( m_a.dot(vdir) <= 0.8 );
            //Now use the cross product to ensure full orthogonality:
            m_a = vdir.cross(m_a).unit();
            //m_b is now easy:
            m_b = vdir.cross(m_a).unit();
            nc_assert_always( vdir.isOrthogonal( m_a ) );
            nc_assert_always( vdir.isOrthogonal( m_b ) );
            nc_assert_always( m_a.isOrthogonal( m_b ) );
            //Finally apply radius
            m_a *= m_radius.get();
            m_b *= m_radius.get();
          }
          m_md.concurrent = true;
          m_md.fixedDirection = m_dir;
          m_md.meanDirection = m_dir;
          setMetaData_EnergyInfo( m_einfo, m_md );
          NCRYSTAL_DEBUGMMCMSG("SourceUniformCircularBeam(center="
                               <<m_center<<", dir="<<m_dir
                               <<", energy="<<m_einfo.description
                               <<", r="<<m_radius
                               <<", a="<<m_a<<", b="<<m_b);
        }

        void toJSONDecodedCfgItems(std::ostream& os) const
        {
          sourceJSONHelper_namexyzwn(os,"circular",
                                     Length{m_center[0]},
                                     Length{m_center[1]},
                                     Length{m_center[2]},
                                     m_w,m_stat);
          sourceJSONHelper_energyItems(os,m_einfo);
          sourceJSONHelper_direction(os,m_dir);
          os << ",\"r\":";
          streamJSON(os,m_radius.get());
        }

        void toJSON(std::ostream& os) const override
        {
          return sourceJSONHelper(os,*this);
        }

        void toString(std::ostream& os) const override
        {
          sourceCfgStrHelper_namexyzwn(os, "circular",
                                       Length{m_center[0]},
                                       Length{m_center[1]},
                                       Length{m_center[2]},
                                       m_w, m_stat );
          sourceCfgStrHelper_energyItems( os, m_einfo );
          sourceCfgStrHelper_direction( os, m_dir );
          //radius is always required, hence must always be here:
          os << ";r="<<fmt(m_radius.get());
        }

        ParticleCountSum particlesProvided() const override {
          ParticleCountSum s;
          s.count = m_stat.nUsed();
          s.weight = m_stat.nUsed() * m_w;
          return s;
        }

        const SourceMetaData& metaData() const override
        {
          return m_md;
        }

        bool particlesMightBeOutside( const Geometry& geom ) const override
        {
          if ( m_radius.get() > 0.0 ) {
            //TODO: Could improve this if geometries were able to return an
            //      inner bounding box/sphere.
            return true;
          } else {
            //special case r=0:
            return !geom.pointIsInside( m_center );
          }
        }

        void fillBasket( RNG& rng, NeutronBasket& nb ) override
        {
          nc_assert( !nb.full() );
          auto counts = m_stat.updateBasketCounts( nb );
          NCRYSTAL_DEBUGMMCMSG("Source Filling n="<<(counts.N-counts.i0));
          if ( counts.N == counts.i0 )
            return;
          nc_assert( counts.N > counts.i0 );
          auto& f = nb.fields;
          if ( m_radius.get() > 0.0 ) {
            //Randomised positions:
            double a, b;
            for ( std::size_t i = counts.i0; i < counts.N; ++i ) {
              std::tie(a,b) = randPointInUnitCircle(rng);
              auto p = m_center + m_a * a + m_b * b;
              f.x.data[i] = p[0];
              f.y.data[i] = p[1];
              f.z.data[i] = p[2];
            }
          } else {
            //special case r=0:
            for ( std::size_t i = counts.i0; i < counts.N; ++i )
              f.x.data[i] = m_center[0];
            for ( std::size_t i = counts.i0; i < counts.N; ++i )
              f.y.data[i] = m_center[1];
            for ( std::size_t i = counts.i0; i < counts.N; ++i )
              f.z.data[i] = m_center[2];
          }
          //Direction, weight, and energy are all constant:
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.ux.data[i] = m_dir[0];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.uy.data[i] = m_dir[1];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.uz.data[i] = m_dir[2];
          for ( std::size_t i = counts.i0; i < counts.N; ++i )
            f.w.data[i] = m_w;
          setEnergy( m_einfo, counts, rng, f.ekin );
        }
      };

      SourcePtr createSourceImpl( const StrView& raw_srcstr )
      {
        namespace PMC = parseMMCCfg;

        auto tokeninfo = PMC::tokenize( raw_srcstr );
        auto& tokens = tokeninfo.tokens;
        auto src_name =  tokeninfo.mainName;
        if ( !src_name.has_value() )
          NCRYSTAL_THROW2(BadInput,"Invalid src cfg: \""<<raw_srcstr<<"\"");

        //Energy is the same for all sources:
        auto energy = PMC::getValue_Energy( tokens,
                                            StrView("1.8") );//default 1.8Aa

        //Other common defaults:
        const char * common_defaults = "n=1e6";//;w=1.0";

        if ( src_name == "constant" ) {
          PMC::applyDefaults( tokens, common_defaults );
          PMC::applyDefaults( tokens, "x=0;y=0;z=0;ux=0;uy=0;uz=1;n=1e6" );
          PMC::checkNoUnknown(tokens,"ekin;wl;n;w;;x;y;z;ux;uy;uz","source");
          return makeSO<SourceConstant>
            ( PMC::getValue_sizet(tokens,"n"),
              std::move(energy),
              1.0,//PMC::getValue_weight( tokens, "w" ),
              Length{ PMC::getValue_dbl(tokens,"x") },
              Length{ PMC::getValue_dbl(tokens,"y") },
              Length{ PMC::getValue_dbl(tokens,"z") },
              NeutronDirection{ PMC::getValue_dbl(tokens,"ux"),
                                PMC::getValue_dbl(tokens,"uy"),
                                PMC::getValue_dbl(tokens,"uz") }
              );
        } else if ( src_name == "circular" ) {
          PMC::applyDefaults( tokens, common_defaults );
          PMC::applyDefaults( tokens, "x=0;y=0;z=0;ux=0;uy=0;uz=1;n=1e6" );
          PMC::checkNoUnknown(tokens,"ekin;wl;n;w;;x;y;z;ux;uy;uz;r","source");
          return makeSO<SourceUniformCircularBeam>
            ( PMC::getValue_sizet(tokens,"n"),
              std::move(energy),
              1.0,//PMC::getValue_weight( tokens, "w" ),
              Length{ PMC::getValue_dbl(tokens,"r") },
              Length{ PMC::getValue_dbl(tokens,"x") },
              Length{ PMC::getValue_dbl(tokens,"y") },
              Length{ PMC::getValue_dbl(tokens,"z") },
              NeutronDirection{ PMC::getValue_dbl(tokens,"ux"),
                                PMC::getValue_dbl(tokens,"uy"),
                                PMC::getValue_dbl(tokens,"uz") }
              );
        } else if ( src_name == "isotropic" ) {
          PMC::applyDefaults( tokens, common_defaults );
          PMC::applyDefaults( tokens, "x=0;y=0;z=0;r=0;n=1e6" );
          PMC::checkNoUnknown(tokens,"ekin;wl;n;w;;x;y;z;r","source");
          return makeSO<SourceIsotropic>
            ( PMC::getValue_sizet(tokens,"n"),
              std::move(energy),
              1.0,//PMC::getValue_weight( tokens, "w" ),
              Length{ PMC::getValue_dbl(tokens,"x") },
              Length{ PMC::getValue_dbl(tokens,"y") },
              Length{ PMC::getValue_dbl(tokens,"z") },
              PMC::getValue_dbl(tokens,"r")
              );
        } else {
          NCRYSTAL_THROW2(BadInput,"Unknown source type requested: \""<<src_name<<"\"");
        }
      }
    }
  }
}

NCMMC::SourcePtr NCMMC::createSource( const StrView& raw_srcstr )
{
  return createSourceImpl( raw_srcstr );
}

void NCMMC::sourceOptsDocsToJSON( std::ostream& os )
{
  os << "{\"intro_text\":";
  streamJSON(os, "Three MiniMC neutron source types are currently available,"
             " differing only in the way they generate the initial position"
             " and direction of the neutrons. Other parameters are common"
             " to all sources, including those related to neutron energies"
             " or how many neutrons to generate.");
  os << ",\"src_list\":[";
  using VV3 = SmallVector<std::array<StrView,3>,9>;
  unsigned sort_key = 0;
  auto addSrc = [&os,&sort_key]( StrView name,
                                 const VV3& specific_params,
                                 StrView descr )
  {
    os << "{\"name\":";
    streamJSON(os,name);
    os<< ",\"descr\":";
    streamJSON(os,descr);
    os<< ",\"specific_params\":";
    streamJSON(os,specific_params);
    os<< ",\"sort_key\":";
    streamJSON(os,sort_key++);
    os<<'}';
  };
  {
    VV3 v;
    v.push_back({"x","0","x coordinate of position [m]."});
    v.push_back({"y","0","y coordinate of position [m]."});
    v.push_back({"z","0","z coordinate of position [m]."});
    v.push_back({"ux","0","x coordinate of direction."});
    v.push_back({"uy","0","y coordinate of direction."});
    v.push_back({"uz","1","z coordinate of direction."});
    addSrc("constant",v,
           "A source which always provide the same position and"
           " direction of the initial neutron.");
  }
  os << ',';
  {
    VV3 v;
    v.push_back({"x","0","x coordinate of disk center [m]."});
    v.push_back({"y","0","y coordinate of disk center [m]."});
    v.push_back({"z","0","z coordinate of disk center [m]."});
    v.push_back({"ux","0","x coordinate of direction."});
    v.push_back({"uy","0","y coordinate of direction."});
    v.push_back({"uz","1","z coordinate of direction."});
    addSrc("circular",v,
           "A source which allows modelling of a beam of uniform"
           " circular cross section. Specifically, all neutrons are"
           " generated with a fixed direction, and with a position"
           " sampled uniformly within a disk of radius r. The"
           " center of the disk is configurable, while its normal"
           " is always parallel to the direction of the neutrons.");
  }
  os << ',';
  {
    VV3 v;
    v.push_back({"x","0","x coordinate of position [m]."});
    v.push_back({"y","0","y coordinate of position [m]."});
    v.push_back({"z","0","z coordinate of position [m]."});
    v.push_back({"r","0",
          "Signed radius of sphere. Neutrons are"
          " always first placed on the point (x,y,z)"
          " and assigned an isotropically sampled"
          " direction. If r!=0, the position"
          " is then translated a distance of -r"
          " along that direction."});
    addSrc("isotropic",v,
           "A source which emits neutrons isotropically from a"
           " particular point or spherical surface.");
  }
  os << "],\"commonpars_descr\":";
  streamJSON(os,
             //All except ekin/wl:
             "For all sources, the parameter \"n\" (default 1e6) is used"
             " to specify the number of neutrons to generate.");
             // " For completeness, the parameter \"w\" (default 1.0) sets the"
             // " initial weight of each neutron, but should in general not be"
             // " modified."
  os << ",\"commonpars_energy_descr\":";
  streamJSON(os,
             "Parameters \"ekin\" and \"wl\" are also common to all sources,"
             " and allow for a flexible specification of the generated"
             " neutron energies, as shown by these examples:");
  os << ",\"commonpars_energy_examples\":";
  {
    SmallVector<std::pair<StrView,StrView>,9> v;
    v.emplace_back("ekin=0.025",
                   "All neutrons are generated with an energy of 0.025eV.");
    v.emplace_back("ekin=0.025-0.05",
                   "Neutron energies are sampled uniformly in the"
                   " interval [0.025eV,0.050eV].");
    v.emplace_back("ekin=0.025+-0.001",
                   "Neutron energies are sampled from a log-normal"
                   " distribution with mean 0.025eV and RMS 0.001eV.");
    v.emplace_back("ekin=thermal:77",
                   "Neutron energies are sampled from a thermal"
                   " (Maxwell) distribution at 77K.");
    v.emplace_back("wl=1.8",
                   "All neutrons are generated with a wavelength of 1.8Aa.");
    v.emplace_back("wl=1.8-5",
                   "Neutron wavelengths are sampled uniformly in"
                   " the interval [1.8Aa,5Aa].");
    v.emplace_back("wl=1.8+-0.01",
                   "Neutron wavelengths are sampled from a log-normal"
                   " distribution with mean 1.8Aa and RMS 0.01Aa.");
    streamJSON(os,v);
  }
  os << ",\"commonpars_list\":[\"ekin\",\"wl\",\"n\"]";


  os << ",\"examples\":";
  {
    SmallVector<std::pair<StrView,StrView>,9> v;
    v.emplace_back("constant;wl=1.8;z=-0.1",
                   "1e6 neutrons are generated with a wavelength of 1.8Aa"
                   " as (x,y,z)=(0m,0m,-0.1m) with direction (0,0,1).");
    v.emplace_back("constant;ekin=0.001-0.1;x=2;z=-0.1;ux=0;uy=1;uz=1",
                   "1e6 neutrons are generated at (x,y,z)=(2m,0m,-0.1m) with"
                   " direction (0,1,1) and wavelengths are sampled"
                   " uniformly in the interval [0.001eV,0.1eV].");
    v.emplace_back("circular;wl=1.8;z=-0.1,r=0.1;n=1e8",
                   "1e8 neutrons are generated with a wavelength of 1.8Aa"
                   " with direction (0,0,1) and with a position sampled"
                   " uniformly on the disk with center (0m,0m,-0.1m) "
                   "and normal (0,0,1).");
    v.emplace_back("isotropic;ekin=0.025+-0.001;n=10000",
                   "10000 neutrons are generated at the point (0m,0m,0m) with an"
                   " isotropically sampled direction and an energy sampled"
                   " from a log-normal distribution with mean 0.025eV and"
                   " variance (0.001eV)^2.");
    streamJSON(os,v);
  }

  os << '}';
}
