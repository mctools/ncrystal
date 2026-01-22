#ifndef NCrystal_ExtnHelper_hh
#define NCrystal_ExtnHelper_hh

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

#include "NCrystal/internal/utils/NCExtraTypes.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //
  // Templated helper class for modelling Bragg diffraction based on a
  // particular extinction model.
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace Extn {

    template<class TModel>
    class ExtnHelper final : MoveOnly {
    public:
      using ModelData = typename TModel::ModelData;
      using PlaneData = typename TModel::PlaneData;

      ExtnHelper( ModelData&& md,
                  const PowderBraggInput::CellData& cell,
                  const std::vector<PowderBraggInput::Plane>& planes )
        : m_modelData( std::move(md) )
      {
        const double v0_times_natoms = cell.volume * cell.n_atoms;
        if ( !std::isfinite(v0_times_natoms) || !(v0_times_natoms>0) )
          NCRYSTAL_THROW(BadInput,"v0_times_natoms is not a positive number.");//fixme: copied from NCPowderBragg.cc

        m_xsectfact = 0.5/v0_times_natoms;
        m_xsectfact *= wl2ekin(1.0);//Adjust units so we can get cross sections
                                    //through multiplication with 1/ekin instead
                                    //of wl^2.

        if ( !planes.empty() ) {
          m_planes.reserve( planes.size() );
          for ( auto& p : planes )
            m_planes.emplace_back( TModel::initPlaneData( p ) );
          m_threshold = NeutronEnergy{ wl2ekin(2.0*planes.front().dsp) };
        }
      }

      EnergyDomain domain() const noexcept
      {
        return { m_threshold, NeutronEnergy{kInfinity} };
      }

      CrossSect xsect_NoCache( NeutronEnergy ekin ) const
      {
        if ( ekin < m_threshold || ncisinf(ekin.get()) )
          return CrossSect{ 0.0 };
        auto neutronData = TModel::initNeutronData( m_modelData, ekin );
        auto wlhalf = ekin.wavelength().get() * 0.5;
        StableSum contrib;
        for ( auto& plane : m_planes ) {
          if ( plane.dsp < wlhalf )
            break;
          const double E = TModel::extinctionFactor( m_modelData,
                                                     neutronData,
                                                     plane );
          contrib.add( plane.fdm * ncclamp( E, 0.0, 1.0 ) );
        }
        return CrossSect{ contrib.sum() * m_xsectfact / ekin.get() };
      }

      //For use in caching or sampling we always need a vector of commulative
      //contributions:
      void updatePlaneCommulContribs( NeutronEnergy ekin,
                                      std::vector<double>& out ) const
      {
        if ( ekin < m_threshold || ncisinf(ekin.get()) ) {
          out.clear();
          return;
        }
        std::vector<double> tmp;
        std::swap( out, tmp );
        auto neutronData = TModel::initNeutronData( m_modelData, ekin );
        auto wlhalf = ekin.wavelength().get() * 0.5;
        StableSum contrib;
        for ( auto& plane : m_planes ) {
          if ( plane.dsp < wlhalf )
            break;
          const double E = TModel::extinctionFactor( m_modelData,
                                                     neutronData,
                                                     plane );
          contrib.add( plane.fdm * ncclamp( E, 0.0, 1.0 ) );
          tmp.push_back( contrib.sum() );
        }
        std::swap( out, tmp );
      }

      CrossSect xsect( NeutronEnergy ekin,
                       const std::vector<double>& contribs ) const
      {
         return CrossSect{ contribs.empty()
                          ? 0.0
                          : ( contribs.back() * m_xsectfact / ekin.get() ) };
      }

      CosineScatAngle sampleScatMu( RNG& rng,
                                    NeutronEnergy ekin,
                                    const std::vector<double>& contribs ) const
      {
        if ( contribs.empty() )
          return CosineScatAngle{1.0};//do not scatter

        //Sample plane by contribution:
        auto idx = pickRandIdxByWeight( rng.generate(), contribs );

        //Look up d-spacing and convert to cos(scat_angle) via Bragg condition:

        const double sin_theta_bragg_squared
          = wl2ekin(2.0*vectAt(m_planes,idx).dsp) / ekin.get();

        //scatter angle A=2*theta_bragg, so with x=sin^2(theta_bragg), we have:
        //   x = sin^2(A/2)= (1-cosA)/2 => 1-2x = cosA = mu

        const double mu = 1.0 - 2.0 * sin_theta_bragg_squared;
        nc_assert_always(ncabs(mu)<=1.0);//fixme: _always
        return CosineScatAngle{ mu };
      }

      CosineScatAngle sampleScatMu_NoCache( RNG& rng, NeutronEnergy ekin ) const
      {
        if ( ekin < m_threshold || ncisinf(ekin.get()) )
          return CosineScatAngle{1.0};
        //Must do it the hard way and use temporary vector.
        std::vector<double> contribs;
        contribs.reserve(1024);
        this->updatePlaneCommulContribs( ekin, contribs );
        return this->sampleScatMu( rng, ekin, contribs );
      }

    private:
      ModelData m_modelData;
      std::vector<PlaneData> m_planes;
      NeutronEnergy m_threshold = NeutronEnergy{kInfinity};
      double m_xsectfact = 0.0;
    };

  }

}

#endif
