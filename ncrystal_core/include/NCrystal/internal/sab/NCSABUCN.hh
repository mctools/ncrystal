#ifndef NCrystal_SABUCN_hh
#define NCrystal_SABUCN_hh

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

#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCSpline.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace UCN {

    //Utilities and infrastructure needed for modelling UCN (ultra cold neutron)
    //production based on a given S(alpha,beta) scattering kernel. This includes
    //a utility class (UCNHElper) for implementing the actual model, as well as
    //two Scatter class implementations which can be used to inject this model
    //into NCrystal as needed.

    class UCNHelper final {

      //////////////////////////////////////////////////////////////////////////////
      //                                                                          //
      // Helper class which analyses a given S(alpha,beta) kernel and provides a  //
      // dedicated model to directly model UCN (ultra cold neutron) production    //
      // events. In addition to a kernel, it needs a ucn_threshold, which is the  //
      // energy below which neutrons are considered UCN.                          //
      //                                                                          //
      // More specifically, at initialisation it carefully integrates the kernel  //
      // over the relevant phase space (i.e. within the usual kinematic boundary  //
      // curves and for beta values which leave the neutron with                  //
      // E<=ucn_threshold). This provides a cross section curve as a function of  //
      // E, which is simply used for fast interpolation to later provide cross    //
      // sections. Similarly, an "S_overlay(E)" curve is also prepared at         //
      // initialisation time, which makes it possible to provide fast-but-exact   //
      // sampling of scattering outcomes.                                         //
      //                                                                          //
      // To avoid various numerical issues, the UCN process will always have 0    //
      // cross sections for neutron energies less than ucn_threshold.             //
      //                                                                          //
      // Note that the model implemented is mathematically and numerically able   //
      // to work with very high ucn_threshold values, also much higher than what  //
      // is usually considered "ultra cold". The only caveat is that with higher  //
      // ucn_threshold values, the computational efficiency for sampling will     //
      // degrade. It might nonetheless be useful to use such a high ucn_threshold //
      // value forvalidation studies, due to the artifact-free model implemented  //
      // (i.e. it can be used as a good benchmark reference).                     //
      //                                                                          //
      //////////////////////////////////////////////////////////////////////////////

    public:

      UCNHelper( shared_obj<const SABData>, NeutronEnergy ucn_threshold );

      NeutronEnergy ucnThreshold() const { return m_data.ucnthr; }

      //Cross section, sampling, and domain::
      CrossSect crossSection( NeutronEnergy ) const;
      EnergyDomain domain() const noexcept;
      ScatterOutcomeIsotropic sampleScatterIsotropic( RNG&, NeutronEnergy ) const;

      //Check if cross section is vanishing everywhere (same as
      //domain().isNull() but more efficient):
      bool isNull() const;

      //Direct access to XS curve (nb: not valid if domain().isNull()):
      const PiecewiseLinearFct1D& accessXSCurve() const { return m_data.xs; }


      //Internal fct, exposed for unit tests:
      static std::pair<std::vector<StableSum>,VectD> getSIntegralsAndOverlayVals( const SABData& sabData,
                                                                                  NeutronEnergy ucn_threshold,
                                                                                  const VectD& eVals );

      double worstSamplingAcceptanceRate() const { return m_data.worst_AR; }
      double averageSamplingAcceptanceRate() const { return m_data.average_AR; }
    private:
      struct Data {
        PiecewiseLinearFct1D xs, overlay;//xs(E) and Soverlay(E)
        NeutronEnergy ucnthr;
        shared_obj<const SABData> m_sabData;
        double worst_AR, average_AR;
      };
      Data m_data;
    };

    class UCNScatter final : public ProcImpl::ScatterIsotropicMat {
      ////////////////////////////////////////////
      // Scatter class which wraps a UCNHelper. //
      ////////////////////////////////////////////
    public:

      //Construct directly:
      UCNScatter( shared_obj<const SABData>, NeutronEnergy ucn_threshold );

      //Or create shared object with cache:
      static shared_obj<const UCNScatter> createWithCache( shared_obj<const SABData>,
                                                           NeutronEnergy ucn_threshold);

      //Required methods + UCNHelper access:
      const char * name() const noexcept override { return "UCNScatter"; }
      CrossSect crossSectionIsotropic( CachePtr&, NeutronEnergy ) const override;
      ScatterOutcomeIsotropic sampleScatterIsotropic( CachePtr&, RNG&, NeutronEnergy ) const override;
      EnergyDomain domain() const noexcept override;
      const UCNHelper& ucnHelper() const { return m_helper; }

    protected:
      Optional<std::string> specificJSONDescription() const override;
    private:
      UCNHelper m_helper;
    };

    class ExcludeUCNScatter final : public ProcImpl::ScatterIsotropicMat {
    public:
      ////////////////////////////////////////////////////////////////////////////
      // Scatter composition class which essentially "subtracts" a UCNScatter   //
      // (B) from any other (isotropic) scatter (A). It does this by directly   //
      // subtracting the cross sections A-B (returning 0 barn if it would lead  //
      // to negative values), and for sampling it simply keeps sampling A until //
      // the result does not have a lower energy than that given by the         //
      // ucnThreshold() of B.                                                   //
      //                                                                        //
      // It's intended usage is to provide usual inelastic models *except* UCN  //
      // production. This thus allows UCN processes to be split out from the    //
      // rest, which is useful for biasing strategies.                          //
      ////////////////////////////////////////////////////////////////////////////
      ExcludeUCNScatter( shared_obj<const ProcImpl::ScatterIsotropicMat>,
                         shared_obj<const UCNScatter> );
      const char * name() const noexcept override { return "ExcludeUCNScatter"; }
      CrossSect crossSectionIsotropic( CachePtr&, NeutronEnergy ) const override;
      ScatterOutcomeIsotropic sampleScatterIsotropic( CachePtr&, RNG&, NeutronEnergy ) const override;
    protected:
      Optional<std::string> specificJSONDescription() const override;
    private:
      shared_obj<const ProcImpl::ScatterIsotropicMat> m_wrappedScatter;
      shared_obj<const UCNScatter> m_ucnScatter;
      EnergyDomain m_ucnDomain;
    };

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::CrossSect NCrystal::UCN::UCNHelper::crossSection( NeutronEnergy ekin ) const
{
  if ( isNull() )
    return CrossSect{ 0.0 };
  const double emin = m_data.xs.xValues().front();
  if ( ekin.dbl() < emin ) {
    auto uf_val = m_data.xs.outOfBoundsYValues().underflowYValue;
    if ( uf_val.has_value() ) {
      nc_assert(uf_val.value()==0.0);
      return CrossSect{ 0.0 };
    }
    //Extrapolate by 1/sqrt(E):
    return CrossSect{ m_data.xs.yValues().front() * std::sqrt( emin / ekin.dbl() ) };
  } else {
    //Standard:
    return CrossSect{ m_data.xs( ekin.dbl() ) };
  }
}

inline bool NCrystal::UCN::UCNHelper::isNull() const
{
  return m_data.xs.yValues().size()==2 && !m_data.xs.yValues().front() && !m_data.xs.yValues().back();
}

#endif
