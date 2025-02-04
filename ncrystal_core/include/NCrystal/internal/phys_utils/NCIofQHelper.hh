#ifndef NCrystal_IofQHelper_hh
#define NCrystal_IofQHelper_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCPointwiseDist.hh"

namespace NCRYSTAL_NAMESPACE {

  class IofQHelper {
  public:

    ///////////////////////////////////////////////////////////////////////////////
    //                                                                           //
    // Helper class which can be used to implement elastic physics described via //
    // an I(Q) curve. The class takes a user provided I(Q) curve and makes it    //
    // possible to either evaluate the integral of or sample according to, the   //
    // function Q*I(Q), over the interval evaluate the integral from 0 to        //
    // 2k. Here k is the wavenumber of a neutron of some particular energy.      //
    //                                                                           //
    // The I(Q) curve must be defined by a number of (Q,I(Q)) points, which      //
    // internally are translated to (Q,Q*I(Q)) points that are then connected    //
    // linearly to form a complete curve. If the Q value of the first provided   //
    // point is not 0.0, a point at Q=0.0 is added, so that the Q*I(Q) curve is  //
    // flat up until the first actual provided point. This choice of             //
    // extrapolation behaviour towards Q=0 is chosen because it ensures that for //
    // small energies, the integrand becomes a constant, which again ensures     //
    // that cross sections will ultimately have a standard 1/sqrt(E) (aka "1/v") //
    // behaviour at low energies. Finally, the I(Q) curve is assumed to be       //
    // vanishing at Q values above the last provided point. If a particular      //
    // use-case necessitates a different extrapolation behaviour, the caller can //
    // simply add additional points at either end of the curve as needed.        //
    //                                                                           //
    ///////////////////////////////////////////////////////////////////////////////

    //The class is constructed from Q and I(Q) values:

    IofQHelper( const VectD& Q, const VectD& IofQ );
    IofQHelper( const std::pair<VectD,VectD>& Q_and_IofQ );

    //Calculate the integral of Q*I(Q) from Q=0 to Qmax=2k, where k is the
    //wavenumber of the neutron of the provided energy. Note that to convert it
    //to a cross section one must still multiply it with a factor of c/E where c
    //is an appropriate constant and E is the neutron energy:
    double calcQIofQIntegral( NeutronEnergy ) const;

    //Sample a Q value according to Q*I(Q) over the interval from Q=0 to
    //Qmax=2k, where k is the wavenumber of the neutron of the provided energy:
    double sampleQValue( RNG&, NeutronEnergy ) const;

    //QMax value:
    double getQMax() const;

  private:
    PointwiseDist m_pwdist;
    NeutronEnergy m_ekinMax;
    double m_normFact;
    struct internal_t;
    IofQHelper( internal_t );
  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::IofQHelper::IofQHelper( const std::pair<VectD,VectD>& QI )
  : IofQHelper(QI.first,QI.second)
{
}

inline double NCrystal::IofQHelper::calcQIofQIntegral( NeutronEnergy ekin ) const
{
  if ( ekin >= m_ekinMax )
    return m_normFact;
  constexpr double kkk = 4.0 * ekin2ksq(1.0);
  const double twok = std::sqrt( kkk * ekin.dbl() );
  return m_pwdist.commulIntegral( twok ) * m_normFact;
}

inline double NCrystal::IofQHelper::sampleQValue( RNG& rng, NeutronEnergy ekin ) const
{
  constexpr double kkk = 4.0 * ekin2ksq(1.0);
  const double twok = std::sqrt( kkk * std::min<double>(m_ekinMax.dbl(),ekin.dbl()) );
  return m_pwdist.sampleBelow( rng, twok );
}

#endif
