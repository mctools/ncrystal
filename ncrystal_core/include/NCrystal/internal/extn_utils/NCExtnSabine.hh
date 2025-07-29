#ifndef NCrystal_ExtnSabineMdl_hh
#define NCrystal_ExtnSabineMdl_hh

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

//#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCExtraTypes.hh"
#include "NCrystal/internal/extn_utils/NCExtnUtils.hh"
//#include "NCrystal/internal/extn_utils/NCExtnEval.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //
  // Models implementing models described in:
  //
  // Sabine, T. M. "The flow of radiation in a real crystal." (2006): 609-616.
  // doi: 10.1107/97809553602060000603
  //
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace Extn {


    struct SabineMdlPurePrimary final {//fixme namespace instead?

      struct ModelData {
        double squared_blksize_div_v0;//units ???
      };

      static ModelData initModelData( const PowderBraggInput::CellData& cell,
                                      Length blockSize )
      {
        ModelData res;
        double blksz_Aa = blockSize.get() / Length::angstrom;
        res.squared_blksize_div_v0 = ncsquare( blksz_Aa / cell.volume ) * 1e-8;
        return res;
      }

      struct NeutronData {
        double wlsq;
        double factor_fsq2x;
      };

      static NeutronData initNeutronData( const ModelData& md,
                                          NeutronEnergy ekin )
      {
        NeutronData res;
        res.wlsq = ekin2wlsq( ekin.get() );
        res.factor_fsq2x = res.wlsq * md.squared_blksize_div_v0;
        return res;
      }

      struct PlaneData {
        //Fields required by all models:
        double dsp;
        double fdm;
        //Other fields, as needed by this particular model:
        double fsq;//structure factor squared [barn]
        double inv2dsp_squared; // = (2*dspacing)^(-2)
      };

      static PlaneData initPlaneData( const PowderBraggInput::Plane& p )
      {
      // double dsp ;//d-spacing [Angstrom]
      // double fsq;//structure factor squared [barn]
      // double mult;//multiplicity [integral, in double for efficiency]
        PlaneData res;
        res.dsp = p.dsp;
        res.fdm = p.fsq * p.dsp * p.mult;
        res.fsq = p.fsq;
        res.inv2dsp_squared = ncsquare( 1.0 / ( 2.0 * p.dsp ) );
        return res;
      }

      static double extinctionFactor( const ModelData& /*used through neutron data*/,
                                      const NeutronData& n,
                                      const PlaneData& p )
      {
        //DblBasket sabine_x, El, Eb, sinth_sq, costh_sq;
        const double sabine_x = n.factor_fsq2x * p.fsq;
        const double El = calcSabineEl_y0( sabine_x );
        const double Eb = calcSabineEb_y0( sabine_x );
        double sinth_sq = n.wlsq * p.inv2dsp_squared;
        //Fixme: just do max(0.0,...) on final E value just before returning?
        //Or have the framework clip the values to [0,1] for absolute certainty!
        sinth_sq = std::min<double>(1.0,sinth_sq);
        double costh_sq = 1.0 - sinth_sq;
        costh_sq = std::max<double>(costh_sq,0.0);
        return El*costh_sq + Eb*sinth_sq;
      }
    };
  }

  //fixme: cleanup:
    // // Vect(fsq), Vect(dsp) (mult can be handled by caller)

    // struct DblBasket {
    //   //Data structure optimized for SIMD.
    //   static constexpr unsigned count = 16;
    //   align(alignof(double)*count) double data[count];
    // };

    // struct FDMBasket {
    //   using DblBasket::count;
    //   DblBasket fsq;
    //   DblBasket inv2dsp;// = 1/(2*dspacing)
    //   DblBasket fdm;// = fsq * dsp * mult
    // };
    // std::vector<FDMBasket> planes;
    //   // std::vector<DblBasket> fsq;
    //   // std::vector<DblBasket> dsp;
    //   // std::vector<DblBasket> mult;
    //   //+ invdsp ?  (1/2d)^2... or 1/2d
    //   //+ d*f*m (instead of mult)

    // //};

    // std::vector<FDMBasket> vectorizeFDM( const PowderBraggInput::Data& data )
    // {
    //   std::vector<FDMBasket> res;
    //   unsigned i = FDMBasket::count;
    //   for ( auto& e : data.planes ) {
    //     if ( i >= FDMBasket::count ) {
    //       res.emplace_back();
    //       i = 0;
    //     }
    //     auto& b = res.back();
    //     b.fsq[i] = e.fsq;
    //     b.inv2dsp[i] = 1.0 / (2.0 * e.dsp);
    //     b.fdm[i] = e.fsq * e.dsp * e.mult;
    //     ++i;
    //   }
    //   //Top off last basket with dummy entries having Fsq=0:
    //   while ( i < FDMBasket::count ) {
    //     auto& b = res.back();
    //     b.fsq[i] = 0.0;
    //     b.inv2dsp[i] = 1.0;
    //     b.fdm[i] = 0.0;
    //     ++i;
    //   }
      // }

//       template<class TDblBasket, class TFDMBasket>
//       static TDblBasket extinctionFactor( const ModelData& md,
//                                           const NeutronData& nd,
//                                           const FDMBasket& fdm )
//       {
//         constexpr auto N = DblBasket::count;
//         DblBasket sabine_x, El, Eb, sinth_sq, costh_sq;
//         for ( auto i : ncrange( N ) )
//           sabine_x[i] = nd.factor_fsq2x * fdm.fsq[i];

//         for ( auto i : ncrange( N ) )
//           El[i] = calcSabineEl_y0( sabine_x[i] );

//         for ( auto i : ncrange( N ) )
//           Eb[i] = calcSabineEb_y0( sabine_x[i] );

//         for ( auto i : ncrange( N ) )
//           sinth_sq[i] = nd.wlsq * ncsquare( fdm.inv2dsp[i] );//fixme: cache inv2dspsq ?

//         //Fixme: just do max(0.0,...) on final E value just before returning?
//         //Or have the framework clip the values to [0,1] for absolute certainty!
//         for ( auto i : ncrange( N ) )
//           sinth_sq[i] = std::min<double>(1.0,sinth_sq[i]);

//         for ( auto i : ncrange( N ) )
//           costh_sq[i] = 1.0 - sinth_sq[i];

//         for ( auto i : ncrange( N ) )
//           costh_sq[i] = std::max<double>(costh_sq[i],0.0);

//         //E = El*cos^2(th) + Eb*sin^2(th):
//         for ( auto i : ncrange( N ) )
//           El[i] *= costh_sq[i];
//         for ( auto i : ncrange( N ) )
//           Eb[i] *= sinth_sq[i];
//         for ( auto i : ncrange( N ) )
//           El[i] += Eb[i];
//         return El;
//       }


//       template<TAccumulator>
//       void accumulatePlaneContribs( NeutronEnergy, TAccumulator& acc )
//       {
//         if ( ekin < m_threshold )
//           return;
//         const auto wl = ekin.wavelength();
//         const double v0 = m_data.cell.volume;
//         const double kkk = ncsquare(wl.get()*m_domainSizeAa/v0)*1e-8;



//         const double wlhalf = wl.get()*0.5;
//         const double factor
//           = ncsquare(wl.get())/(2.0*m_data.cell.volume*m_data.cell.n_atoms);


//         StableSum contrib;
//         for ( auto& e : m_data.planes ) {
//           if ( e.dsp < wlhalf )
//             break;



//           double sabine_x = kkk * e.fsq;
//           double El = calcSabineEl_y0( sabine_x );
//           double Eb = calcSabineEb_y0( sabine_x );
//           double sinth_sq = std::min<double>(1.0,ncsquare(0.5 * wl.get() / e.dsp));
//           double costh_sq = std::max<double>(1.0 - sinth_sq,0.0);
//           double extinction_correction = El * costh_sq + Eb * sinth_sq;
//           contrib.add( e.dsp * e.fsq * e.mult * extinction_correction  );
//         }
//         return CrossSect{ factor * contrib.sum()  };


//       }

// std::vector<double>&

// // Length blockSize );
// //     };
//     };
//   }
}

#endif
