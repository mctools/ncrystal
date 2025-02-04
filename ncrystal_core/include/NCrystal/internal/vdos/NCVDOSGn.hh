#ifndef NCrystal_VDOSGn_hh
#define NCrystal_VDOSGn_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  class VDOSEval;

  class VDOSGn final : private MoveOnly {

  public:
    //Construct Sjolanders *asymmetric* Gn functions for a given VDOS
    //(Vibrational Density Of States) phonon spectrum as defined in Sjolander,
    //Arkiv for Fysik., Bd 14, nr 21, 1958 (chapter 2.). This facilitates the
    //attainment of higher order phonon contributions to incoherent inelastic
    //neutron scatterings, and relies on an iterative procedure (see Sjolander's
    //eq. II.27) in which the first (single-phonon) tabulated G1 function -
    //provided by an VDOSEval instance - can essentially be self-convolved n
    //times in order to obtain the Gn function with information about n'th order
    //phonon contributions.
    //
    //The actual technical implementation has been tuned to be more efficient
    //and accurate than the direct application of Sjolander's model would
    //be. For reasons numerical stability, the G{2n} and G{2n+1} functions are
    //respectively obtained by convolutions G{n}*G{n} and G{n}*G{n+1}, so a
    //given Gn function is obtained by just ~log2(n) convolutions. This
    //significantly reduces accumulation of numerical errors. The convolutions
    //themselves are carried out by a Fast-Fourier-Transform method with fixed
    //binning, and to enable efficient calculations of *very* high orders (up to
    //at least 10000 orders have been tested), the binning and range of high
    //order spectra are continuously reduced through a customized truncation &
    //thinning procedure which has tremendous impact on the computational
    //efficiency in terms of CPU and memory usage, but essentially no impact on
    //the validity of the results (at least when using the default
    //truncation/thinning options below).
    //
    //Upon construction, the VDOSGn class only contains G1, but calling code can
    //call the growMaxOrder function to dynamically expand this to higher
    //orders. It is an error to query results of Gn functions unless
    //growMaxOrder has first been called with at least n.
    //
    //Note that for consistency between the present interface and Sjolander's
    //notation, we index the orders, n=1,2,3,4,... Although this is contrary to
    //the 0-based indexing normally used in C++, it would be highly confusing to
    //use n=0 to indicate single-phonon scattering, not the least because it is
    //customary in neutron scattering to occasionally identify elastic
    //scattering with "0-order phonon scattering". For added safety, the Order
    //type used in the VDOSGn interface below is a special class which provides
    //extra sanity checks in debug builds of provided numbers.

    class Order;

    ///////////////////
    // Constructor: //
    ///////////////////

    //Initialise based on VDOS (in form of VDOSEval instance) and optionally
    //choices for truncation/thinning.
    enum class TruncAndThinningChoices { Default, Disabled };
    struct TruncAndThinningParams {
      int minOrder = 5;//Below this order, no truncation or thinning takes place
                       //(0=always, -1=never)
      unsigned thinNBins = 1000;//double binwidth whenever number of bins
                                //exceeds this value (0 disables)
      double truncationThreshold = 1e-14;//trim ranges to remove negligible
                                         //noise at edges (0 disables)
      TruncAndThinningParams(TruncAndThinningChoices);
      TruncAndThinningParams() = default;
    };

    VDOSGn( const VDOSEval&,
            const TruncAndThinningParams ttpars = TruncAndThinningChoices::Default );

    /////////////////
    // Destructor: //
    /////////////////
    ~VDOSGn();

    ///////////////////////////
    // Material temperature: //
    ///////////////////////////
    double kT() const { return m_kT; }

    ////////////////////////////////////////////
    //Check or increase maximum available Gn: //
    ////////////////////////////////////////////
    void growMaxOrder( Order );//Not MT-safe
    Order maxOrder() const;

    //////////////////////////////
    // Access properties of Gn: //
    //////////////////////////////

    //Evaluate Gn, for n in 1..maxOrder, at given energy point:
    double eval( Order n, double energy ) const;
    //get energy range where spectrum is above relthreshold of the maximal value:
    PairDD eRange( Order n, double relthreshold ) const;
    //Full energy range, bin width and spectrum:
    PairDD eRange( Order n) const;
    double binWidth( Order) const;
    const VectD& getRawSpectrum( Order ) const;

    ///////////////////////////////////////////////////////////
    // Enable verbose output (default is disabled unless the //
    // NCRYSTAL_DEBUG_PHONON environment variable is set.    //
    ///////////////////////////////////////////////////////////

    static void enableVerboseOutput(bool status = true);
    static bool verboseOutputEnabled();

  private:
    struct Impl;
    Pimpl<Impl> m_impl;
    double m_kT;
  };

  class VDOSGn::Order {
    //Essentially an unsigned integer which (in debug builds) guards against
    //assignment from negative numbers, 0, or numbers too high to be an actual
    //order number.
  public:
    Order(unsigned o) : m_order(o) { checkValid(); }
    Order(int o) : m_order(o) { nc_assert(o>0); checkValid(); }
    Order(const Order& o) = default;
    Order& operator=(const Order& o) { m_order = o.m_order; return *this; }
    Order& operator=(const int& o) { nc_assert(o>0); m_order = o; checkValid(); return *this; }
    Order& operator=(const unsigned& o) { m_order = o; checkValid(); return *this; }
    bool operator<(const Order& o) { return m_order < o.m_order; }
    bool operator<=(const Order& o) { return m_order <= o.m_order; }
    Order& operator++() { ++m_order; checkValid(); return *this; }
    unsigned value() const { return m_order; }
  private:
    unsigned m_order;
    void checkValid() { nc_assert( m_order>0 && m_order<100000 ); }
  };

}

#endif
