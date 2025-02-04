#ifndef NCrystal_Hists_hh
#define NCrystal_Hists_hh

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
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace Hists {

    //Utility classes for implementing histograms, both the bin contents, the
    //running statistics, and a Hist1D class combining both:

    enum class AllowWeights { YES, NO, YES_BUT_DELAY_ALLOC };
    enum class OverflowHandling { Clamp, Record };

    namespace detail {
      struct EmptyVect {
        void resize(std::size_t, double) {};
        constexpr bool empty() const { return true; }
        constexpr double operator[](std::size_t) const { return 0.0; }
        ncconstexpr17 EmptyVect& operator=( const EmptyVect&) { return *this; }
      };
      struct EmptyPairDD {
        void swap( PairDD& ) {}
        PairDD operator()() const { return PairDD{0.0,0.0}; }
      };
    }

    ///////////////////////////////////////////////////////////////////////////
    // HistBinData1D class                                                   //
    ///////////////////////////////////////////////////////////////////////////

    template<
      AllowWeights ALLOW_WEIGHTS = AllowWeights::YES,
      OverflowHandling OF_HANDLING = OverflowHandling::Record,
      class TStorage = std::vector<double>
      >
    class HistBinData1D final {
    public:
      using size_t = std::uint_least32_t;
      constexpr static AllowWeights opt_allow_weights = ALLOW_WEIGHTS;
      constexpr static OverflowHandling opt_of_handling = OF_HANDLING;
      using opt_storate_type = TStorage;
    private:
      static constexpr size_t
      nbins_overflow = ( opt_of_handling == OverflowHandling::Record ? 2 : 0 );
    public:

      HistBinData1D( size_t nbins, double xmin, double xmax)
        : m_xmin(xmin),
          m_xmax(xmax),
          m_nbins(nbins)
      {
        nc_assert_always( std::isfinite( m_xmin ) );
        nc_assert_always( std::isfinite( m_xmax ) );
        nc_assert_always( m_xmax > m_xmin );
        nc_assert_always(nbins>=1 && nbins < 1000000000 );
        m_content.resize(m_nbins + nbins_overflow,0.0);
        if ( opt_allow_weights == AllowWeights::YES )
          m_errors.resize(m_nbins + nbins_overflow,0.0);
        const double delta = (m_xmax-m_xmin) / m_nbins;
        nc_assert_always(delta>0.0);
        m_invDelta = 1.0 / delta;
        if ( nbins_overflow == 0 ) {
          PairDD cr( 0.99*delta, m_xmax - (0.99*delta+m_xmin) );
          m_clampRange.swap(cr);
        }
      }

      void toJSON( std::ostringstream& os, bool include_data_arrays ) const
      {
        streamJSONDictEntry(os, "xmin", m_xmin,
                            JSONDictPos::FIRST);
        streamJSONDictEntry(os, "xmax", m_xmax);
        streamJSONDictEntry(os, "nbins", m_nbins);

        if ( include_data_arrays ) {
          streamJSONDictEntry(os, "content", this->getContents());
          if ( m_errors.empty() ) {
            VectD dummy;
            streamJSONDictEntry(os, "errorsq", dummy);
          } else {
            streamJSONDictEntry(os, "errorsq", this->getErrorsSquared());
          }
        }

        if ( opt_of_handling == OverflowHandling::Record ) {
          nc_assert_always( nbins_overflow == 2 );
          streamJSONDictEntry(os, "underflow", vectAt(m_content,0) );
          streamJSONDictEntry(os, "overflow", vectAt(m_content,m_nbins+1));
          if ( opt_allow_weights == AllowWeights::NO || m_errors.empty() ) {
            streamJSONDictEntry(os, "underflow_errorsq", vectAt(m_content,0) );
            streamJSONDictEntry(os, "overflow_errorsq", vectAt(m_content,m_nbins+1));
          } else {
            nc_assert_always(m_errors.size()==m_nbins+2);
            streamJSONDictEntry(os, "underflow_errorsq", vectAt(m_errors,0) );
            streamJSONDictEntry(os, "overflow_errorsq", vectAt(m_errors,m_nbins+1));
          }
        }
        os << '}';
      }

      bool mergeCompatible( const HistBinData1D& o ) const
      {
        return ( m_xmin == o.m_xmin
                 && m_xmax == o.m_xmax
                 && m_nbins == o.m_nbins );
      }

      void merge( const HistBinData1D& o )
      {
        if (!mergeCompatible(o))
          NCRYSTAL_THROW(CalcError,"Attempting to merge incompatible 1D"
                         " histogram data");
        //contents:
        {
          nc_assert_always(m_content.size() == o.m_content.size());
          double * it = m_content.data();
          double * itE = it + m_content.size();
          const double * itO = o.m_content.data();
          while (it!=itE)
            *(it++) += *(itO++);
        }

        //errors:
        if ( opt_allow_weights == AllowWeights::NO
             || o.m_errors.empty() ) {
          return;
        }

        if ( opt_allow_weights == AllowWeights::YES_BUT_DELAY_ALLOC
             && m_errors.empty() && o.m_errors.empty() ) {
          return;
        }

        if ( opt_allow_weights == AllowWeights::YES_BUT_DELAY_ALLOC
             && m_errors.empty() ) {
          m_errors = m_content;
        }

        {
          nc_assert_always(m_errors.size() == o.m_errors.size());
          double * it = m_errors.data();
          double * itE = it + m_errors.size();
          const double * itO = o.m_errors.data();
          while (it!=itE)
            *(it++) += *(itO++);
        }
      }

      double getBinContent( size_t ibin) const
      {
        return vectAt(m_content,(nbins_overflow==0?ibin:ibin+1));
      }

      double getBinError( size_t ibin) const
      {
        constexpr size_t offset = ( nbins_overflow==0 ? 0 : 1 );
        if ( opt_allow_weights == AllowWeights::YES )
          return std::sqrt( vectAt(m_errors,ibin+offset) );
        if ( opt_allow_weights == AllowWeights::NO )
          return std::sqrt( vectAt(m_content,ibin+offset) );
        return ( m_errors.empty()
                 ? std::sqrt( vectAt(m_content,ibin+offset) )
                 : std::sqrt( vectAt(m_errors,ibin+offset) ) );
      }

      double getXMin() const { return m_xmin; }
      double getXMax() const { return m_xmax; }
      unsigned getNBins() const { return m_nbins; }

      template<OverflowHandling U = opt_of_handling>
      double getUnderflow( typename std::enable_if<
                           U==OverflowHandling::Record>::type* = nullptr ) const
      {
        return vectAt(m_content,0);
      }

      template<OverflowHandling U = opt_of_handling>
      double getOverflow( typename std::enable_if<
                          U==OverflowHandling::Record>::type* = nullptr ) const
      {
        return vectAt(m_content,m_nbins+1);
      }

      //TODO: Provide access to uncertainty in underflow/overflow counts!?

      Span<const double> getContents() const
      {
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        return ( nbins_overflow==0
                 ? Span<const double>{ m_content.data(),
                                       m_content.data() + m_nbins }
                 : Span<const double>{ m_content.data() + 1,
                     m_content.data() + m_nbins + 2 } );
      }

      Span<const double> getErrorsSquared() const
      {
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        if ( opt_allow_weights == AllowWeights::NO || m_errors.empty() )
          return getContents();
        return ( nbins_overflow==0
                 ? Span<const double>{ m_errors.data(),
                                       m_errors.data() + m_nbins }
                 : Span<const double>{ m_errors.data() + 1,
                     m_errors.data() + m_nbins + 2 } );
      }

      void fill(double val) {
        size_t idx = this->valueToBin( val );
        vectAt(m_content,idx) += 1.0;
        if ( opt_allow_weights == AllowWeights::YES)
          vectAt(m_errors,idx) += 1.0;
        if ( opt_allow_weights == AllowWeights::YES_BUT_DELAY_ALLOC
             && !m_errors.empty() )
          vectAt(m_errors,idx) += 1.0;
      }

      void fillN(double val, size_t N )
      {
        if ( N == 0 )
          return;
        const double dN( N );
        const size_t idx = this->valueToBin( val );
        vectAt(m_content,idx) += dN;
        if ( opt_allow_weights == AllowWeights::YES)
          vectAt(m_errors,idx) += dN;
        if ( opt_allow_weights == AllowWeights::YES_BUT_DELAY_ALLOC
             && !m_errors.empty() )
          vectAt(m_errors,idx) += dN;
      }

      template<AllowWeights U = opt_allow_weights>
      void fill( double val, double weight,
                 typename std::enable_if<U!=AllowWeights::NO>::type* = nullptr )
      {
        if ( !(weight>0.0) )
          return;
        if ( opt_allow_weights == AllowWeights::YES_BUT_DELAY_ALLOC ) {
          if ( weight == 1.0 )
            this->fill(val);
          if ( m_errors.empty() )
            m_errors = m_content;
        }
        const size_t idx = this->valueToBin( val );
        vectAt(m_content,idx) += weight;
        vectAt(m_errors,idx) += weight*weight;
      }

      void dump_metadata( std::ostream& os ) const
      {
        os << "HistBinData1D(nbins="<<m_nbins
           <<",xmin="<<fmtg(m_xmin)
           <<",xmax="<<fmtg(m_xmax)<<"):\n";
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        if ( nbins_overflow==2 ) {
          const double underflow = vectAt(m_content,0);
          const double overflow =  vectAt(m_content,m_nbins+1);
          os << "  overflow handling: RECORD\n"
            "  underflow : "<<fmtg(underflow)
             << "\n  overflow  : "<<fmtg(overflow)<<'\n';
        } else {
          os << "  overflow handling: CLAMP\n";
        }
      }

      void dump_contents( std::ostream& os) const
      {
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        for (size_t ibin = 0; ibin < m_nbins; ++ibin)
          os << "  content[ibin="<<ibin<<"] : "
             <<fmtg(this->getBinContent(ibin))
             <<" +- "<<fmtg(this->getBinError(ibin))<<'\n';
      }

      void dump( std::ostream& os ) const
      {
        os << this->dump_metadata()
           << this->dump_contents();
      }

    private:

      size_t valueToBin( double val ) {
        nc_assert( !std::isnan(val) );
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        if ( nbins_overflow == 0 ) {
          auto idx = static_cast<size_t>( m_invDelta
                                          * ncclamp(val-m_xmin,m_clampRange) );
          nc_assert( idx < m_nbins );
          return idx;
        } else  {
          if ( val < m_xmin ) return 0;
          if ( val >= m_xmax ) return ( val == m_xmax ? m_nbins : m_nbins+1 );
          return 1 + std::min<size_t>(m_nbins,
                                      static_cast<size_t>( m_invDelta
                                                           * ( val-m_xmin ) ));
        }
      }

      opt_storate_type m_content;//sum weight for each bin
      //sum weight-squares for each bin, if weights are allowed:
      typename std::conditional<opt_allow_weights==AllowWeights::NO,
                                detail::EmptyVect,
                                opt_storate_type>::type m_errors;
      double m_xmin;
      double m_xmax;
      double m_invDelta;
      std::size_t m_nbins;
      typename std::conditional<opt_allow_weights==AllowWeights::NO,
                                detail::EmptyPairDD,
                                PairDD>::type m_clampRange;
    };

    ///////////////////////////////////////////////////////////////////////////
    // RunningStats1D class                                                  //
    ///////////////////////////////////////////////////////////////////////////

    class RunningStats1D final {
    public:

      // Calculate one-pass mean, sum and RMS of values.
      //
      // This could in principle be achieved by using the same simple approach
      // to one-pass calculation of RMS used for example in ROOT (using sumW,
      // sumWX and sumWX2 variables). However, this is highly numerically
      // unstable, since one in the end calculates RMS by
      // sqrt(bignum^2-otherbignum^2).
      //
      // In the implementation here we instead define a different variable:
      //
      // T = sumWX2 - (sumWX)^2/sumW  => RMS = sqrt(T/sumW)
      //
      // One can (carefully) derive formulas for how T should be updated when
      // respectively registering new values and merging two T values. These
      // formulas are numerically stable, and are represented in the
      // implementation below.
      //
      // Below we keep T in the variable called "m_rms_state", but for reference
      // we also provide the simple ROOT-style approach by defining
      // NCRYSTAL_HIST_ROOT_STYLE_RMS (***FOR TESTING ONLY***). In that case,
      // rms_state corresponds to the usual sumWX2, and one can expect
      // significant levels of incorrectness in the RMS estimations when
      // distributions have true rms << |mean|.
      //
      // Thomas Kittelmann, 2014, adapted for NCrystal 2024.

      //Record on-the-fly values, or merge data from another instance:
      void registerValue( double val );
      void registerValue( double val, double weight );
      void registerNValues( double val, std::size_t N );
      void merge( const RunningStats1D& o );

      //Access data (note that Mean, RMS, maxFilled, and minFilled can only be
      //calculated if the histogram has data):
      bool hasData() const;
      double getIntegral() const;
      double calcMean() const;
      double calcRMS() const;
      double calcRMSSq() const;//faster by a sqrt than RMS
      double maxFilled() const;
      double minFilled() const;

      void toJSON(std::ostringstream&) const;//all stats as JSON encoded dict

    private:
      double m_sumw = 0.0;
      double m_sumwx = 0.0;
      double m_rms_state = 0.0;
      double m_maxfilled = -1.0;//use max<min to indicate no fills yet
      double m_minfilled = 1.0;
      void update_maxmin( double val );
    };

    ///////////////////////////////////////////////////////////////////////////
    // Hist1D class                                                          //
    ///////////////////////////////////////////////////////////////////////////

    template<
      AllowWeights ALLOW_WEIGHTS = AllowWeights::YES,
      OverflowHandling OF_HANDLING = OverflowHandling::Record,
      class TStorage = std::vector<double>
      >
    class Hist1D final {
    public:
      constexpr static auto opt_allow_weights = ALLOW_WEIGHTS;
      constexpr static auto opt_of_handling = OF_HANDLING;
      using opt_storate_type = TStorage;
      using bindata1d_t = HistBinData1D<opt_allow_weights,
                                        opt_of_handling,
                                        opt_storate_type>;
    private:
      bindata1d_t m_bindata;
      RunningStats1D m_stats;
      std::string m_title;
    public:
      Hist1D( size_t nbins, double xmin, double xmax, std::string title = {} )
        : m_bindata(nbins,xmin,xmax),
          m_title(std::move(title))
      {
      }

      const std::string& title() const { return m_title; }
      const bindata1d_t& binData() const { return m_bindata; }
      const RunningStats1D& stats() const { return m_stats; }

      void fill(double val) {
        m_bindata.fill(val);
        m_stats.registerValue( val );
      }

      void fillN(double val, size_t N )
      {
        m_bindata.fillN(val,N);
        m_stats.registerNValues( val, N );
      }

      template<AllowWeights U = opt_allow_weights>
      void fill( double val, double weight,
                 typename std::enable_if<U!=AllowWeights::NO>::type* = nullptr )
      {
        m_bindata.fill(val,weight);
        m_stats.registerValue( val, weight );
      }

      void merge( const Hist1D& o )
      {
        m_bindata.merge(o.m_bindata);
        if ( m_title != o.m_title )
          NCRYSTAL_THROW(CalcError,"Attempting to merge incompatible 1D"
                         " histograms (titles are different)");
        m_stats.merge(o.m_stats);
      }

      void toJSON( std::ostringstream& os,
                   bool include_data_arrays = true ) const
      {
        streamJSONDictEntry( os, "title", m_title, JSONDictPos::FIRST);
        os<< ',';
        streamJSON(os,"stats");
        os << ':';
        m_stats.toJSON(os);
        os << ',';
        streamJSON(os,"bindata");
        os << ':';
        m_bindata.toJSON(os,include_data_arrays);
        os << '}';
      }
    };

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace Hists {

    inline bool RunningStats1D::hasData() const {
      nc_assert( bool(m_sumw) == bool(m_maxfilled >= m_minfilled) );
      return bool(m_sumw);
    }

    inline double RunningStats1D::getIntegral() const
    {
      return m_sumw;
    }

    inline double RunningStats1D::calcRMS() const
    {
#ifdef NCRYSTAL_HIST_ROOT_STYLE_RMS
      //Unstable calc needs ncabs protection:
      return std::sqrt( ncabs( calcRMSSq() ) );
#else
      double rms2 = calcRMSSq();
      nc_assert( rms2 >= 0.0 );
      return std::sqrt( rms2 );
#endif
    }

    inline double RunningStats1D::maxFilled() const
    {
      if ( m_maxfilled < m_minfilled )
        NCRYSTAL_THROW(CalcError,"maxFilled not well defined in empty histograms");
      return m_maxfilled;
    }

    inline double RunningStats1D::minFilled() const
    {
      if ( m_maxfilled < m_minfilled )
        NCRYSTAL_THROW(CalcError,"minFilled not well defined in empty histograms");
      return m_minfilled;
    }

    inline double RunningStats1D::calcMean() const
    {
      if (!hasData())
        NCRYSTAL_THROW(CalcError,"Mean not well defined in empty histograms");
      return m_sumwx / m_sumw;
    }

    inline void RunningStats1D::registerValue( double val )
    {
      update_maxmin(val);
#ifdef NCRYSTAL_HIST_ROOT_STYLE_RMS
      m_rms_state += val*val;
#else
      const double d1 = m_sumw * val - m_sumwx;
      const double d2 = m_sumw * ( 1.0 + m_sumw );
      if ( d2 )
        m_rms_state += d1 * d1 / d2;
#endif
      ++m_sumw;
      m_sumwx += val;
    }

    inline void RunningStats1D::registerValue( double val, double weight )
    {
      if ( ! (weight>0.0) )
        return;
      update_maxmin(val);
#if NCRYSTAL_HIST_ROOT_STYLE_RMS
      m_rms_state += weight*val*val;
#else
      const double d1 = m_sumw * val - m_sumwx;
      const double d2 = m_sumw * ( weight + m_sumw );
      if ( d2 )
        m_rms_state += weight * ( d1 * d1 / d2 );
#endif
      m_sumw += weight;
      m_sumwx += val * weight;
    }

    inline void RunningStats1D::update_maxmin( double val )
    {
      if ( m_maxfilled < m_minfilled ) {
        //first fill:
        m_minfilled = m_maxfilled = val;
      } else {
        //subsequent fill:
        m_minfilled = std::min<double>(val,m_minfilled);
        m_maxfilled = std::max<double>(val,m_maxfilled);
      }
    }

  }
}

#endif
