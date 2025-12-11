#ifndef NCrystal_Hists_hh
#define NCrystal_Hists_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCMath.hh"
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
        using value_type = double;
        using size_type = std::size_t;
        void resize(size_type, double) {};
        constexpr bool empty() const { return true; }
        constexpr size_type size() const { return 0; }
        const double* data() const { return nullptr; }
        const double& at(size_type) const
        {
          nc_assert_always(false);
          static double tmp; return tmp;
        }
        double& at(size_type)
        {
          nc_assert_always(false);
          static double tmp; return tmp;
        }
        const double& operator[](size_type i) const { return this->at(i); }
        double& operator[](size_type i) { return this->at(i); }
        ncconstexpr17 EmptyVect& operator=( const EmptyVect&) { return *this; }
      };
      struct EmptyPairDD { void swap( PairDD& ) {} };
      inline double pair_first( const PairDD& p ) { return p.first; }
      inline double pair_second( const PairDD& p ) { return p.second; }
      inline double pair_first( const EmptyPairDD& ) { return 0.0; }
      inline double pair_second( const EmptyPairDD& ) { return 0.0; }
    }

    ///////////////////////////////////////////////////////////////////////////
    // HistBinData1D class                                                   //
    ///////////////////////////////////////////////////////////////////////////

    struct Binning {
      using size_type = std::uint_least32_t;
      Binning( size_type nbins_, double xmin_, double xmax_ ) ncnoexceptndebug
        : nbins(nbins_), xmin(xmin_), xmax(xmax_)
      {
#ifndef NDEBUG
        validate();
#endif
      }
      bool operator==( const Binning& o ) const ncnoexceptndebug
      {
#ifndef NDEBUG
        validate();
        o.validate();
#endif
        nc_assert( std::isfinite( xmin ) && std::isfinite( xmax ) );
        return nbins == o.nbins && xmin == o.xmin && xmax == o.xmax;
      }
      bool operator<( const Binning& o ) const ncnoexceptndebug
      {
#ifndef NDEBUG
        validate();
        o.validate();
#endif
        return ( nbins != o.nbins ? ( nbins < o.nbins )
                 : ( xmin != o.xmin ? xmin < o.xmin : xmax < o.xmax ) );
      }
      void validate() const;
      size_type nbins;
      double xmin;
      double xmax;
    };
    std::ostream& operator<<( std::ostream&, const Binning& );

    template<
      AllowWeights ALLOW_WEIGHTS = AllowWeights::YES,
      OverflowHandling OF_HANDLING = OverflowHandling::Record,
      class TStorage = std::vector<double>
      >
    class HistBinData1D final {
      static_assert( std::is_same<typename TStorage::value_type,double>::value, "");
    public:
      using size_type = Binning::size_type;
      constexpr static AllowWeights opt_allow_weights = ALLOW_WEIGHTS;
      constexpr static OverflowHandling opt_of_handling = OF_HANDLING;
      using opt_storage_type = TStorage;
    private:
      static constexpr size_type
      nbins_overflow = ( opt_of_handling == OverflowHandling::Record ? 2 : 0 );
    public:

      const Binning& getBinning() const
      {
        return m_b;
      }

      HistBinData1D( size_type nbins, double xmin, double xmax )
        : HistBinData1D( Binning( nbins, xmin, xmax ) )
      {
      }

      HistBinData1D( const Binning& b )
        : m_b(b)
      {
        b.validate();
        m_content.resize(m_b.nbins + nbins_overflow,0.0);
        if ( opt_allow_weights == AllowWeights::YES )
          m_errors.resize(m_b.nbins + nbins_overflow,0.0);
        const double delta = (m_b.xmax-m_b.xmin) / m_b.nbins;
        nc_assert_always(delta>0.0);
        m_invDelta = 1.0 / delta;
        if ( nbins_overflow == 0 ) {
          PairDD cr( 0.99*delta, m_b.xmax - (0.99*delta+m_b.xmin) );
          m_clampRange.swap(cr);
        }
      }

      void toJSON( std::ostream& os, bool include_data_arrays ) const
      {
        streamJSONDictEntry(os, "xmin", m_b.xmin,
                            JSONDictPos::FIRST);
        streamJSONDictEntry(os, "xmax", m_b.xmax);
        streamJSONDictEntry(os, "nbins", m_b.nbins);

        if ( include_data_arrays ) {
          streamJSONDictEntry(os, "content", this->getContents());
          if ( m_errors.empty() ) {
            streamJSONDictEntry(os, "errorsq", json_null_t{});
          } else {
            streamJSONDictEntry(os, "errorsq", this->getErrorsSquared());
          }
        }

        if ( opt_of_handling == OverflowHandling::Record ) {
          nc_assert_always( nbins_overflow == 2 );
          streamJSONDictEntry(os, "underflow", vectAt(m_content,0) );
          streamJSONDictEntry(os, "overflow", vectAt(m_content,m_b.nbins+1));
          if ( opt_allow_weights == AllowWeights::NO || m_errors.empty() ) {
            streamJSONDictEntry(os, "underflow_errorsq", vectAt(m_content,0) );
            streamJSONDictEntry(os, "overflow_errorsq", vectAt(m_content,m_b.nbins+1));
          } else {
            nc_assert_always(m_errors.size()==m_b.nbins+2);
            streamJSONDictEntry(os, "underflow_errorsq", vectAt(m_errors,0) );
            streamJSONDictEntry(os, "overflow_errorsq", vectAt(m_errors,m_b.nbins+1));
          }
        }
        os << '}';
      }

      bool mergeCompatible( const HistBinData1D& o ) const
      {
        return m_b == o.m_b;
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

      double getBinContent( size_type ibin) const
      {
        return vectAt(m_content,(nbins_overflow==0?ibin:ibin+1));
      }

      double getBinError( size_type ibin) const
      {
        constexpr size_type offset = ( nbins_overflow==0 ? 0 : 1 );
        if ( opt_allow_weights == AllowWeights::YES )
          return std::sqrt( vectAt(m_errors,ibin+offset) );
        if ( opt_allow_weights == AllowWeights::NO )
          return std::sqrt( vectAt(m_content,ibin+offset) );
        return ( m_errors.empty()
                 ? std::sqrt( vectAt(m_content,ibin+offset) )
                 : std::sqrt( vectAt(m_errors,ibin+offset) ) );
      }

      double getXMin() const { return m_b.xmin; }
      double getXMax() const { return m_b.xmax; }
      unsigned getNBins() const { return m_b.nbins; }

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
        return vectAt(m_content,m_b.nbins+1);
      }

      template<OverflowHandling U = opt_of_handling>
      double getOverflowErrorSquared( typename std::enable_if<
                                      U==OverflowHandling::Record>::type* = nullptr ) const
      {
        if ( opt_allow_weights == AllowWeights::NO || m_errors.empty() )
          return vectAt(m_content,m_b.nbins+1);
        else
          return vectAt(m_errors,m_b.nbins+1);
      }

      template<OverflowHandling U = opt_of_handling>
      double getUnderflowErrorSquared( typename std::enable_if<
                                       U==OverflowHandling::Record>::type* = nullptr ) const
      {
        if ( opt_allow_weights == AllowWeights::NO || m_errors.empty() )
          return vectAt(m_content,0);
        else
          return vectAt(m_errors,0);
      }

      Span<const double> getContents() const
      {
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        return ( nbins_overflow==0
                 ? Span<const double>{ m_content.data(),
                                       m_content.data() + m_b.nbins }
                 : Span<const double>{ m_content.data() + 1,
                     m_content.data() + m_b.nbins + 1 } );
      }

      Span<const double> getErrorsSquared() const
      {
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        if ( opt_allow_weights == AllowWeights::NO || m_errors.empty() )
          return getContents();
        return ( nbins_overflow==0
                 ? Span<const double>{ m_errors.data(),
                                       m_errors.data() + m_b.nbins }
                 : Span<const double>{ m_errors.data() + 1,
                     m_errors.data() + m_b.nbins + 1 } );
      }

      void fill(double val) {
        size_type idx = this->valueToBin( val );
        vectAt(m_content,idx) += 1.0;
        if ( opt_allow_weights == AllowWeights::YES)
          vectAt(m_errors,idx) += 1.0;
        if ( opt_allow_weights == AllowWeights::YES_BUT_DELAY_ALLOC
             && !m_errors.empty() )
          vectAt(m_errors,idx) += 1.0;
      }

      void fillN(double val, double N )
      {
        if ( !(N>0) )
          return;
        nc_assert(std::isfinite(N));
        const size_type idx = this->valueToBin( val );
        vectAt(m_content,idx) += N;
        if ( opt_allow_weights == AllowWeights::YES)
          vectAt(m_errors,idx) += N;
        if ( opt_allow_weights == AllowWeights::YES_BUT_DELAY_ALLOC
             && !m_errors.empty() )
          vectAt(m_errors,idx) += N;
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
        const size_type idx = this->valueToBin( val );
        vectAt(m_content,idx) += weight;
        vectAt(m_errors,idx) += weight*weight;
      }

      void dump_metadata( std::ostream& os ) const
      {
        os << "HistBinData1D(nbins="<<m_b.nbins
           <<",xmin="<<fmtg(m_b.xmin)
           <<",xmax="<<fmtg(m_b.xmax)<<"):\n";
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        if ( nbins_overflow==2 ) {
          const double underflow = vectAt(m_content,0);
          const double overflow =  vectAt(m_content,m_b.nbins+1);
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
        for (size_type ibin = 0; ibin < m_b.nbins; ++ibin)
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

      size_type valueToBin( double val ) {
        nc_assert( !std::isnan(val) );
        static_assert( nbins_overflow==0 || nbins_overflow==2,"");
        if ( nbins_overflow == 0 ) {
          const double clampval = ncclamp(val-m_b.xmin,
                                          detail::pair_first(m_clampRange),
                                          detail::pair_second(m_clampRange));
          auto idx = static_cast<size_type>(m_invDelta * clampval);
          nc_assert( idx < m_b.nbins );
          return idx;
        } else  {
          if ( val < m_b.xmin )
            return 0;
          if ( val >= m_b.xmax )
            return ( val == m_b.xmax ? m_b.nbins : m_b.nbins+1 );
          return 1 + std::min<size_type>(m_b.nbins,
                                      static_cast<size_type>( m_invDelta
                                                           * ( val-m_b.xmin ) ));
        }
      }

      opt_storage_type m_content;//sum weight for each bin
      //sum weight-squares for each bin, if weights are allowed:
      typename std::conditional<opt_allow_weights==AllowWeights::NO,
                                detail::EmptyVect,
                                opt_storage_type>::type m_errors;
      Binning m_b;
      double m_invDelta;
      typename std::conditional<opt_of_handling==OverflowHandling::Record,
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
      void registerNValues( double val, double N );
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

      void toJSON(std::ostream&) const;//all stats as JSON encoded dict

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
      static_assert( std::is_same<typename TStorage::value_type,double>::value, "");
    public:
      constexpr static auto opt_allow_weights = ALLOW_WEIGHTS;
      constexpr static auto opt_of_handling = OF_HANDLING;
      using opt_storage_type = TStorage;
      using bindata1d_t = HistBinData1D<opt_allow_weights,
                                        opt_of_handling,
                                        opt_storage_type>;
      using size_type = typename bindata1d_t::size_type;

    private:
      bindata1d_t m_bindata;
      RunningStats1D m_stats;
      std::string m_title;
    public:

      Hist1D( const Binning& b, std::string title = {} )
        : m_bindata(b),
          m_title(std::move(title))
      {}

      Hist1D( size_type nbins, double xmin, double xmax, std::string title = {} )
        : m_bindata(nbins,xmin,xmax),
          m_title(std::move(title))
      {
      }

      const std::string& title() const { return m_title; }
      const bindata1d_t& binData() const { return m_bindata; }
      const RunningStats1D& stats() const { return m_stats; }
      const Binning& binning() const { return m_bindata.getBinning(); }

      void fill(double val) {
        m_bindata.fill(val);
        m_stats.registerValue( val );
      }

      void fillN(double val, double N )
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

      void toJSON( std::ostream& os,
                   bool include_data_arrays = true ) const
      {
        streamJSONDictEntry( os, "datatype", "NCrystalHist1D_v1",
                             JSONDictPos::FIRST);
        streamJSONDictEntry( os, "title", m_title );
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
#ifdef NCRYSTAL_HIST_ROOT_STYLE_RMS
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
