#ifndef NCrystal_SABEval_hh
#define NCrystal_SABEval_hh

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

#include "NCrystal/internal/sab/NCSABUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    //Provide utilities for evaluation S(alpha,beta) tables, either within a
    //single grid cell (SABCellEval), or an entire table (SABEval).

    enum class SABInterpolationOrder { BETA_FIRST, ALPHA_FIRST };
    enum class KBStatus { FullyInside, FullyOutside, Crossing };
    enum class IgnorePart { Above, Below };

    template< InterpolationScheme alphaInterpType_ = InterpolationScheme::LOGLIN,
              SABInterpolationOrder interpOrder_ = SABInterpolationOrder::ALPHA_FIRST>
    class SABCellEval {
    public:
      struct SCE_Data;

      static constexpr auto alphaInterpType = SCE_Data::alphaInterpType;
      static constexpr auto betaInterpType = SCE_Data::betaInterpType;
      static constexpr auto interpOrder = SCE_Data::interpOrder;
      typedef double SValues[4];

      SABCellEval( PairDD alpha, PairDD beta, const double (&svals)[4] );
      //svals: { S(alpha0,beta0), S(alpha1,beta0), S(alpha0,beta1), S(alpha1,beta1) };

      //Base parameters:
      const PairDD& alpha() const noexcept { return m_data.alpha; }
      const PairDD& beta() const noexcept { return m_data.beta; }
      typedef double DblArr4[4];
      const DblArr4& sValues() const noexcept { return m_data.S; }
      double sMax() const;

      //Evaluate at a point (the point must belong to the cell!):
      double eval( PairDD alpha_beta ) const;
      double eval( double alpha, double beta ) const;

      //Sample point within cell according to height of eval(alpha,beta).
      PairDD sample( RNG& ) const;

      //Get integral of entire cell (directly or appended to a StableSum
      //instance):
      double integral() const;
      void addIntegral( StableSum& ) const;

      //Check if cell is fully inside or outside the kinematic reach (kinematic
      //boundary, KB) of a given neutron:
      KBStatus kbStatus( double ekin_div_kT ) const;

      //Only consider parts within kinematic reach of a given neutron:
      double integralWithinKinematicBounds( double ekin_div_kT ) const;
      void addIntegralWithinKinematicBounds( StableSum&, double ekin_div_kT ) const;
      double sOverlayValueWithinKinematicBounds( double ekin_div_kT ) const;//returns 0.0 if E/kT=0

      //Or further restrict the integral to the region below a given betamax
      //value (e.g. to find the contribution to ultra cold neutron production):
      double integralWithinKinematicBoundsBelowBetamax( double ekin_div_kT, double betamax ) const;
      void addIntegralWithinKinematicBoundsBelowBetamax( StableSum&, double ekin_div_kT, double betamax ) const;
      double sOverlayValueWithinKinematicBoundsBelowBetamax( double ekin_div_kT, double betamax ) const;//returns 0.0 if no area

      //Split cell:
      SABCellEval splitAtBeta( double betaval, IgnorePart ip ) const { return { m_data.splitAtBeta(betaval,ip) }; }
      SABCellEval splitAtAlpha( double alphaval, IgnorePart ip ) const { return { m_data.splitAtAlpha(alphaval,ip) }; }

    private:
      SCE_Data m_data;
      SABCellEval( SCE_Data&& dd ) : m_data(std::move(dd)) {}
    };


    //strongly type alphaGrid.size()-1:
    class NAlphaCells {
    public:
      NAlphaCells( const VectD& alphaGrid ) ncnoexceptndebug;
      explicit constexpr NAlphaCells( unsigned ) noexcept;
      constexpr unsigned value() const noexcept { return m_value; }
    private:
      unsigned m_value;
    };

    //Evaluate at a point (returns 0.0 outside the grid!):
    using rawcellidx_t = unsigned;
    class CellIndex;
    class PackedCellIndex;
    class NAlphaCells;

    CellIndex getCellIndex(const SABData&, double alpha, double beta );

    class PackedCellIndex {
    public:
      //Cell indices are in unpacked form (ialpha,ibeta) where the indices are
      //the ones of the lower grid point. Thus, if the grid has nalpha*nbeta
      //grid points, there are only (nalpha-1)*(nbeta-1) cells. We thus define
      //nalphacells=nalpha-1 and ncellbeta=nbeta-1, and let ialpha run from
      //0..(nalphacells-1), and similarly for ibeta.

      using index_t = rawcellidx_t;
      static constexpr index_t invalid = std::numeric_limits<index_t>::max();
      constexpr PackedCellIndex() noexcept : m_idx(invalid) {}
      constexpr PackedCellIndex( NullOptType ) noexcept : m_idx(invalid) {}
      constexpr PackedCellIndex( NAlphaCells, CellIndex ) noexcept;
      CellIndex unpack( NAlphaCells ) const ncnoexceptndebug;
      constexpr bool isValid() const noexcept { return m_idx != invalid; }
      constexpr index_t value() const noexcept { return m_idx; }
      static index_t nRawMax( NAlphaCells, std::size_t nbeta );
      struct from_raw_t {};
      constexpr PackedCellIndex( from_raw_t, index_t raw ) noexcept : m_idx(raw) { }
      index_t rawValue() const noexcept { return m_idx; }
    private:
      index_t m_idx;
    };

    class CellIndex {
    public:
      using index_t = rawcellidx_t;
    private:
      static constexpr index_t invalid = std::numeric_limits<index_t>::max();
    public:
      constexpr CellIndex( NullOptType ) noexcept : m_ia{ invalid }, m_ib{ invalid } {}
      CellIndex( index_t ia, index_t ib )  ncnoexceptndebug;
      PackedCellIndex pack( NAlphaCells nac ) const noexcept;
      index_t packRaw( NAlphaCells nac ) const noexcept;
      constexpr bool isValid() const noexcept { return m_ia != invalid; }
      constexpr index_t ia() const noexcept { return m_ia; }
      constexpr index_t ib() const noexcept { return m_ib; }
    private:
      index_t m_ia, m_ib;
    };


    template< InterpolationScheme alphaInterpType_ = InterpolationScheme::LOGLIN,
              SABInterpolationOrder interpOrder_ = SABInterpolationOrder::ALPHA_FIRST>
    class SABEval {
    public:
      using celleval_t = SABCellEval<alphaInterpType_,interpOrder_>;

      SABEval( shared_obj<const SABData> );
      SABEval( const SABData* );//Warning: lifetime of SABData must outlive SABEval object!!
      const SABData& sabData() const;
      rawcellidx_t nRawCellIdxMax() const;
      celleval_t getCell( rawcellidx_t cell_rawidx ) const;
      celleval_t getCell( PackedCellIndex ) const;
      celleval_t getCell( CellIndex ) const;
      Optional<celleval_t> getCellEval(double alpha, double beta) const
      {
        //Not inlined below due to VisualStudio issues.
        Optional<celleval_t> cell;
        auto idx = getCellIndex(*m_sab,alpha,beta);
        if ( idx.isValid() ) nclikely {
          //inside grid
          cell.emplace( getCell( idx ) );
        }
        return cell;
      }

      double eval( PairDD alpha_beta ) const;
      double eval( double alpha, double beta ) const;

    private:
      const SABData* m_sab;
      std::shared_ptr<const SABData> m_sabptr;//to keep alive
      NAlphaCells m_nalphacells;
    };

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    namespace detail_sce {
      struct SCE_SLogVals {
        static constexpr bool hasLogValues = true;
        double logS[4];
        void update_logS( int i, double sval ) ncnoexceptndebug
        {
          nc_assert( sval >= 0.0 && !ncisnanorinf(sval) );
          nc_assert(i<4);
          this->logS[i] = ( sval > 0.0 ? std::log( sval ) : -kInfinity );
        }
      };
      struct SCE_Empty {
        static constexpr bool hasLogValues = false;
        void update_logS( int, double ) noexcept {}
      };
    }

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    struct SABCellEval<AIT,IO>::SCE_Data
      : public std::conditional<( IO==SABInterpolationOrder::ALPHA_FIRST && AIT==InterpolationScheme::LOGLIN ),
                                detail_sce::SCE_SLogVals, detail_sce::SCE_Empty >::type {
      static constexpr auto alphaInterpType = AIT;
      static constexpr auto betaInterpType = InterpolationScheme::LINLIN;//always linear (for now?)
      static constexpr auto interpOrder = IO;
      PairDD alpha, beta;//todo: encapsulate data members?
      double S[4];
      void validate() const;
      SCE_Data(PairDD alpha_, PairDD beta_, const double (&svals)[4])
        : alpha(alpha_), beta(beta_)
      {
        for ( auto i : ncrange(4) ) {
          S[i] = svals[i];
          this->update_logS( i, svals[i] );
        }
#ifndef NDEBUG
        validate();
#endif
      }
      double eval( double alpha, double beta ) const;
      double integral() const;
      void integral(StableSum&) const;
      SCE_Data splitAtBeta( double beta, IgnorePart ) const;
      SCE_Data splitAtAlpha( double alpha, IgnorePart ) const;
      PairDD sample( RNG& ) const;
      void integralWKB( StableSum&, double ekin_div_kT ) const;
      double sOverlayWKB( double ekin_div_kT ) const;
      KBStatus kbStatus( double ekin_div_kT ) const;
      double sMax() const { return ncmax(S[0],S[1],S[2],S[3]); }
    };

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline SABCellEval<AI,IO>::SABCellEval( PairDD alpha, PairDD beta, const double (&svals)[4] ) : m_data(alpha, beta, svals) {}
    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline double SABCellEval<AI,IO>::eval( PairDD alpha_beta ) const { return m_data.eval(alpha_beta.first,alpha_beta.second); }
    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline double SABCellEval<AI,IO>::eval( double alpha, double beta ) const { return m_data.eval(alpha,beta); }
    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline double SABCellEval<AI,IO>::integral() const { return m_data.integral(); }
    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline double SABCellEval<AI,IO>::sMax() const { return m_data.sMax(); }
    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline void SABCellEval<AI,IO>::addIntegral( StableSum& sum ) const { m_data.integral(sum); }
    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline PairDD SABCellEval<AI,IO>::sample( RNG& rng ) const { return m_data.sample(rng); }

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    KBStatus SABCellEval<AI,IO>::kbStatus( double ekin_div_kT ) const
    {
      return m_data.kbStatus(ekin_div_kT);
    }

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline double SABCellEval<AI,IO>::integralWithinKinematicBoundsBelowBetamax( double ekin_div_kT, double betamax ) const
    {
      if ( m_data.beta.first >= betamax || m_data.beta.second <= -ekin_div_kT )
        return 0.0;//likely excludes most cells if [ekin_div_kT,betamax] is a narrow threshold.
      StableSum sum;
      if ( m_data.beta.second <= betamax ) {
        m_data.integralWKB( sum, ekin_div_kT );
      } else {
        m_data.splitAtBeta( betamax, IgnorePart::Above ).integralWKB( sum, ekin_div_kT );
      }
      return sum.sum();
    }

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    double SABCellEval<AI,IO>::sOverlayValueWithinKinematicBounds( double ekin_div_kT ) const
    {
      if ( !(ekin_div_kT>0.0) || m_data.beta.second <= -ekin_div_kT )
        return 0.0;
      return m_data.sOverlayWKB( ekin_div_kT );
    }

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    double SABCellEval<AI,IO>::sOverlayValueWithinKinematicBoundsBelowBetamax( double ekin_div_kT, double betamax ) const
    {
      if ( m_data.beta.first >= betamax || m_data.beta.second <= -ekin_div_kT )
        return 0.0;
      if ( m_data.beta.second <= betamax ) {
        return m_data.sOverlayWKB( ekin_div_kT );
      } else {
        return m_data.splitAtBeta( betamax, IgnorePart::Above ).sOverlayWKB( ekin_div_kT );
      }
    }

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline void SABCellEval<AI,IO>::addIntegralWithinKinematicBoundsBelowBetamax( StableSum& sum, double ekin_div_kT, double betamax ) const
    {
      if ( m_data.beta.first >= betamax || m_data.beta.second <= -ekin_div_kT )
        return;//likely excludes most cells if [ekin_div_kT,betamax] is a narrow threshold.
      if ( m_data.beta.second <= betamax ) {
        m_data.integralWKB( sum, ekin_div_kT );
      } else {
        m_data.splitAtBeta( betamax, IgnorePart::Above ).integralWKB( sum, ekin_div_kT );
      }
    }

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    inline double SABCellEval<AI,IO>::integralWithinKinematicBounds( double ekin_div_kT ) const
    {
      if ( m_data.beta.second <= -ekin_div_kT )
        return 0.0;
      StableSum sum;
      m_data.integralWKB( sum, ekin_div_kT );
      return sum.sum();
    }

    template< InterpolationScheme AI, SABInterpolationOrder IO>
    void SABCellEval<AI,IO>::addIntegralWithinKinematicBounds( StableSum& sum, double ekin_div_kT ) const
    {
      if ( !(m_data.beta.second <= -ekin_div_kT) )
        m_data.integralWKB( sum, ekin_div_kT );
    }

    inline NAlphaCells::NAlphaCells( const VectD& alphaGrid ) ncnoexceptndebug
    : m_value( static_cast<unsigned>(alphaGrid.size()-1) )
    {
      nc_assert( alphaGrid.size()>1);
      nc_assert( alphaGrid.size()<std::numeric_limits<unsigned>::max() );
    }

    inline constexpr NAlphaCells::NAlphaCells( unsigned value ) noexcept
      : m_value( value )
    {
    }

    inline CellIndex::CellIndex( index_t ia, index_t ib ) ncnoexceptndebug
      : m_ia{ ia }, m_ib{ ib }
    {
      nc_assert( bool(ia==invalid) == bool(ib==invalid) );
    }

    inline PackedCellIndex CellIndex::pack( NAlphaCells nac ) const noexcept
    {
      return PackedCellIndex{ nac, *this };
    }

    inline PackedCellIndex::index_t CellIndex::packRaw( NAlphaCells nac ) const noexcept
    {
      return PackedCellIndex{ nac, *this }.rawValue();
    }

    inline PackedCellIndex::index_t PackedCellIndex::nRawMax( NAlphaCells nac, std::size_t nbeta )
    {
      nc_assert(nac.value()>0 && nbeta>1 );
      return nac.value()*(nbeta-1);
    }

    inline constexpr PackedCellIndex::PackedCellIndex( NAlphaCells nac, CellIndex ci ) noexcept
      : m_idx{ ci.isValid() ? static_cast<index_t>( ci.ia() + nac.value() * ci.ib() ) : invalid }
    {
    }

    inline CellIndex PackedCellIndex::unpack( NAlphaCells nac ) const ncnoexceptndebug
    {
      return CellIndex{ static_cast<index_t>(m_idx%nac.value()),
                        static_cast<index_t>(m_idx/nac.value()) };
    }

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline SABEval<AIT,IO>::SABEval( shared_obj<const SABData> sab )
      : m_sab(sab.get()),
        m_sabptr(std::move(sab)),
        m_nalphacells(m_sab->alphaGrid())
    {}

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline SABEval<AIT,IO>::SABEval( const SABData* sab )
      : m_sab(sab),
        m_sabptr(nullptr),
        m_nalphacells(m_sab->alphaGrid())
    {}

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline const SABData& SABEval<AIT,IO>::sabData() const { return *m_sab; }

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline PackedCellIndex::index_t SABEval<AIT,IO>::nRawCellIdxMax() const
    {
      return PackedCellIndex::nRawMax( m_nalphacells,
                                       m_sab->betaGrid().size() );
    }

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline typename SABEval<AIT,IO>::celleval_t SABEval<AIT,IO>::getCell( PackedCellIndex::index_t cell_rawidx ) const
    {
      return getCell( PackedCellIndex( PackedCellIndex::from_raw_t{}, cell_rawidx ) );
    }

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline typename SABEval<AIT,IO>::celleval_t SABEval<AIT,IO>::getCell( PackedCellIndex idx ) const
    {
      return getCell( idx.unpack( m_nalphacells ) );
    }

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline double SABEval<AIT,IO>::eval( PairDD alpha_beta ) const { return eval(alpha_beta.first,alpha_beta.second); }

    template< InterpolationScheme AIT, SABInterpolationOrder IO>
    inline double SABEval<AIT,IO>::eval( double alpha, double beta ) const
    {
      auto ce = getCellEval(alpha,beta);
      return ( ce.has_value() ? ce.value().eval(alpha,beta) : 0.0 );//0.0 if outside grid
    }


  }
}

#endif
