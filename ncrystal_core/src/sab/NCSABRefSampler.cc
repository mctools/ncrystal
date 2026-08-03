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

#include "NCrystal/internal/sab/NCSABRefSampler.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace {

    //TODO: These randShuffle utilities could be moved to somewhere common, if
    //needed.
    template<class TSwapIJ>
    void randShuffleFisherYates(RNG& rng, std::size_t n, TSwapIJ&& swapIJ) {
      //A simple Fisher-Yates (Knuth) shuffle, but written to accomodate our
      //particular RNG interface which typically provides uint32_t values more
      //cheaply than uint64_t values, and with a flexible swap function which
      //for instance allows multiple containers to be shuffled at once in
      //randShuffleMultiple.
      if (n < 2)
        return;

      constexpr std::uint64_t lim_u32 =
        static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max());

      std::uint64_t i64 = static_cast<std::uint64_t>(n);

      // 64-bit part: i64 > lim_u32, so we need the more expensive
      // rng.generateInt64:
      for (; i64 > lim_u32; --i64) {
        std::uint64_t j64 = rng.generateInt64(i64); // [0, i64)
        const std::size_t a = static_cast<std::size_t>(j64);
        const std::size_t b = static_cast<std::size_t>(i64 - 1);
        swapIJ(a, b);
      }

      // 32-bit part: the index now fits in uint32_t so we can use the cheaper
      // rng.generateInt32:
      nc_assert(i64 <= lim_u32);
      std::uint32_t i32 = static_cast<std::uint32_t>(i64);
      for (; i32 > 1; --i32) {
        std::uint32_t j32 = rng.generateInt(i32); // [0, i32)
        swapIJ(static_cast<std::size_t>(j32),
               static_cast<std::size_t>(i32 - 1));
      }
    }

    namespace details {
      template<class V>
      void checkAllSizesEq(std::size_t n, const V& v) {
        if (v.size() != n)
          NCRYSTAL_THROW(BadInput,"randShuffleMultiple: size mismatch");
      }
      template<class V0, class... Vs>
      void checkAllSizesEq(std::size_t n, const V0& v0, const Vs&... vs) {
        if (v0.size() != n)
          NCRYSTAL_THROW(BadInput,"randShuffleMultiple: size mismatch");
        checkAllSizesEq(n, vs...);
      }
      template<class... TVectors>
      void randShuffleMultipleImpl(RNG& rng, std::size_t n, TVectors&... vecs) {
        checkAllSizesEq(n, vecs...);
        randShuffleFisherYates(rng, n, [&](std::size_t a, std::size_t b) {
          using std::swap;
          int dummy[] = {0, (swap(vecs[a], vecs[b]), 0)...};
          (void)dummy;
        });
      }
    }
    // Shuffle one or more containers in lockstep (must have same size):
    template<class TVFirst, class... TVRest>
    void randShuffle(RNG& rng, TVFirst& first, TVRest&... rest) {
      details::randShuffleMultipleImpl(rng, first.size(), first, rest...);
    }

    //If CellData fields a1,a2,b1,b2 have already been set, this can be used to
    //finish up, based on lower corner alpha/beta indices.
    using CellData = SABUtils::CellData;
    void initCellDataSvals( CellData& c, const VectD& sab,
                            std::size_t na, std::size_t ia, std::size_t ib )
    {
      nc_assert( c.a2 > c.a1 );
      nc_assert( c.a1 >= 0.0 );
      nc_assert( c.b2 > c.b1 );
      std::size_t i = na*ib+ia;
      c.S[0] = vectAt(sab,i++);
      c.S[1] = vectAt(sab,i);
      i += (na-1);
      c.S[2] = vectAt(sab,i++);
      c.S[3] = vectAt(sab,i);
      for ( int j = 0; j < 4; ++j )
        c.logS[j] = ( c.S[j] > 0.0 ? std::log(c.S[j]) : 0.0 );
    }

    //This can initialise a CellData object from scratch based on lower corner
    //alpha/beta indices.
    CellData createCellData( std::size_t ia, std::size_t ib, const SABData& sd )
    {
      std::size_t na = sd.alphaGrid().size();
      CellData c;
      c.a1 = vectAt(sd.alphaGrid(),ia);
      c.a2 = vectAt(sd.alphaGrid(),ia+1);
      c.b1 = vectAt(sd.betaGrid(),ib);
      c.b2 = vectAt(sd.betaGrid(),ib+1);
      initCellDataSvals( c, sd.sab(), na, ia, ib );
      return c;
    }
    static constexpr std::int64_t packfactor = 1048576;
    static_assert( NC::ncconstexpr_ispow2( packfactor ), "" );
    static_assert(sizeof(std::size_t)>=8,"");
  }
}

namespace NCRYSTAL_NAMESPACE {

  namespace SABRef {
    namespace {

      struct EPtContribInfo {
        using CCells = std::vector<std::pair<double,std::int64_t>>;
        CCells contrib_cells;//{ contrib, cellidx }
        VectD cumul;

        struct Decoded { bool fully_inside; CellData cellData; };
        Decoded decodeCell( const SABData& sab, std::size_t i ) const
        {
          std::int64_t cellidx = vectAt(contrib_cells,i).second;
          Decoded res;
          if ( cellidx < 0 ) {
            res.fully_inside = true;
            cellidx = -cellidx;
          } else {
            res.fully_inside = false;
          }
          nc_assert(cellidx>=1);
          --cellidx;//remove +1, left with packfactor*ib+ia
          nc_assert(cellidx>=0);
          std::size_t ia = static_cast<std::size_t>(cellidx % packfactor);
          std::size_t ib = static_cast<std::size_t>(cellidx / packfactor);
          res.cellData = createCellData( ia, ib, sab );
          return res;
        }
      };

      EPtContribInfo analyseEPt( const SABData& sab, double E_div_kT )
      {
        EPtContribInfo res;

        const auto& alphas = sab.alphaGrid();
        const auto& betas = sab.betaGrid();
        const auto na = alphas.size();
        const auto nb = betas.size();
        nc_assert_always( na < packfactor && nb < packfactor );
        nc_assert_always( na>=2 && nb >= 2 && E_div_kT > 0.0
                          && std::isfinite(E_div_kT) );

        //First figure out which cells have which contributions, and if they are
        //completely covered or not (indicated by negative cellidx):
        res.contrib_cells.reserve(8192);

        const double foure = 4*E_div_kT;
        const double minusE = -E_div_kT;
        auto ptInsideKB = [foure]( double a, double b )
        {
          return ncsquare( a-b ) < foure * a;//fixme: <= ?
        };

        for ( auto ia : ncrange(na-1) ) {
          const double a1 = vectAt(alphas,ia);
          const double a2 = vectAt(alphas,ia+1);
          for ( auto ib : ncrange(nb-1) ) {
            const double b2 = vectAt(betas,ib+1);
            if ( b2 <= minusE )
              continue;//no overlap (cheap check)
            const double b1 = vectAt(betas,ib);
            //Check if cell is fully above beta+(alpha):
            const double wb2 = std::sqrt( foure * a2 );
            if ( b1 >= a2 + wb2 )
              continue;
            //Check against beta-(alpha) is a bit more complicated:
            if ( a1 >= E_div_kT && b2 <= a1 - std::sqrt( foure * a1 ) )
              continue;
            if ( a2 <= E_div_kT && b2 <= a2 - wb2 )
              continue;

            //Ok, we have an overlap. Let us check that S is not zero in entire
            //cell:
            CellData cell;
            cell.a1 = a1;
            cell.a2 = a2;
            cell.b1 = b1;
            cell.b2 = b2;
            initCellDataSvals( cell, sab.sab(), na, ia, ib );
            if ( cell.S[0]+cell.S[1]+cell.S[2]+cell.S[3]==0.0 )
              continue;//no actual contribution since S=0 in entire cell

            //We next need to determine if the cell is fully inside the
            //kinematic boundary or not. We simply check if the corners are
            //outside (and we only need to check 3 corners, since (a2,b2) will
            //be inside if the other 3 corners are):

            const bool inside11 = ptInsideKB( a1, b1 );
            const bool inside12 = ptInsideKB( a1, b2 );
            const bool inside21 = ptInsideKB( a2, b1 );
            //const bool inside22 = ptInsideKB( a2, b2 );

            const std::int64_t cellidx = packfactor*ib+ia + 1;//+1 so the sign
                                                              //of cellidx can
                                                              //encode info
            nc_assert( cellidx >= 1 );
            if ( inside11 & inside12 & inside21/* & inside22*/ ) {
              //cell is fully inside boundary
              nc_assert(ptInsideKB( a2, b2 ));//fourth corner is also inside
              StableSumKahan sum;
              SABUtils::StdLogLinCellIntegrator::integrateFullCell( cell, sum );
              res.contrib_cells.emplace_back(sum.sum(),-cellidx);//-cellidx
            } else {
              //cell crosses boundary
              StableSumKahan sum;
              SABUtils::StdLogLinCellIntegrator
                ::integrateWithinKB( cell, E_div_kT,
                                     SABCfg::IntegrationScheme::MaxPrec, sum );
              res.contrib_cells.emplace_back(sum.sum(),cellidx);//+cellidx
            }
          }
        }
        //Now sort by contribution, small to large. That will make our
        //cumulative sum more robust (otherwise tiny contributions after large
        //ones might get imprecise).
        std::sort( res.contrib_cells.begin(), res.contrib_cells.end() );

        //Create cumul vect:
        res.cumul.reserve(res.contrib_cells.size());
        StableSum contribsum;
        for ( auto& e : res.contrib_cells ) {
          nc_assert( e.first >= 0.0 );
          contribsum.add( e.first );
          res.cumul.push_back( contribsum.sum() );
        }
        return res;
      }


    }
  }
}

std::pair<NC::VectD,NC::VectD>
NC::SABRef::refSampleAlphaBeta( RNG& rng,
                                const SABData& sab,
                                double E_div_kT,
                                std::uint64_t nsample )
{
  auto analysedEpt = analyseEPt( sab, E_div_kT );
  const std::size_t ncells = analysedEpt.contrib_cells.size();

  //Now count how many are sampled in each cell, so we can afterwards do the
  //actual sampling one cell at a time:
  std::vector<std::size_t> cellcount( ncells, 0 );
  for ( std::uint64_t i = 0; i < nsample; ++i )
    ++vectAt( cellcount, pickRandIdxByWeight( rng, analysedEpt.cumul ) );

  //Now do the actual sampling:
  std::pair<NC::VectD,NC::VectD> res;
  auto& a = res.first;
  auto& b = res.second;
  a.reserve( nsample );
  b.reserve( nsample );
  for ( auto i : ncrange( ncells ) ) {
    auto n = vectAt(cellcount,i);
    if (!n)
      continue;//nothing in this cell

    auto decodedCell = analysedEpt.decodeCell( sab, i );
    bool fully_inside = decodedCell.fully_inside;
    auto& cell = decodedCell.cellData;

    //Now sample:
    if ( fully_inside ) {
      //full cell sampling => just use simple alg even for this reference
      //(fixme: we could try to only use RecFellSampler as a check)
      SABUtils::FullCellSampler fc(&cell);
      while ( n ) {
        auto ab = fc.sampleAlphaBeta( rng );
        a.push_back(ab.first);
        b.push_back(ab.second);
        --n;
      }
    } else {
      //crossing => use RefCellSampler
      SABUtils::RefCellSampler rcs( cell, E_div_kT );
      while ( n ) {
        auto ab = rcs.sampleAlphaBeta( rng );
        a.push_back(ab.alpha);
        b.push_back(ab.beta);
        --n;
      }
    }
  }

  //Finally, shuffle results to not get them ordered by cell:
  randShuffle( rng, a, b);
  return res;
}
