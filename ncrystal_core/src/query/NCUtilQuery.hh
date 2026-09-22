#ifndef NCrystal_UtilQuery_hh
#define NCrystal_UtilQuery_hh

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

#include "NCrystal/internal/query/NCQuery.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include <chrono>

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //queryimpl_fmadiagnose_cmulXXX: Version with NCRYSTAL_FMADISPATCH_ATTR
    NCRYSTAL_FMADISPATCH_ATTR
    void queryimpl_fmadiagnose_cmulDispatch( double* re, double* im,
                                             const double* cre,
                                             const double* cim,
                                             std::size_t n )
    {
      for ( std::size_t i = 0; i < n; ++i ) {
        double a = re[i], b = im[i], c = cre[i], d = cim[i];
        re[i] = std::fma( a, c, -(b*d) );
        im[i] = std::fma( a, d, b*c );
      }
    }

    //queryimpl_fmadiagnose_cmulXXX: Version without NCRYSTAL_FMADISPATCH_ATTR
    void queryimpl_fmadiagnose_cmulReference( double* re, double* im,
                                              const double* cre,
                                              const double* cim,
                                              std::size_t n )
    {
      for ( std::size_t i = 0; i < n; ++i ) {
        double a = re[i], b = im[i], c = cre[i], d = cim[i];
        re[i] = std::fma( a, c, -(b*d) );
        im[i] = std::fma( a, d, b*c );
      }
    }

    void queryimpl_fmadiagnose( std::ostream& os )
    {
      const std::size_t n = 32768;
      VectD re0(n), im0(n), cre(n), cim(n);
      for ( auto i : ncrange(n) ) {
        const double x = static_cast<double>(i);
        re0[i] = 1.0/(1.0+0.001*x);
        im0[i] = 0.3*std::sin(0.01*x);
        cre[i] = std::cos(0.02*x);
        cim[i] = std::sin(0.02*x);
      }

      VectD reA(re0), imA(im0), reB(re0), imB(im0);
      queryimpl_fmadiagnose_cmulDispatch( reA.data(), imA.data(),
                                          cre.data(), cim.data(), n );
      queryimpl_fmadiagnose_cmulReference( reB.data(), imB.data(),
                                           cre.data(), cim.data(), n );
      const bool correct = ( reA == reB && imA == imB );

      //Best-of-ntrials timing, in ns/element (kept fast: this is a diagnostic
      //query, not a benchmark):
      auto timeIt = []( void(*f)(double*,double*,const double*,
                                 const double*,std::size_t),
                        const VectD& re0_, const VectD& im0_,
                        const VectD& cre_, const VectD& cim_ )
      {
        const std::size_t nn = re0_.size();
        const int ntrials = 5, nrep = 20;
        double best = 1e300;
        for ( auto trial : ncrange(ntrials) ) {
          (void)trial;
          VectD re(re0_), im(im0_);
          const auto t0 = std::chrono::steady_clock::now();
          for ( auto r : ncrange(nrep) ) {
            (void)r;
            f( re.data(), im.data(), cre_.data(), cim_.data(), nn );
          }
          const auto t1 = std::chrono::steady_clock::now();
          volatile double sink = re[0] + im[0];
          (void)sink;
          const double ns
            = std::chrono::duration<double,std::nano>(t1-t0).count();
          best = ncmin( best, ns / ( static_cast<double>(nrep)
                                    * static_cast<double>(nn) ) );
        }
        return best;
      };
      const double tDispatch = timeIt( queryimpl_fmadiagnose_cmulDispatch,
                                       re0, im0, cre, cim );
      const double tReference = timeIt( queryimpl_fmadiagnose_cmulReference,
                                        re0, im0, cre, cim );
      const double speedup = tReference / tDispatch;

      //Whether this is an optimised (NDEBUG) build: in an unoptimised
      //(e.g. simplebuild debug cache) build, target_clones clones do not
      //get the optimisation needed to show their benefit, so callers
      //should not enforce a speedup requirement unless this is true:
#ifdef NDEBUG
      constexpr bool optimised = true;
#else
      constexpr bool optimised = false;
#endif
      streamJSONDictEntry( os, "enabled", bool(NCRYSTAL_FMADISPATCH_ENABLED),
                           JSONDictPos::FIRST );
      streamJSONDictEntry( os, "correct", correct );
      streamJSONDictEntry( os, "optimised", optimised );
      streamJSONDictEntry( os, "speedup", speedup, JSONDictPos::LAST );
    }

    void queryimpl_version( std::ostream& os, const Query& query )
    {
      if ( query.size() != 1 || query.front() != "version" )
        NCRYSTAL_THROW2(BadInput, "Invalid JSON query (no arguments"
                        " should follow key \"version\")");
      int vmajor = NCRYSTAL_VERSION_MAJOR;
      int vminor = NCRYSTAL_VERSION_MINOR;
      int vpatch = NCRYSTAL_VERSION_PATCH;
      int vint = NCRYSTAL_VERSION;
      if ( ncgetenv_bool("FIX_QUERY_VERSION_FOR_TESTS") ) {
        vmajor = 123;
        vminor = 45;
        vpatch = 6;
        vint = 1000000*vmajor+1000*vminor+vpatch;
      }
      std::array<int,3> vlist = { vmajor, vminor, vpatch };
      streamJSONDictEntry( os, "major", vmajor, JSONDictPos::FIRST );
      streamJSONDictEntry( os, "minor", vminor );
      streamJSONDictEntry( os, "patch", vpatch );
      streamJSONDictEntry( os, "int", vint );
      streamJSONDictEntry( os, "list", vlist );
      std::ostringstream vstr;
      vstr << vmajor << '.' << vminor << '.' << vpatch;
      streamJSONDictEntry( os, "str", vstr.str(), JSONDictPos::LAST );
    }

    void queryimpl_util( std::ostream& os, const Query& query )
    {
      auto invalid = [&query](const char * reason){
        std::ostringstream ss;
        ss << "Invalid util query: ";
        streamJSON( ss, query );//trick: use streamJSON for easy format.
        if ( reason )
          ss<< " ("<<reason<<')';
        NCRYSTAL_THROW( BadInput, ss.str() );
      };
      constexpr auto sv_util = StrView::make("util");
      if ( query.size() < 2 || query.front() != sv_util ) {
        invalid(nullptr);
        return;
      }
      //shift off the "util" and key entries:
      const auto& key = query.at(1).trimmed();
      auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
      //auto argstr = [&arg]( std::size_t i ) { return arg(i).to_string(); };
      const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

      constexpr auto sv_list = StrView::make("list");
      constexpr auto sv_wl2ekin = StrView::make("wl2ekin");
      constexpr auto sv_ekin2wl = StrView::make("ekin2wl");
      constexpr auto sv_mathval = StrView::make("mathval");
      constexpr auto sv_fmadiagnose = StrView::make("fmadiagnose");
      if ( key == sv_list ) {
        if ( nargs != 0 )
          invalid("no arguments should come after: [\"util\",\"list\"]");
        os<<"[\"wl2ekin\", \"ekin2wl\", \"mathval\", \"fmadiagnose\"]";
      } else if ( key == sv_fmadiagnose ) {
        if ( nargs != 0 )
          invalid("no arguments should come after:"
                  " [\"util\",\"fmadiagnose\"]");
        queryimpl_fmadiagnose( os );
      } else if ( isOneOf(key,sv_wl2ekin,sv_ekin2wl) ) {
        double val = ( nargs == 1
                       ? arg(0).toDbl().value_or(-1.0)
                       : -1.0 );
        if ( !(val>=0.0) )
          invalid( key==sv_wl2ekin
                   ? "correct usage: [\"util\",\"wl2ekin\",VALUE_ANGSTROM]"
                   : "correct usage: [\"util\",\"ekin2wl\",VALUE_EV]" );
        streamJSON( os, ( key==sv_wl2ekin
                          ? wl2ekin( val )
                          : ekin2wl( val ) ) );
      } else if ( key == sv_mathval ) {
        if ( nargs != 1 || arg(0) != "kpowx" )
          invalid("Only supported mathval arg is for now [\"kpowx\"]");
        VectD testsamplek = { 1e-300, 1e-200, 1e-30, 1e-10, 1e-3, 0.5,
          0.95, 0.99,
          1.0-1.1e-4, 1.0-1e-4, 1.0-0.9e-4, 1.0,
          1.0+1.1e-4, 1.0+1e-4, 1.0+0.9e-4,
          1.01, 2.0, 1e2, 1e5, 1e30, 1e200, 1e300 };
        VectD testsampleR = { 1e-300,1e-20,1e-5,1e-2,0.1,0.5,0.9,
                              0.9899,0.9901,1.0-1e-5,1.0-1e-10,1.0 };
        VectD testsample;
        testsample.reserve( testsamplek.size() * testsampleR.size() );
        for ( auto& k : testsamplek ) {
          for ( auto& R : testsampleR )
            testsample.push_back(randKPowX( k, std::log(k), R ));
        }

        VectD testk = { 1e-300, 1e-10, 1e-3, 0.5, 0.74,0.76,0.779,0.781,0.79,
                        0.81, 0.89, 0.91, 0.95, 0.99, 1.0-1e-5, 1.0-1e-14, 1.0,
                        1.0+1e-14, 1.0+1e-5, 1.01, 1.05, 1.09, 1.11, 1.19,
                        1.21, 1.219,1.221, 1.24,1.26,3.0, 400.0, 1e4, 1e10,
                        1e300 };
        VectD kpowx, xkpowx;
        for ( auto k : testk ) {
          Optional<double> lnk = std::log(k);
          double v1 = integrate01_kpowx( k, lnk );
          double v2 = integrate01_kpowx( k );
          nc_assert_always( v1==v2 );
          kpowx.push_back( v1 );
          v1 = integrate01_xkpowx( k, lnk );
          v2 = integrate01_xkpowx( k );
          nc_assert_always( v1==v2 );
          xkpowx.push_back( v1 );
        }
        os << "{\"k\":";
        streamJSON(os,testk);
        os << ",\"integral01_kpowx\":";
        streamJSON(os,kpowx);
        os << ",\"integral01_xkpowx\":";
        streamJSON(os,xkpowx);
        os << ",\"sample_kpowx\":{";
        os << "\"k\":";
        streamJSON(os,testsamplek);
        os << ",\"R\":";
        streamJSON(os,testsampleR);
        os << ",\"samples\":";
        streamJSON(os,testsample);
        os<<"}}";
      } else {
        invalid(nullptr);
      }

    }
  }
}

#endif
