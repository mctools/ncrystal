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

#include "NCrystal/internal/phys_utils/NCIofQHelper.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

struct NC::IofQHelper::internal_t
{
  VectD Q;
  VectD QtimesIofQ;
  double normFact;
};

NC::IofQHelper::IofQHelper( internal_t data )
  : m_pwdist( std::move(data.Q),std::move(data.QtimesIofQ) ),
    m_ekinMax( NeutronEnergy{ ksq2ekin( ncsquare( 0.5 * m_pwdist.getXVals().back() ) ) } ),
    m_normFact(data.normFact)
{
}

NC::IofQHelper::IofQHelper( const VectD& Q, const VectD& IofQ )
  : IofQHelper([&Q,&IofQ]() -> internal_t
  {
    auto n = Q.size();
    if ( !nc_is_grid(Q) )
      NCRYSTAL_THROW(BadInput,"IofQHelper got invalid Q grid");
    if ( n != IofQ.size() )
      NCRYSTAL_THROW(BadInput,"IofQHelper got IofQ vector of invalid length");
    if ( ! (Q.front()>0.0) )
      NCRYSTAL_THROW(BadInput,"IofQHelper got Q vector whose first element is not >= 0");

    double emax(0.0);
    for ( auto& e : IofQ ) {
      if ( ! ( e>=0.0 ) )
        NCRYSTAL_THROW(BadInput,"IofQHelper: I(Q) values must be >= 0");
      emax = std::max<double>(e,emax);
    }
    if (!(emax>0.0))
      NCRYSTAL_THROW(BadInput,"IofQHelper: I(Q) must have some values >= 0");

    //Shave off excess trailing zeroes:
    while ( n > 2 && !(vectAt(IofQ,n-1)>0.0) && !(vectAt(IofQ,n-2)>0.0) )
      --n;

    internal_t res;
    auto& q = res.Q;
    auto& f = res.QtimesIofQ;
    if ( Q.front() > 0 ) {
      q.reserve( n + 1 );
      f.reserve( n + 1 );
      q.push_back( 0.0 );
      f.push_back( Q.front() * IofQ.front() );
    } else {
      q.reserve( n );
      f.reserve( n );
    }
    for( auto i : ncrange( n ) ) {
      q.push_back( vectAt(Q,i) );
      nc_assert( vectAt(IofQ,i) >= 0.0 );
      f.push_back( vectAt(Q,i) * vectAt(IofQ,i) );
    }

    StableSum sum;
    for ( auto i : ncrange( q.size()-1 ) )
      sum.add( ( vectAt(q,i+1)-vectAt(q,i) ) * ( vectAt(f,i+1)+vectAt(f,i) ) );
    res.normFact = 0.5 * sum.sum();
    return res;
  }())
{
}
