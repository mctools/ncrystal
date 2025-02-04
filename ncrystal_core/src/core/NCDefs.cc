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

  static_assert(std::is_trivially_copyable<ThreeVector>::value,"");
  static_assert(std::is_trivially_destructible<ThreeVector>::value,"");
  static_assert(std::is_nothrow_destructible<ThreeVector>::value,"");
  static_assert(std::is_trivially_constructible<std::array<double,3>>::value,"");
  static_assert(std::is_trivially_constructible<ThreeVector::stdarray_type>::value,"");
  static_assert( ncconstexpr_roundToNextMultipleOf(0,3) == 3, "" );
  static_assert( ncconstexpr_roundToNextMultipleOf(1,3) == 3, "" );
  static_assert( ncconstexpr_roundToNextMultipleOf(3,3) == 3, "" );
  static_assert( ncconstexpr_roundToNextMultipleOf(32,8) == 32, "" );
  static_assert( ncconstexpr_max( 3,7 ) == 7, "" );
  static_assert( ncconstexpr_max( 3,7,2,-3 ) == 7, "" );
  static_assert( ncconstexpr_gcd( 3,7 ) == 1, "" );
  static_assert( ncconstexpr_gcd( 10,15 ) == 5, "" );
  static_assert( ncconstexpr_gcd( 5*7,5*7*11 ) == 5*7, "" );
  static_assert( ncconstexpr_gcd( 2*5*7,5*7*7 ) == 5*7, "" );
  static_assert( ncconstexpr_lcm( 1,3 ) == 3, "" );
  static_assert( ncconstexpr_lcm( 10,15 ) == 30, "" );
  static_assert( ncconstexpr_lcm( 3,3 ) == 3, "" );
  static_assert( ncconstexpr_lcm( 3,7 ) == 21, "" );
  static_assert( ncconstexpr_lcm( 8,8 ) == 8, "" );
  static_assert( ncconstexpr_lcm( 3,7,1 ) == 21, "" );
  static_assert( ncconstexpr_lcm( 3,7,2 ) == 42, "" );
  static_assert( ncconstexpr_lcm( 3,7,14 ) == 42, "" );

  namespace {
    std::atomic<uint64_t>& getGlobalUIDCounter() {
      static std::atomic<uint64_t> s_global_uid_counter(1);
      return s_global_uid_counter;
    }
  }

  //For unit testing:
  NCRYSTAL_API uint64_t detail_peakNextGlobalUIDValue() { return getGlobalUIDCounter().load(); }
}

NCrystal::UniqueID::UniqueID()
  : m_uid(getGlobalUIDCounter()++)
{
}

NCrystal::RNG::~RNG() = default;

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
void NCrystal::RNG::generateMany( std::size_t n, double* tgt )
{
  for ( std::size_t i = 0; i < n; ++i )
    *tgt++ = this->generate();
}

void NCrystal::RNG::generateRandomBits( std::size_t nbytes, uint8_t* data )
{
  static_assert(sizeof(uint8_t)==1,"");
  constexpr auto nlarge_n64bit = 16;
  constexpr auto nlarge = nlarge_n64bit * sizeof(std::uint64_t);
  while ( nbytes > nlarge ) {
    std::uint64_t buf[nlarge/sizeof(std::uint64_t)];
    static_assert(sizeof(buf)==nlarge,"");
    for ( auto i = 0; i < nlarge_n64bit; ++i )
      buf[i] = this->generate64RndmBits();
    std::memcpy(data,buf,nlarge);
    data += nlarge;
    nbytes -= nlarge;
  }
  while ( nbytes >= sizeof(std::uint64_t) ) {
    std::uint64_t buf = this->generate64RndmBits();
    std::memcpy(data,(void*)&buf,sizeof(std::uint64_t));
    data += sizeof(std::uint64_t);
    nbytes -= sizeof(std::uint64_t);
  }
  if ( nbytes ) {
    //small tail:
    std::uint64_t buf = this->generate64RndmBits();
    std::memcpy(data,(void*)&buf,nbytes);
  }
}
#endif
