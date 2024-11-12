////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/internal/NCCFileUtils.hh"

namespace NC = NCrystal;

#include <stdexcept>
#include <cstring>

namespace {
  void require_equal( const NC::mcu8str& str, const std::string& content )
  {
    if ( content.size() != str.size || content != str.c_str )
      throw std::runtime_error("Content check failed");
    const std::string convstr = NC::mcu8str_to_stdstring(&str);
    if ( content.size() != convstr.size() || content != convstr )
      throw std::runtime_error("Conversion check failed");
  }
  void require_owns( const NC::mcu8str& str, std::size_t buflen )
  {
    if ( !str.owns_memory )
      throw std::runtime_error("Ownership check failed (should own)");
    if ( buflen != str.buflen )
      throw std::runtime_error("Buflen check failed (owned)");
  }

  void require_not_owns( const NC::mcu8str& str, std::size_t buflen )
  {
    if ( str.owns_memory )
      throw std::runtime_error("Ownership check failed (should not own)");
    if ( buflen != str.buflen )
      throw std::runtime_error("Buflen check failed (not owned)");
  }
  void cleanup( NC::mcu8str& str )
  {
    NC::mcu8str_dealloc( &str );
    require_equal(str,"");
    require_not_owns(str,0);
  }
}

int main() {
  auto s = NC::mcu8str_create_empty();
  require_equal(s,"");
  require_not_owns(s,0);

  auto s_hello = NC::mcu8str_create_from_cstr("hello");
  require_equal(s_hello,"hello");
  require_owns(s_hello,6);

  auto s_hello2 = NC::mcu8str_create_from_cstr("hellohellohello");
  require_equal(s_hello2,"hellohellohello");
  require_owns(s_hello2,16);
  cleanup(s_hello2);

  NC::mcu8str_dealloc( &s_hello );
  require_equal(s_hello,"");
  require_not_owns(s_hello,0);

  NC::mcu8str_append_cstr(&s_hello,"hello again");
  require_equal(s_hello,"hello again");
  require_owns(s_hello,12);
  NC::mcu8str_clear( &s_hello );
  require_equal(s_hello,"");
  require_owns(s_hello,12);

  auto s2 = NC::mcu8str_create( 12 );
  require_equal(s2,"");
  require_owns(s2,12+1);

  char buf[16];
  auto s3 = NC::mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
  require_equal(s3,"");
  require_not_owns(s3,sizeof(buf));
  NC::mcu8str_append_cstr( &s3, "I am 13 chars" );
  require_equal(s3,"I am 13 chars");
  require_not_owns(s3,sizeof(buf));
  NC::mcu8str_append_cstr( &s3, "A" );
  require_equal(s3,"I am 13 charsA");//14 chars, 15 char buf
  require_not_owns(s3,sizeof(buf));
  NC::mcu8str_append_cstr( &s3, "B" );
  require_equal(s3,"I am 13 charsAB");//15 chars, 16 char buf
  require_not_owns(s3,sizeof(buf));
  NC::mcu8str_append_cstr( &s3, "C" );
  require_equal(s3,"I am 13 charsABC");//16 chars, 17 char buf
  require_owns(s3,17);

  NC::mcu8str_swap( &s3, &s_hello);
  require_equal(s3,"");
  require_owns(s3,12);
  require_equal(s_hello,"I am 13 charsABC");//16 chars, 17 char buf
  require_owns(s_hello,17);

  //Self append:
  NC::mcu8str_append( &s_hello, &s_hello );
  require_equal(s_hello,"I am 13 charsABCI am 13 charsABC");
  require_owns(s_hello,33);//2*16+1

  //Self assignment:
  NC::mcu8str_assign( &s_hello, &s_hello );
  require_equal(s_hello,"I am 13 charsABCI am 13 charsABC");
  require_owns(s_hello,33);//2*16+1

  //Assignment (smaller string):
  auto s_yo = NC::mcu8str_create_from_cstr("yo");
  require_equal(s_yo,"yo");
  require_owns(s_yo,3);
  NC::mcu8str_assign( &s_hello, &s_yo );
  require_equal(s_hello,"yo");
  require_owns(s_hello,33);//does not reduce buflen

  //Append (other):
  NC::mcu8str_append( &s_yo, &s_hello );
  require_equal(s_hello,"yo");
  require_owns(s_hello,33);
  require_equal(s_yo,"yoyo");
  require_owns(s_yo,5);
  NC::mcu8str_reserve( &s_yo, 2 );
  require_equal(s_yo,"yoyo");
  require_owns(s_yo,5);
  NC::mcu8str_reserve( &s_yo, 1000 );
  require_equal(s_yo,"yoyo");
  require_owns(s_yo,1000+1);

  char buf2[40];
  auto s4 = NC::mcu8str_create_from_staticbuffer( buf2, sizeof(buf2) );
  require_equal(s4,"");
  require_not_owns(s4,40);
  NC::mcu8str_reserve( &s4, 39 );
  require_equal(s4,"");
  require_not_owns(s4,40);
  cleanup(s4);
  NC::mcu8str_reserve( &s4, 39 );
  require_equal(s4,"");
  require_owns(s4,40);

  char buf3[20];
  auto s5 = NC::mcu8str_create_from_staticbuffer( buf3, sizeof(buf3) );
  require_equal(s5,"");
  require_not_owns(s5,20);
  NC::mcu8str_append_cstr( &s5, "yiha" );
  require_equal(s5,"yiha");
  require_not_owns(s5,20);

  auto s6 = NC::mcu8str_copy(&s5);
  require_equal(s6,"yiha");
  require_owns(s6,5);
  cleanup(s5);
  require_equal(s5,"");
  require_not_owns(s5,0);

  //Careful editing in place is allowed (otherwise we could not pass the c_str
  //to other C API functions):
  s6.c_str[1] = 'a';
  require_equal(s6,"yaha");
  require_owns(s6,5);
  NC::mcu8str_reserve( &s6, 35 );
  require_equal(s6,"yaha");
  require_owns(s6,35+1);
  s6.c_str[2] = 'H';
  s6.c_str[4] = 'b';
  s6.c_str[5] = 'c';
  s6.c_str[6] = '\0';
  NC::mcu8str_update_size(&s6);
  require_equal(s6,"yaHabc");
  require_owns(s6,35+1);

  NC::mcu8str s7 = NC::mcu8str_create_empty();
  require_equal(s7,"");
  require_not_owns(s7,0);
  NC::mcu8str_update_size( &s7 );
  require_equal(s7,"");
  require_not_owns(s7,0);

  NC::mcu8str_dealloc(&s6);
  require_equal(s6,"");
  require_not_owns(s6,0);
  NC::mcu8str_update_size( &s6 );
  require_equal(s6,"");
  require_not_owns(s6,0);

  cleanup(s5);
  cleanup(s5);
  cleanup(s5);
  cleanup(s_hello);
  cleanup(s);
  cleanup(s2);
  cleanup(s3);
  cleanup(s_yo);
  cleanup(s4);
  cleanup(s6);

  NC::mcu8str s_never_cleaned = NC::mcu8str_create_empty();
  char buf4[128];
  auto s_never_cleaned2 = NC::mcu8str_create_from_staticbuffer( buf4, sizeof(buf4) );
  (void)s_never_cleaned;
  (void)s_never_cleaned2;

  //FIXME: test mcu8str_is_ascii + file utilities

  return 0;
}

