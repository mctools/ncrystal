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

#include "NCrystal/internal/NCStrView.hh"
#include <iostream>
namespace NC = NCrystal;

int main()
{
  std::string loremipsum = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";

  auto testimpl = [](bool expectbad, NC::StrView text,const NC::WordWrapCfg& cfg)
  {
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
    try {
      NC::streamWrappedText(std::cout,text,cfg);
    } catch ( NC::Error::BadInput& e ) {
      if ( expectbad ) {
        std::cout<<"Got expected BadInput error: "<<e.what()<<std::endl;
        std::cout<<"-----------------------------------------------------------------------\n"<<std::endl;
        return;
      }
      throw;
    }
    if ( expectbad )
      throw std::runtime_error("test did not fail like expected");
    std::cout<<"-----------------------------------------------------------------------\n"<<std::endl;
  };
  auto test = [&testimpl](NC::StrView text,const NC::WordWrapCfg& cfg = NC::WordWrapCfg()) { testimpl(false,text,cfg); };
  auto testBad = [&testimpl](NC::StrView text,const NC::WordWrapCfg& cfg = NC::WordWrapCfg()) { testimpl(true,text,cfg); };

  test(loremipsum);
  NC::WordWrapCfg cfg20;
  cfg20.colwidth=20;
  NC::WordWrapCfg cfg20_prefix{cfg20};
  auto testprefix = NC::StrView::make("PREFIX$>");
  cfg20_prefix.prefix = testprefix;
  cfg20_prefix.colwidth = 20 + testprefix.size();
  NC::WordWrapCfg cfg20_forbidoverflow{cfg20};
  cfg20_forbidoverflow.overflow_is_error = true;

  NC::WordWrapCfg cfg20_initialoffset{cfg20};
  cfg20_initialoffset.initial_offset = 12;

  NC::WordWrapCfg cfg20_prefix_and_initialoffset{cfg20_prefix};
  cfg20_prefix_and_initialoffset.initial_offset = 4;

  test(loremipsum,cfg20);
  test("012345678912345     \t \n \r   abcd    ",cfg20);
  test("   \n012345678912345     \t \n \r   abcd    ",cfg20);
  test("012345678912345     \t \n \r   abcde   ",cfg20);
  test("012345678912345     \t \n \r   abcd asd   ",cfg20);
  test("01234567891234567890     \t \n \r   abcd asd   ",cfg20);
  test("012345678912345678901     \t \n \r   abcd asd   ",cfg20);
  testBad("012345678912345678901     \t \n \r   abcd asd   ",cfg20_forbidoverflow);
  test(" dfdf df df df 0123456789123456789012345     \t \n \r   abcd asd   ",cfg20);
  test(" dfdf df df df 012345678912345678901     \t \n \r   abcd asd   ",cfg20);
  testBad(" dfdf df df df 012345678912345678901     \t \n \r   abcd asd   ",cfg20_forbidoverflow);

  test(loremipsum,cfg20_prefix);
  test("012345678912345     \t \n \r   abcd    ",cfg20_prefix);
  test("   \n012345678912345     \t \n \r   abcd    ",cfg20_prefix);
  test("012345678912345     \t \n \r   abcde   ",cfg20_prefix);
  test("012345678912345     \t \n \r   abcd asd   ",cfg20_prefix);
  test("01234567891234567890     \t \n \r   abcd asd   ",cfg20_prefix);
  test("012345678912345678901     \t \n \r   abcd asd   ",cfg20_prefix);
  test(" dfdf df df df 0123456789123456789012345     \t \n \r   abcd asd   ",cfg20_prefix);
  test(" dfdf df df df 012345678912345678901     \t \n \r   abcd asd   ",cfg20_prefix);

  test(loremipsum,cfg20_initialoffset);
  test("012345678912345     \t \n \r   abcd    ",cfg20_initialoffset);
  test("   \n012345678912345     \t \n \r   abcd    ",cfg20_initialoffset);
  test("012345678912345     \t \n \r   abcde   ",cfg20_initialoffset);
  test("012345678912345     \t \n \r   abcd asd   ",cfg20_initialoffset);
  test("01234567891234567890     \t \n \r   abcd asd   ",cfg20_initialoffset);
  test("012345678912345678901     \t \n \r   abcd asd   ",cfg20_initialoffset);
  test(" dfdf df df df 0123456789123456789012345     \t \n \r   abcd asd   ",cfg20_initialoffset);
  test(" dfdf df df df 012345678912345678901     \t \n \r   abcd asd   ",cfg20_initialoffset);

  test(loremipsum,cfg20_prefix_and_initialoffset);
  test("012345678912345     \t \n \r   abcd    ",cfg20_prefix_and_initialoffset);
  test("   \n012345678912345     \t \n \r   abcd    ",cfg20_prefix_and_initialoffset);
  test("012345678912345     \t \n \r   abcde   ",cfg20_prefix_and_initialoffset);
  test("012345678912345     \t \n \r   abcd asd   ",cfg20_prefix_and_initialoffset);
  test("01234567891234567890     \t \n \r   abcd asd   ",cfg20_prefix_and_initialoffset);
  test("012345678912345678901     \t \n \r   abcd asd   ",cfg20_prefix_and_initialoffset);
  test(" dfdf df df df 0123456789123456789012345     \t \n \r   abcd asd   ",cfg20_prefix_and_initialoffset);
  test(" dfdf df df df 012345678912345678901     \t \n \r   abcd asd   ",cfg20_prefix_and_initialoffset);

  return 0;
}
