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

#include "NCrystal/internal/cfgutils/NCCfgVars.hh"
#include <iostream>
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"

namespace NC = NCrystal;

namespace {

  std::string sanitise_whitespace( std::string s ) {
    std::string s2 = std::move(s);
    NC::strreplace( s2, "\r", "\\r" );
    NC::strreplace( s2, "\n", "\\n" );
    return s2;
  }

  void test_varids() {

    //unsigned enumval_p = 0;
    for (auto& p : { "bad","density","phasechoice","atomdb","dcutoff","dcutoffup","infofactory","temp",
                    "absnfactory", "bkgd", "bragg", "coh_elas", "elas", "incoh_elas", "inelas",
                    "vdoslux", "scatfactory", "dir1", "dir2", "dirtol",
                    "lcaxis", "lcmode", "mos", "mosprec", "sccutoff"} ) {
      auto pstr = NC::StrView(p);
      auto vid = NC::Cfg::varIdFromName( pstr );
      auto str2 = vid.has_value() ? NC::Cfg::varName(vid.value()) : "invalid";
      std::cout<<p<< "  ->  ";
      if ( !vid.has_value() ) {
        std::cout<<"NOTFOUND"<<std::endl;
      } else {
        std::cout<<NC::enumAsInt(vid.value())<< "  ->  "<<str2<<std::endl;
        nc_assert_always( pstr == str2 );
        //nc_assert_always( enumval_p++ == NC::enumAsInt(vid.value()) );
      }
      //nc_assert_always( enumval_p == 25 );
    }
  }

  void test_parsing() {
    using NC::Cfg::CfgManip;

    auto doTest_impl = [](const char*cstr,bool expectbad)
    {
      std::cout<<"\n------> Parsing \""<<sanitise_whitespace(cstr)<<"\":\n";
      try {
        auto l = NC::Cfg::CfgData();
        auto toplvl = CfgManip::applyStrCfg( l, cstr );
        std::cout<<"  => Parsed "<<CfgManip::size(l)<<" entries (+"<<toplvl.size()<<" toplvl entries)\n";
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;
        std::cout<<"  => .set_temp(500)"<<std::endl;
        CfgManip::set_temp(l,NC::Temperature{500});
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;


      } catch ( NC::Error::BadInput& e ) {
        if (expectbad) {
          std::cout<<"  => Got expected ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
          return;
        }
        throw;
      }
      if ( expectbad )
        NCRYSTAL_THROW(LogicError,"Did not fail as expected!!");
    };
    auto doTest = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,false); };
    auto doTestBad = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,true); };

    doTestBad("vdoslux=34");
    doTest("atomdb=Al is Al27");
    doTest("atomdb=Al is :: Al27");
    doTestBad("atomdb2=Al is Al27");
    doTest("; ; ;\t ;atomdb=Al is Al27@B is B10");
    doTest("; ; \r;\t ;atomdb=Al is Al27@B is B10");
    doTest("atomdb = ;");
    doTest("   atomdb =  ; dcutoff  = 0.5   ");
    doTestBad("vdoslux=2;atomdb={Al is Al27@B is B10}");
    doTestBad("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 and lala ; ; ;vdoslux\n\t=34");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux\n\t=3");
    doTestBad("vdoslux=2;atomdb=Al is Al27@B is B10=Al is Al27@B is B10");
    doTestBad("vdoslux=2;atomdb==Al is Al27@B is B10");
    doTestBad("vdoslux=2;atomdb=Al is Al27@B is B10  =");
    //TODO: we could add a lot more tests here...
  }

  void test_parse_and_fill() {
    using NC::Cfg::CfgManip;
    auto doTest_impl = [](const char*cstr,bool expectbad)
    {
      std::cout<<"\n------> Parsing \""<<sanitise_whitespace(cstr)<<"\":\n";
      try {
        auto l = NC::Cfg::CfgData();
        auto toplvl = CfgManip::applyStrCfg(l,cstr);
        std::cout<<"  => Parsed "<<CfgManip::size(l)<<" entries (+"<<toplvl.size()<<" toplvl entries)\n";
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;
        std::cout<<"  => .set_temp(500)"<<std::endl;
        CfgManip::set_temp(l,NC::Temperature{500});
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;
      } catch ( NC::Error::BadInput& e ) {
        if (expectbad) {
          std::cout<<"  => Got expected ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
          return;
        }
        throw;
      }
      if ( expectbad )
        NCRYSTAL_THROW(LogicError,"Did not fail as expected!!");
    };
    auto doTest = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,false); };
    auto doTestBad = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,true); };
    doTest("temp=300.0");
    doTestBad("temp=100.0fahrenheit");
    doTest("temp=100.0000000000000000000000000000000000");
    doTest("temp=0.123456789012345");
    doTest("temp=0.1234567890123456");
    doTest("temp=0.12345678901234567");
    doTest("temp=0.123456789012345678");
    doTest("temp=0.1234567890123456789");
    doTest("temp=0.1234567890123456789000000000000000000000000000000000000000000000000000010");
    doTest("temp=100.123456789123456789123456789");
    //    std::terminate();

    doTestBad("temp=100.123456789123456789123456789e300");
    doTestBad("temp=-2");
    doTest("temp=-1");
    doTestBad("temp=-0");
    doTest("temp=100.0F");
    doTest("temp=300.0K");
    doTest("temp=30.0C");
    doTest("temp=-20.0C");
    doTest("temp= 5 K ");
    doTestBad("temp=-5 K");
    doTestBad("temp=-infC");
    doTestBad("temp=infC");
    doTestBad("temp=inf");
    doTestBad("temp=+inf");
    doTestBad("temp=-inf");
    doTestBad("temp=nanK");
    doTestBad("temp=nan");
    doTest("temp=-400F");
    doTestBad("dcutoff=inf");
    doTestBad("dcutoff=+inf");

    doTest("dcutoff=0");
    doTest("dcutoff=-1");//backwards compat
    doTest("dcutoff=0.4");
    doTest("dcutoff=0.4Aa");
    doTest("dcutoff=0.04nm");
    doTestBad("dcutoff=-1.1");
    doTestBad("dcutoff=-2");
    doTestBad("dcutoff=-0.1");
    doTestBad("dcutoff=1e-4Aa");
    doTest("dcutoff=1e-3Aa");
    doTest("dcutoff=1e5");
    doTestBad("dcutoff=1.1e5");


    doTest("; ; ;\t ;atomdb=Al is Al27@B is B10");
    doTest("; ; \r;\t ;atomdb=Al is Al27@B is B10");
    doTest("atomdb = ;");
    doTest("   atomdb =  ; dcutoff  = 0.5   ");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10  ; ; ;vdoslux\n\t=3");
    doTestBad("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 and ala ; ; ;vdoslux\n\t=34");
    doTestBad("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux\n\t=34");

    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux=3;temp=20.85C;dcutoff=0.6Aa");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux=3;temp=20.85C;dcutoff=0.6Aa;elas=1");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux=3;temp=20.85C;dcutoff=0.6Aa;elas=1;bragg=0");

    //Test smallvector buffer overflow:
    static_assert(NC::Cfg::VarBufVector::nsmall==7,"");
    doTest("temp=20;dcutoff=1;dcutoffup=2;sccutoff=0.5;sans=0;incoh_elas=0;inelas=freegas;absnfactory=stdabs;atomdb=H:is:H2;scatfactory=stdscat");

    {
      NC::Cfg::CfgData cd{};
      std::cout<<"Test C++ interface:"<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";
      std::cout<<"  => .get_infofactory(): "<<CfgManip::get_infofactory(cd)<<std::endl;
      CfgManip::set_infofactory(cd,NC::StrView::make("lala")); std::cout<<"  => .set_infofactory(StrView(\"lala\")) : "<<CfgManip::get_infofactory(cd)<<std::endl;
      CfgManip::set_infofactory(cd,"bla_blu-bla"); std::cout<<"  => .set_infofactory(const char*(\"bla_blu-bla\")) : "<<CfgManip::get_infofactory(cd)<<std::endl;
      CfgManip::set_infofactory_stdstr(cd,std::string("wuhu")); std::cout<<"  => .set_infofactory(std::string(\"wuhu\")) : "<<CfgManip::get_infofactory(cd)<<std::endl;
      std::cout<<"  => .set_mos(1arcmin)"<<std::endl;
      CfgManip::set_mos( cd, NC::MosaicityFWHM{ 1*NC::kArcMin } );
      std::cout<<"  => .get_mos() : "<<CfgManip::get_mos(cd)<<std::endl;

      std::cout<<"  => .set_lcaxis( {0,0,1} )"<<std::endl;
      CfgManip::set_lcaxis( cd, {0,0,1} );
      std::cout<<"  => .get_lcaxis() : "<<CfgManip::get_lcaxis(cd)<<std::endl;
      CfgManip::set_lcaxis( cd, {1,1,1} );
      std::cout<<"  => .get_lcaxis() : "<<CfgManip::get_lcaxis(cd)<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";
      (void)CfgManip::applyStrCfg(cd,"mos=1arcmin");
      std::cout<<"  => .applyStrCfg(\"mos=1arcmin\") : "<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";
      (void)CfgManip::applyStrCfg(cd,"mos= 1 \tarcmin ");
      std::cout<<"  => .applyStrCfg(\"mos= 1 \\tarcmin \") : "<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";

      auto applyAndPrint = [&cd](const char* cfgstr)
      {
        (void)CfgManip::applyStrCfg(cd,cfgstr);
        std::cout<<"  => .applyStrCfg(\""<<cfgstr<<"\") : "<<std::endl;
        std::cout<<"  => cfgstr: "<<cd<<"\n";
      };
      applyAndPrint("dir1=@crys:0,0,1@lab:1,2,3");
      applyAndPrint("dir1=@crys_hkl:4,5,6@lab:1,2,3");
      applyAndPrint(" dir1 = @ crys_hkl : 4 , 5 , 6 @ lab : 1 , \t2 , 3\n ");
      applyAndPrint("dir1=@crys_hkl:4,5,6@lab:1,2,3");

    }
  }


  void test_cfgdata() {

    static constexpr auto name = "vdoslux";
    static constexpr auto var_id = NC::Cfg::constexpr_varIdFromName(name);
    static constexpr auto var_name = NC::Cfg::varName( var_id );
    static_assert(var_id==NC::Cfg::VarId::vdoslux,"");
    static_assert(NC::constexpr_cstrequal(var_name,name),"");
    std::cout<<"\""<<name<<"\" -> "<<NC::enumAsInt(var_id)<<" -> \""<<var_name<<"\""<<std::endl;

    nc_assert_always( !( NC::StrView("temp") < NC::StrView("dcutoff") ) ) ;
    nc_assert_always( NC::StrView("coh") < NC::StrView("dcutoff") ) ;
    nc_assert_always( NC::StrView("coh_elas") < NC::StrView("dcutoff") ) ;
    nc_assert_always( NC::StrView("dcutoff") < NC::StrView("dcutoffup") ) ;
    nc_assert_always( NC::StrView("dcutoffup") < NC::StrView("incoh_elas") ) ;

    // make_varinfo<vardef_coh_elas>(),
    // make_varinfo<vardef_dcutoff>(),
    // make_varinfo<vardef_dcutoffup>(),
    // make_varinfo<vardef_incoh_elas>(),
    // make_varinfo<vardef_infofactory>(),
    // make_varinfo<vardef_temp>(),
    // make_varinfo<vardef_vdoslux>()

    for ( const auto& varinfo : NC::Cfg::varlist ) {
      //std::cout<<"\""<<varinfo.name()<<"\""<<std::endl;
      auto opt_varid = NC::Cfg::varIdFromName( varinfo.nameSV() );
      nc_assert_always(opt_varid.has_value());
      auto varid = opt_varid.value();
      //inline constexpr unsigned constexpr_varName2Id( const char * name, unsigned detail_idx = 0 )
      nc_assert_always( &NC::Cfg::varInfo(varid) == &varinfo );
      std::cout<<"\""<<varinfo.name()<<"\" -> "<<NC::enumAsInt(varid)<<" -> \""<<NC::Cfg::varInfo(varid).name()<<"\""<<std::endl;
    }

    {
      auto opt_varid = NC::Cfg::varIdFromName( "temp" );
      nc_assert_always(opt_varid.has_value());
      const auto & varinfo = NC::Cfg::varInfo(opt_varid.value());
      auto varid = opt_varid.value();
      auto strval = "120F";
      auto buf_temp = varinfo.from_str( varid, strval );
      std::cout<<" setting \""<<varinfo.name()<<"\" to \""<<strval
               <<"\" -> "<<NC::Cfg::vardef_temp::get_val(buf_temp)<<" -> \"";
      varinfo.stream( std::cout, buf_temp );
      std::cout<<"\"\n";
    }
  }
}

int main() {
  test_cfgdata();
  //  std::cout<<sizeof(B)<<std::endl;
  //  std::cout<<sizeof(NC::SmallVector<B,8,NC::SVMode::FASTACCESS>)<<std::endl;
  test_varids();
  test_parsing();
  test_parse_and_fill();
  NC::Cfg::dumpCfgVarList( std::cout, NC::Cfg::CfgVarListMode::TXT_FULL );
  std::cout<<std::endl;
  NC::Cfg::dumpCfgVarList( std::cout, NC::Cfg::CfgVarListMode::TXT_SHORT );
  std::cout<<std::endl;
  NC::Cfg::dumpCfgVarList( std::cout, NC::Cfg::CfgVarListMode::JSON );
  std::cout<<std::endl;
}
