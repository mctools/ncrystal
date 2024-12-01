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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/factories/NCDataSources.hh"
#include "NCrystal/factories/NCFact.hh"
#include "NCrystal/interfaces/NCSCOrientation.hh"
#include "NCrystal/tools/NCDump.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include <iostream>
#include <sstream>
#include "NCrystal/interfaces/NCVersion.hh"

namespace NC = NCrystal;

#define REQUIRE nc_assert_always

void test_strrep(const NC::MatCfg& cfg)
{
  auto a = NC::MatCfg(cfg.toStrCfg()).toStrCfg();
  auto b = cfg.toStrCfg();
  if ( a!=b ) {
    std::cout<<"Matchup failure! MatCfg(cfg.toStrCfg()).toStrCfg() : "<<a<<std::endl;
    std::cout<<"                 cfg.toStrCfg()                    : "<<b<<std::endl;
  }
  REQUIRE(a==b);
}

bool floateq(double a, double b, double rtol=1.0e-6, double atol=1.0e-6)
{
  return std::abs(a-b) <= 0.5 * rtol * (std::abs(a) + std::abs(b)) + atol;
}

#define REQUIREFLTEQ(x,y) nc_assert_always(floateq(x,y))

void test_strrep_and_dump(const NC::MatCfg& cfg) {
  std::cout<< "Dump: ";
  cfg.dump(std::cout);
  test_strrep(cfg);
  NC::MatCfg cfg2(std::string("Al_fake2.ncmat;")+cfg.toStrCfg(false));
  std::string s = cfg.toStrCfg(false);
  std::string s2 = cfg2.toStrCfg(false);
  std::cout << "toStrCfg(false) : \"" <<s<< "\"" << std::endl;
  REQUIRE(s==s2);
}

void test()
{
  NC::MatCfg cfg( "Al_fake.ncmat;temp=77.0" );
  //  REQUIRE( !cfg.ignoredEmbeddedConfig() );
  REQUIRE( !cfg.isSingleCrystal() );
  REQUIRE( !cfg.isLayeredCrystal() );
  //method removed  REQUIRE( cfg.isPolyCrystal() );
  REQUIRE( cfg.get_temp().dbl() == 77.0 );
  test_strrep_and_dump(cfg);

  cfg = NC::MatCfg( "Al_fake.ncmat" );
  REQUIRE( cfg.get_temp().dbl() == -1.0 );//code level default
  cfg.set_temp(190.0);
  REQUIRE( cfg.get_temp().dbl() == 190.0 );
  cfg.applyStrCfg("temp= 500");
  REQUIRE( cfg.get_temp().dbl() == 500.0 );
  //  REQUIRE( !cfg.ignoredEmbeddedConfig() );

  REQUIRE( cfg.get_dcutoff() == 0 );
  cfg.set_dcutoff(0.2);
  REQUIRE( cfg.get_dcutoff() == 0.2 );

  //packfact is now obsolete (The get_packfact function is certainly obsolete
  //and always return 1). But we can use density=0.2x or cfg.scale_density(0.2)
  //which has the same (albeit now commulative) effect..
  REQUIRE( cfg.get_density().type == NC::DensityState::Type::SCALEFACTOR );
  REQUIRE( cfg.get_density().value == 1.0 );
  REQUIRE( cfg.get_packfact() == 1.0 );
  cfg.scale_density(0.2);
  REQUIRE( cfg.get_packfact() == 1.0 );
  REQUIRE( cfg.get_density().type == NC::DensityState::Type::SCALEFACTOR );
  REQUIRE( cfg.get_density().value == 0.2 );

  REQUIRE( cfg.get_coh_elas() == true );
  cfg.set_coh_elas(false);
  REQUIRE( cfg.get_coh_elas() == false );
  cfg.applyStrCfg("coh_elas=1");
  REQUIRE( cfg.get_coh_elas() == true );
  cfg.applyStrCfg("bragg=false");
  REQUIRE( cfg.get_coh_elas() == false );
  cfg.applyStrCfg("bragg=true");
  REQUIRE( cfg.get_coh_elas() == true );

  REQUIRE( cfg.get_incoh_elas() == true );
  cfg.set_incoh_elas(false);
  REQUIRE( cfg.get_incoh_elas() == false );
  cfg.applyStrCfg("incoh_elas=1");
  REQUIRE( cfg.get_incoh_elas() == true );

  cfg.applyStrCfg("elas=0");
  REQUIRE( cfg.get_coh_elas() == false );
  REQUIRE( cfg.get_incoh_elas() == false );
  cfg.applyStrCfg("elas=1");
  REQUIRE( cfg.get_coh_elas() == true );
  REQUIRE( cfg.get_incoh_elas() == true );
  cfg.applyStrCfg("elas=0");
  cfg.applyStrCfg("bragg=1");
  REQUIRE( cfg.get_coh_elas() == true );
  REQUIRE( cfg.get_incoh_elas() == false );

  REQUIRE( cfg.get_inelas() == "auto" );
  cfg.set_inelas("0"); REQUIRE( cfg.get_inelas() == "0" );
  cfg.set_inelas("false"); REQUIRE( cfg.get_inelas() == "0" );
  cfg.set_inelas("sterile"); REQUIRE( cfg.get_inelas() == "0" );
  cfg.set_inelas("none"); REQUIRE( cfg.get_inelas() == "0" );

  cfg.applyStrCfg("elas=1;inelas=auto");
  REQUIRE( cfg.get_coh_elas() == true );
  REQUIRE( cfg.get_incoh_elas() == true );
  REQUIRE( cfg.get_inelas() == "auto" );
  cfg.applyStrCfg("bkgd=0");
  REQUIRE( cfg.get_coh_elas() == true );
  REQUIRE( cfg.get_incoh_elas() == false );
  REQUIRE( cfg.get_inelas() == "0" );
  cfg.applyStrCfg("elas=1;inelas=auto");
  cfg.applyStrCfg("bkgd=none");
  REQUIRE( cfg.get_coh_elas() == true );
  REQUIRE( cfg.get_incoh_elas() == false );
  REQUIRE( cfg.get_inelas() == "0" );
  cfg.applyStrCfg("elas=1;inelas=auto");

  // REQUIRE( cfg.get_bkgd_name()=="bla" );
  // REQUIREFLTEQ( cfg.get_bkgdopt_dbl("hej",17.5) , -0.000011 );
  // REQUIRE( cfg.get_bkgdopt_int("bla",5) == 2 );
  // REQUIREFLTEQ( cfg.get_bkgdopt_dbl("foo",13.4) , 13.4 );
  // REQUIRE( cfg.get_bkgdopt_int("foo",13) == 13 );
  // REQUIRE( cfg.get_bkgdopt_flag("tada") );
  // cfg.set_bkgd("phonondebye");
  // REQUIRE( cfg.get_bkgd_name() == "phonondebye" );
  // REQUIRE( !cfg.get_bkgdopt_flag("elastic") );
  // cfg.set_bkgd("phonondebye:elastic");
  // REQUIRE( cfg.get_bkgd_name() == "phonondebye" );
  // REQUIRE( cfg.get_bkgdopt_flag("elastic") );
  // REQUIRE( !cfg.get_bkgdopt_flag("thermalise") );

  REQUIRE( cfg.getDataType() == "ncmat" );

  REQUIRE( cfg.get_infofactory() == "" );
  cfg.set_infofactory("NCMATFactory");
  REQUIRE( cfg.get_infofactory() == "NCMATFactory" );
  cfg.applyStrCfg("infofactory=");//set to empty string (ending with "=" is only allowed for string parameters)

  cfg.set_infofactory("bla");
  REQUIRE( cfg.get_infofactory()=="bla" );
  // REQUIREFLTEQ( cfg.get_infofactopt_dbl("hej",17.5) , -0.000011 );
  // REQUIRE( cfg.get_infofactopt_int("bla",5) == 2 );
  // REQUIREFLTEQ( cfg.get_infofactopt_dbl("foo",13.4) , 13.4 );
  // REQUIRE( cfg.get_infofactopt_int("foo",13) == 13 );
  // REQUIRE( cfg.get_infofactopt_flag("tada") );
  cfg.set_infofactory("mycoolfact");
  REQUIRE( cfg.get_infofactory() == "mycoolfact" );
  // REQUIRE( !cfg.get_infofactopt_flag("dontbesilly") );
  cfg.set_infofactory("mycoolfact");
  REQUIRE( cfg.get_infofactory() == "mycoolfact" );
  // REQUIRE( cfg.get_infofactopt_flag("dontbesilly") );
  // REQUIRE( !cfg.get_infofactopt_flag("blablafoo") );

  REQUIRE( cfg.get_scatfactory() == "" );
  cfg.set_scatfactory("SomeScatterFact");
  REQUIRE( cfg.get_scatfactory() == "SomeScatterFact" );
  REQUIRE( cfg.get_absnfactory() == "" );
  cfg.set_absnfactory("SomeAbsFact");
  REQUIRE( cfg.get_absnfactory() == "SomeAbsFact" );
  //obsolete: cfg.set_packfact(1.0);//required before can set SC parameters

  cfg.set_mos(NC::MosaicitySigma{0.05});
  REQUIREFLTEQ( NC::MosaicitySigma{cfg.get_mos()}.dbl(), 0.05 );

  cfg.set_mos(NC::MosaicityFWHM{0.02});
  REQUIRE( cfg.get_mos().dbl() == 0.02 );

  REQUIREFLTEQ( cfg.get_dirtol() , 1e-4 );
  cfg.set_dirtol(1e-6);
  REQUIREFLTEQ( cfg.get_dirtol() , 1e-6 );
  cfg.set_dirtol(NC::kPi);
  REQUIREFLTEQ( cfg.get_dirtol() , NC::kPi );

  cfg.applyStrCfg("dir1=@crys:1.0,2.0,3.0@lab:4.0,5.0,6.0");
  cfg.applyStrCfg("dir2=@crys_hkl:1.1,2.0,3.0@lab:4.01,5.0,6.0");

  std::cout<<"toStrCfg: "<<cfg.toStrCfg()<<std::endl;
  {
    bool failed = false;
    try { cfg.toEmbeddableCfg(); } catch (NC::Error::BadInput&) { failed = true;  }
    REQUIRE( failed );
  };
  cfg.set_density(NC::Density{1.2});
  std::cout<<"toStrCfg: "<<cfg.toStrCfg()<<std::endl;
  std::cout<<"toEmbeddableCfg: "<<cfg.toEmbeddableCfg()<<std::endl;
  NCrystal::SCOrientation sco = cfg.createSCOrientation();
  REQUIRE( sco.isComplete() );

  test_strrep_and_dump(cfg);

  cfg = NC::MatCfg( "Al_fake.ncmat" );
  cfg.applyStrCfg("temp=80");
  test_strrep_and_dump(cfg);

  cfg = NC::MatCfg( "Al_fake.ncmat" );
  cfg.applyStrCfg("mos=0.08");
  test_strrep_and_dump(cfg);

  cfg = NC::MatCfg( "Al_fake.ncmat" );
  REQUIRE( !cfg.isLayeredCrystal() );
  REQUIRE( cfg.get_lcmode()== 0 );
  cfg.set_lcaxis({0,0.5,0.666666666666666666666666666666666666666});
  REQUIRE( cfg.isLayeredCrystal() );
  REQUIRE( cfg.get_lcmode() == 0 );
  test_strrep_and_dump(cfg);
  cfg = NC::MatCfg( "Al_fake.ncmat;lcaxis=0,0.5,0.666666666666666666666666666666666666666" );
  REQUIRE( cfg.isLayeredCrystal() );
  test_strrep_and_dump(cfg);
  cfg.set_lcmode(-3);
  REQUIRE( cfg.get_lcmode() == -3 );
  cfg.set_lcmode(10023);
  REQUIRE( cfg.get_lcmode() == 10023 );
  test_strrep_and_dump(cfg);

  cfg = NC::MatCfg( "Al_fake.ncmat" );
  try { cfg.get_mos(); }
  catch (NC::Error::MissingInfo&e) { std::cout<<"NCrystal::"<<e.getTypeName()<< ": "<<e.what() << std::endl; }
  try {
    cfg.get_dir1();
  } catch (NC::Error::MissingInfo&e) { std::cout<<"NCrystal::"<<e.getTypeName()<< ": "<<e.what() << std::endl; }
  cfg.applyStrCfg("mos=0.5");
  REQUIRE(cfg.get_mos()==NC::MosaicityFWHM{0.5});
  test_strrep_and_dump(cfg);

  cfg = NC::MatCfg( "Al_fake.ncmat" );
  test_strrep_and_dump(cfg);
  cfg.set_mos(NC::MosaicityFWHM{0.001});
  NC::HKLPoint c1 = {5,1,1};
  NC::LabAxis l1 = {1,0,0};
  cfg.set_dir1( c1,l1 );
  NC::CrystalAxis c2 = {0,0,1/3.};
  NC::LabAxis l2 = {0,1,1};
  cfg.set_dir2( c2, l2 );
  test_strrep_and_dump(cfg);

  // std::cout<<"A:"<<NC::MatCfg( "Al_sg225.ncmat;atomdb= @ @ X :  :::: is Al  @  @@ @ Al is 0.99 Al 0.0100000000000000000000000000000000001  ::Cr"
  //                              "@Be10:0.1e-15u:1fm:1b:1e+5b").get_atomdb()<<std::endl;
  // std::cout<<"B:"<<"X:is:Al@Al:is:0.99:Al:0.0100000000000000000000000000000000001:Cr@Be10:0.1e-15u:1fm:1b:1e+5b"<<std::endl;

  REQUIRE( NC::MatCfg( "Al_sg225.ncmat;atomdb= @ @ X :  :::: is Al  @  @@ @ Al is 0.99 Al 0.0100000000000000000000000000000000001  ::Cr"
                       "@Be10:0.1e-15u:1fm:1b:1e+5b").get_atomdb()
           == "X:is:Al@Al:is:0.99:Al:0.0100000000000000000000000000000000001:Cr@Be10:0.1e-15u:1fm:1b:1e+5b" );

  //Need UCNModeConstruct since REQUIRE macro gets confused otherwise!
  auto UCNModeConstruct = []( NC::UCNMode::Mode mode , double eucn = NC::UCNMode::default_threshold().dbl() )
  {
    return NC::UCNMode{ mode, NC::NeutronEnergy{eucn} };
  };

  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=refine:2meV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=only:0.02eV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:2e-7" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:200e-9" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=only:200neV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:2.12e-3" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:212.0e-5" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=refine:2.12meV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:2.12e-3eV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:212.0e-5eV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=only:0.123eV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:123e-3eV" ).get_ucnmode().value() <<std::endl;
  std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:1.0" ).get_ucnmode().value() <<std::endl;
  //This gives ugly result "200.00000000000003neV" :   std::cout << NC::MatCfg( "Al_sg225.ncmat;ucnmode=refine:2e-4meV" ).get_ucnmode().value() <<std::endl;

  {
    auto cfg2 = NC::MatCfg( "Al_sg225.ncmat" );
    REQUIRE( !cfg2.get_ucnmode().has_value() );
  }
  {
    auto cfg2 = NC::MatCfg( "Al_sg225.ncmat;ucnmode=refine" );
    REQUIRE( cfg2.get_ucnmode() == UCNModeConstruct( NC::UCNMode::Mode::Refine ) );
    REQUIRE( cfg2.get_ucnmode().has_value() );
  }
  {
    auto cfg2 = NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove" );
    REQUIRE( cfg2.get_ucnmode() == UCNModeConstruct( NC::UCNMode::Mode::Remove ) );
  }
  {
    auto cfg2 = NC::MatCfg( "Al_sg225.ncmat;ucnmode=only" );
    REQUIRE( cfg2.get_ucnmode() == UCNModeConstruct( NC::UCNMode::Mode::Only ) );
  }
  {
    auto cfg2 = NC::MatCfg( "Al_sg225.ncmat;ucnmode=refine:1e-4" );
    REQUIRE( cfg2.get_ucnmode() == UCNModeConstruct( NC::UCNMode::Mode::Refine, 1e-4 ) );
  }
  {
    auto cfg2 = NC::MatCfg( "Al_sg225.ncmat;ucnmode=remove:2meV" );
    REQUIRE( cfg2.get_ucnmode() == UCNModeConstruct( NC::UCNMode::Mode::Remove, 2e-3 ) );
  }
  {
    auto cfg2 = NC::MatCfg( "Al_sg225.ncmat;ucnmode=only:4neV" );
    REQUIRE( cfg2.get_ucnmode() == UCNModeConstruct( NC::UCNMode::Mode::Only, 4e-9 ) );
  }
}

void test_copyonwrite() {
  std::cout<<"======= COW ====================================="<<std::endl;
  {
    std::cout<<"===  1 new, 3 copies:"<<std::endl;
    auto cfg = NC::MatCfg( "Al_fake.ncmat" );
    auto cfg2 = cfg;
    auto cfg3(cfg);
    auto cfg4 = cfg3.clone();
    std::cout<<"===  trigger cow:"<<std::endl;
    cfg2.set_temp(100);
    std::cout<<"===  operations which does not trigger cow:"<<std::endl;
    cfg2.scale_density(0.5);//does not trigger since cfg2's data is not shared
    cfg2.get_temp();
    std::cout<<"===  trigger cow different way:"<<std::endl;
    cfg.applyStrCfg( "bragg=false" );
    std::cout<<"===  end of scope - trigger all destructors:"<<std::endl;
  }
  //Test Impl(*other_impl) constructor:
  const NC::MatCfg cfg( "Al_fake.ncmat;temp=77.0;bragg=false;infofactory=SomeFact;dir1=@crys:1.0,2.0,3.0@lab:4.0,5.0,6.0" );
  cfg.dump(std::cout);
  auto cfg2 = cfg;
  cfg2.scale_density(0.5);
  cfg2.dump(std::cout);
  auto cfg3 = cfg;
  cfg3.dump(std::cout);
  std::cout<<"======= EOCOW ==================================="<<std::endl;
}

void test_units()
{
  std::cout<<"======= UNITS ====================================="<<std::endl;
  std::string c("Al_sg225.ncmat;dcutoff=2.0;");
  //temperatures:
  REQUIREFLTEQ(NC::MatCfg(c+"temp=200").get_temp().dbl(),200);
  REQUIREFLTEQ(NC::MatCfg(c+"temp=200 K").get_temp().dbl(),200);
  REQUIREFLTEQ(NC::MatCfg(c+"temp=200 C ").get_temp().dbl(),273.15+200);
  REQUIREFLTEQ(NC::MatCfg(c+"temp= 32   F ").get_temp().dbl(),273.15);
  REQUIREFLTEQ(NC::MatCfg(c+"temp=212F").get_temp().dbl(),373.15);
  //angles
  REQUIREFLTEQ(NC::MatCfg(c+"mos=0.001").get_mos().dbl(),0.001);
  REQUIREFLTEQ(NC::MatCfg(c+"mos=0.001 rad").get_mos().dbl(),0.001);
  REQUIREFLTEQ(NC::MatCfg(c+"mos=2deg").get_mos().dbl(),2*NC::kDeg);
  REQUIREFLTEQ(NC::MatCfg(c+"mos=0.17arcsec").get_mos().dbl(),0.17*NC::kArcSec);
  REQUIREFLTEQ(NC::MatCfg(c+"mos=0.17arcmin").get_mos().dbl(),0.17*NC::kArcMin);
  REQUIREFLTEQ(NC::MatCfg(c+"dirtol=0.001").get_dirtol(),0.001);
  REQUIREFLTEQ(NC::MatCfg(c+"dirtol=0.001 rad").get_dirtol(),0.001);
  REQUIREFLTEQ(NC::MatCfg(c+"dirtol=2deg").get_dirtol(),2*NC::kDeg);
  REQUIREFLTEQ(NC::MatCfg(c+"dirtol=0.17arcsec").get_dirtol(),0.17*NC::kArcSec);
  REQUIREFLTEQ(NC::MatCfg(c+"dirtol=0.17arcmin").get_dirtol(),0.17*NC::kArcMin);

  //lengths:
  REQUIREFLTEQ(NC::MatCfg(c+"dcutoff=2").get_dcutoff(),2.0);
  REQUIREFLTEQ(NC::MatCfg(c+"dcutoff=2Aa").get_dcutoff(),2.0);
  REQUIREFLTEQ(NC::MatCfg(c+"dcutoff=0.2nm").get_dcutoff(),2.0);
  REQUIREFLTEQ(NC::MatCfg(c+"dcutoff=2e-7mm").get_dcutoff(),2.0);
  REQUIREFLTEQ(NC::MatCfg(c+"dcutoff=2e-8cm").get_dcutoff(),2.0);
  REQUIREFLTEQ(NC::MatCfg(c+"dcutoff=2e-10m").get_dcutoff(),2.0);
  REQUIRE(   NC::MatCfg(c+"dcutoff=2.123456789123456789123456789e-8cm").toStrCfg()
          == NC::MatCfg(c+"   ;  ;;  dcutoff  =   2.123456789123456789123456789e-8\n   cm   \n").toStrCfg() );
  //No longer correct (only use origstrrep if gives shorter str):
  //  REQUIRE(   NC::MatCfg(c+"dcutoff=2.123456789123456789123456789e-8cm").toStrCfg() == "Al_sg225.ncmat;dcutoff=2.123456789123456789123456789e-8cm" );

  //Check that two identical strings actually gives identical results:
  auto cfg_aa = NC::MatCfg("   Al_sg225.ncmat ; temp   =    3393.15 K ");
  auto cfg_bb = NC::MatCfg("Al_sg225.ncmat;temp=3393.15K");
  nc_assert_always( !(cfg_aa<cfg_bb) );
  nc_assert_always( !(cfg_bb<cfg_aa) );

  //Check that temp=(20-epsilon)C and temp=293.15 does not lead to a need for different NCInfo objects:
  auto cfg_a = NC::MatCfg("Al_sg225.ncmat;temp=19.99999999999999C");
  auto cfg_b = NC::MatCfg("Al_sg225.ncmat;temp=293.15K");
  nc_assert_always( !(cfg_a<cfg_b) );
  nc_assert_always( !(cfg_b<cfg_a) );

  //but small differences should be diffent:
  auto cfg_c = NC::MatCfg("Al_sg225.ncmat;temp=20.000000001");
  auto cfg_d = NC::MatCfg("Al_sg225.ncmat;temp=293.15K");
  nc_assert_always( int(cfg_c<cfg_d)+int(cfg_d<cfg_c)==1 );

  std::cout<<"======= EOUNITS ==================================="<<std::endl;
}

void test_density () {
  std::cout<<"======= DENSITY ====================================="<<std::endl;
  auto testdata_impl = [](const char * embeddedcfg,
                          const char * params,
                          bool expectbad,
                          bool applyParamsAfter )
  {
    std::cout<<std::endl;
    std::ostringstream ss;
    ss<<"NCMAT v999\n";
    if ( embeddedcfg && !std::string(embeddedcfg).empty() )
      ss << "# NCRYSTALMATCFG["<<embeddedcfg<<"]\n";
    auto textData = NC::makeSO<const NC::TextData>( NC::RawStrData( std::string(ss.str()) ),
                                                    NC::TextData::DataType{"unknown"} );
    std::cout<<"=> Trying to load";
    if ( !applyParamsAfter && params )
      std::cout<<" (applying params \""<<params<<"\" in constructor)";
    std::cout<<":"<<std::endl;
    for ( auto& l : *textData )
      std::cout<<"   -> "<<l<<std::endl;
    static_assert(std::is_copy_constructible<NC::MatCfg>::value, "");
    static_assert(std::is_copy_assignable<NC::MatCfg>::value, "");
    static_assert(NC::Optional<NC::MatCfg>::has_copy_semantics,"");
    NC::Optional<NC::MatCfg> opt_cfg, opt_cfg2;
    try {
      opt_cfg.emplace( textData, (params&&!applyParamsAfter) ? params : "" );
      static int i = 0;
      if ( i++ %2 == 0 )
        opt_cfg2 = opt_cfg.value();//trigger cowpimpl refcounts>1 once in a while
      if (params&&applyParamsAfter) {
        std::cout<<"=> Now calling .applyStrCfg(\""<<params<<"\")"<<std::endl;
        opt_cfg.value().applyStrCfg(params);
      }
      opt_cfg.value().checkConsistency();
    } catch ( NC::Error::BadInput& e ) {
      if (expectbad) {
        std::cout<<"=> Got expected MatCfg ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
        return;
      }
      throw;
    }
    if ( expectbad )
      NCRYSTAL_THROW(LogicError,"Did not fail as expected!!");
    nc_assert_always(opt_cfg.has_value());
    auto& cfg = opt_cfg.value();
    std::cout<<"=> Went OK! Resulting toStrCfg(): \""<<cfg.toStrCfg()<<"\""<<std::endl;
  };
  auto testdata_impl2 = [&testdata_impl](const char * embeddedcfg,
                                         const char * params,
                                         bool expectbad)
  {
    std::vector<bool> applyAfterChoices = { true };
    if ( params )
      applyAfterChoices.push_back(false);
    for ( auto applyAfter : applyAfterChoices )
      testdata_impl(embeddedcfg,params,expectbad,applyAfter);
  };

  auto testdata = [&testdata_impl2]( const char * embeddedcfg = nullptr,
                                    const char * params = nullptr)
  {
    testdata_impl2(embeddedcfg,params,false);
  };
  auto testdata_bad = [&testdata_impl2]( const char * embeddedcfg = nullptr,
                                         const char * params = nullptr)
  {
    testdata_impl2(embeddedcfg,params,true);
  };


  testdata();
  testdata_bad("dcutoff=1;dcutoffup=0.5");
  testdata_bad("mos=2deg");
  testdata_bad("dir1=@crys:1,1,1@lab:0,0,1");
  testdata_bad("dir2=@crys:1,1,1@lab:0,0,1");
  testdata_bad("dir1=@crys:1,1,1@lab:0,0,1;mos=2deg;dirtol=2deg");
  testdata_bad("mos=2deg");
  testdata_bad("dirtol=2deg");

  testdata_bad("dir1=@crys:1,1,1@lab:0,0,1"
               ";dir2=@crys:1,1,1@lab:0,0,1"
               ";mos=2deg;dirtol=180deg");
  testdata("dir1=@crys:1,1,1@lab:0,0,1"
           ";dir2=@crys:1,0,1@lab:0,1,0"
           ";mos=2deg;dirtol=180deg");
  testdata_bad("dir1=@crys:1,1,1@lab:0,1,0"
               ";dir2=@crys:1,0,1@lab:0,1,0"
               ";mos=2deg;dirtol=180deg");
  testdata_bad("dir1=@crys:1,1,1@lab:0,1,0"
               ";dir2=@crys:1,0,1@lab:0,0.999999999999,0"
               ";mos=2deg;dirtol=180deg");
  testdata_bad("dir1=@crys:1,1,1@lab:0,1,0"
               ";dir2=@crys:2,2,1.9999999@lab:1,0,0"
               ";mos=2deg;dirtol=180deg");




  testdata_bad("density=1.2x");//embedded cfgs can not have scale factors != 1
  testdata_bad("packfact=1.2");
  testdata("density=1x");
  testdata_bad("packfact=1.0");//packfact obsolete
  testdata_bad("","packfact=1.2;packfact=2.0");//packfact obsolete
  testdata_bad("","packfact=1.0;packfact=2.0");//packfact obsolete
  testdata("density=1.23e6gcm3");

  testdata("density=1.23e4perAa3");
  testdata_bad("density=1.23e4peraa3");
  testdata("density=1.23e163gcm3");
  testdata_bad("density=1.0000000001e200gcm3");
  testdata_bad("density=0gcm3");
  testdata_bad("density=0x");
  testdata_bad("density=-0.1gcm3");
  testdata_bad("density=1.23e163gcm3;density=1e40x");
  testdata_bad("density=-1.23gcm3");
  testdata_bad("density=0.00gcm3");
  testdata("density=1.23e6kgm3");
  testdata("density=1.23e163kgm3");
  testdata_bad("density=-1.23kgm3");
  testdata_bad("density=0.00kgm3");

  testdata("","density=1.23kgm3");
  testdata("density=1.23kgm3","density=2.2x");
  testdata_bad("density=1.23kgm3","packfact=2.2");//packfact obsolete
  testdata("density=1.23kgm3","density=1gcm3");
  testdata("density=1.23kgm3;density=2x","density=2x");
  testdata_bad("density=1.23kgm3;packfact=2","density=2x");
  testdata_bad("density=1.23kgm3;density=2x","packfact=2");
  testdata_bad("density=1.23kgm3;packfact=2","packfact=2");
  testdata_bad("density=1.23kgm3;packfact=1","packfact=2");
  testdata_bad("density=1.23kgm3;packfact=2","packfact=1");
  testdata("density=1.23kgm3;density=1.5x");
  testdata("","temp=20.0K;density=2x");
  testdata("","temp=20.0K;density=2x;density=2.0perAa3;density=0.25x");
  testdata_bad("density=2x","density=1gcm3");
  testdata_bad("","temp=20.0K;packfact=2.0");
  testdata_bad("","temp=20.0K;packfact=2.0;density=2.0perAa3;density=0.25x");
  testdata_bad("packfact=2.0","density=1gcm3");

  std::cout<<"======= EODENSITY ====================================="<<std::endl;
}


std::string filterDGCODEPath( std::string s ) {
  NC::strreplace(s, NC::ncgetcwd(), "${PWD}");
  return s;
}

void openBadInput(const std::string& filename, bool doCreateInfo = false )
{
  // auto cfg = NC::MatCfg(filename);//+";ignorefilecfg");
  // //  REQUIRE(cfg.ignoredEmbeddedConfig());
  // REQUIRE(cfg.get_temp().dbl()==-1.0);
  // REQUIRE(cfg.get_dcutoff()==0);
  try
    {
      auto cfg2 = NC::MatCfg(filename);
      if ( doCreateInfo )
        NC::createInfo(cfg2);
    }
  catch (NC::Error::BadInput&e)
    {
      std::cout<<"NCrystal::"<<e.getTypeName()<< ": "<<filterDGCODEPath(e.what()) << std::endl;
      return;
    }
  NCRYSTAL_THROW2(LogicError,"Did not end in BadInput exception as expected (\""<<filename<<"\")")
}

void reginmem(const std::string& filename,const std::ostringstream& content ) {
  NC::registerInMemoryFileData(filename,content.str());
}

void test_filecfg()
{
  using Stream = std::ostringstream;
  {
    Stream ss;
    ss << "#Just some garbage\n"
       << "#File\n"
       << "#with an NCRYSTALMATCFG[temp=2.4;dcutoff=0.6] string\n"
       << "#somewhere.\n";
    reginmem("fakefile.bla",ss);
  }
  auto cfg = NC::MatCfg("fakefile.bla");
  cfg.dump(std::cout);
  REQUIRE(cfg.get_temp().dbl()==2.4);
  REQUIRE(cfg.get_dcutoff()==0.6);
  REQUIRE(cfg.get_dcutoffup()>1.0e300);
  REQUIRE(std::isinf(cfg.get_dcutoffup()));

  REQUIRE(cfg.get_absnfactory()=="");
  REQUIRE(cfg.get_coh_elas()==true);
  REQUIRE(cfg.get_inelas()=="auto");
  REQUIRE(cfg.get_infofactory()=="");
  REQUIREFLTEQ(cfg.get_dirtol(),0.0001);
  REQUIRE(cfg.get_packfact()==1.0);
  REQUIRE(cfg.get_scatfactory()=="");
  cfg = NC::MatCfg("fakefile.bla");//;ignorefilecfg");
  //  REQUIRE( cfg.ignoredEmbeddedConfig() );
  cfg.dump(std::cout);
  //  REQUIRE(cfg.get_temp().dbl()==-1.0);
  //REQUIRE(cfg.get_dcutoff()==0);
  cfg = NC::MatCfg("fakefile.bla;temp=40");//had ignorefilecfg
  //  REQUIRE( cfg.ignoredEmbeddedConfig() );
  cfg.dump(std::cout);
  REQUIRE(cfg.get_temp().dbl()==40);
  //REQUIRE(cfg.get_dcutoff()==0);
  cfg = NC::MatCfg("fakefile.bla;temp=60");
  cfg.dump(std::cout);
  REQUIRE(cfg.get_temp().dbl()==60);
  REQUIRE(cfg.get_dcutoff()==0.6);

  {
    Stream ss;
    ss << "#Same but with some allowed spaces\n"
       << "NCRYSTALMATCFG[ temp = 2.4 ; dcutoff = 0.6 ] .\n";
    reginmem("fakefile.bla",ss);
  }

  cfg = NC::MatCfg("fakefile.bla");
  cfg.dump(std::cout);
  REQUIRE(cfg.get_temp().dbl()==2.4);
  REQUIRE(cfg.get_dcutoff()==0.6);

  {
    Stream ss;
    ss << "#Same but empty\n"
       << "  NCRYSTALMATCFG[] .\n";
    reginmem("fakefile.bla",ss);
  }

  cfg = NC::MatCfg("fakefile.bla");
  cfg.dump(std::cout);
  REQUIRE(cfg.get_temp().dbl()==-1.0);
  REQUIRE(cfg.get_dcutoff()==0);

  {
    Stream ss;
    ss << "#Forbidden space\n"
       << "NCRYSTALMATCFG [temp=4].\n";//note the forbidden space before [
    reginmem("fakefile.bla",ss);
  }

  openBadInput("fakefile.bla");
  {
    Stream ss;
    ss << "#Just some garbage\n"
       << "#File\n"
       << "#with two NCRYSTALMATCFG[temp=2.4;dcutoff=0.6] strings\n"
       << "#which is illegal NCRYSTALMATCFG\n";
    reginmem("fakefile.bla",ss);
  }

  openBadInput("fakefile.bla");
  {
    Stream ss;
    ss << "#Just some garbage\n"
       << "#File\n"
       << "#with two NCRYSTALMATCFG[temp=2.4;dcutoff=0.6] strings on the same line NCRYSTALMATCFG\n"
       << "#which is illegal.\n";
    reginmem("fakefile.bla",ss);
  }
  openBadInput("fakefile.bla");


  //fixed @TEMPERATURE should give error if trying to override via temp keyword.
  {
    Stream ss;
    ss << "NCMAT v7\n";
    ss << "@STATEOFMATTER\n";
    ss << "  gas\n";
    ss << "@DYNINFO\n";
    ss << "  element  H\n";
    ss << "  fraction 1\n";
    ss << "  type     sterile\n";
    ss << "@DENSITY\n";
    ss << "  1 kg_per_m3\n";
    ss << "@TEMPERATURE\n";
    ss << "  400.0\n";
    reginmem("fakefile.bla",ss);
  }
  {
    NC::MatCfg ok_cfg("fakefile.bla");
    REQUIRE( NC::createInfo( ok_cfg )->getTemperature().dbl() == 400.0 );
  }
  {
    NC::MatCfg ok_cfg("fakefile.bla;temp=-1");
    REQUIRE( NC::createInfo( ok_cfg )->getTemperature().dbl() == 400.0 );
  }
  {
    NC::MatCfg ok_cfg("fakefile.bla;temp=400K");
    REQUIRE( NC::createInfo( ok_cfg )->getTemperature().dbl() == 400.0 );
  }
  {
    NC::MatCfg ok_cfg("fakefile.bla;temp=399.9K");//No exception because only createInfo triggers the error!!!
  }
  openBadInput("fakefile.bla;temp=399.9K",true);//true means createInfo => trigger error.

  try
    {
      auto cfg2 = NC::MatCfg("doesnotexists.ncmat");
    }
  catch (NC::Error::FileNotFound&e)
    {
      std::cout<<"NCrystal::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
    }

  //ok file, but ignorefilecfg comes too late (update: ignorefilecfg is now obsolete):
  try
    {
      auto cfg2 = NC::MatCfg("Al_fake.ncmat;temp=4;ignorefilecfg");
    }
  catch (NC::Error::BadInput&e)
    {
      std::cout<<"NCrystal::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
    }


  REQUIRE( NC::createScatter("Al_sg225.ncmat;bragg=0;comp=bragg").isNull() );
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp=bragg");
    REQUIRE( cfg2.get_coh_elas() );
    REQUIRE( !cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="0" );
    REQUIRE( !cfg2.get_sans() );
  }
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp=coh_elas");
    REQUIRE( cfg2.get_coh_elas() );
    REQUIRE( !cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="0" );
    REQUIRE( !cfg2.get_sans() );
  }
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp=bragg,sans");
    REQUIRE( cfg2.get_coh_elas() );
    REQUIRE( !cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="0" );
    REQUIRE( cfg2.get_sans() );
  }
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp =   bragg  ,   bragg  ,,,,sans");
    REQUIRE( cfg2.get_coh_elas() );
    REQUIRE( !cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="0" );
    REQUIRE( cfg2.get_sans() );
  }
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp=elas");
    REQUIRE( cfg2.get_coh_elas() );
    REQUIRE( cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="0" );
    REQUIRE( !cfg2.get_sans() );
  }
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp=inelas");
    REQUIRE( !cfg2.get_coh_elas() );
    REQUIRE( !cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="auto" );
    REQUIRE( !cfg2.get_sans() );
  }
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp=incoh_elas");
    REQUIRE( !cfg2.get_coh_elas() );
    REQUIRE( cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="0" );
    REQUIRE( !cfg2.get_sans() );
  }
  {
    auto cfg2 = NC::MatCfg("Al_sg225.ncmat;comp=sans");
    REQUIRE( !cfg2.get_coh_elas() );
    REQUIRE( !cfg2.get_incoh_elas() );
    REQUIRE( cfg2.get_inelas()=="0" );
    REQUIRE( cfg2.get_sans() );
  }
}

void test_atomdb() {
  static const std::string hegasncmat("NCMAT v3\n@DYNINFO\nelement He\nfraction 1\ntype freegas\n@DENSITY\n0.178577 kg_per_m3\n");
  static const std::string hegasncmat_atomdb = hegasncmat + "@ATOMDB\nHe is He3\n";
  static const std::string hegasncmat_inconsistentatomdb = hegasncmat + "@ATOMDB\n\n   \nnodefaults\nHe is He3\n";
  static const std::string hegasncmat_atomdbsyntaxerror = hegasncmat + "@ATOMDB\nHe is He3\nnodefaults\n";

  auto testOK = [](const std::string& ncmatdata, const std::string& cfgstr="",
                   bool expectbadcfg = false, bool expectbadload = false ) {
    std::cout<<"\n---------------------------------------------------------------------------------"<<std::endl;
    std::string fn = "fake.ncmat";
    NC::registerInMemoryFileData(fn,std::string(ncmatdata));
    std::unique_ptr<NC::MatCfg> cfg;
    try {
      cfg = std::make_unique<NC::MatCfg>(fn+";inelas=freegas;"+cfgstr);
    } catch ( NC::Error::BadInput& e ) {
      if (expectbadcfg) {
        std::cout<<"Got expected MatCfg ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
        return;
      }
      throw;
    }
    nc_assert_always(!expectbadcfg);
    nc_assert_always(!!cfg);
    cfg->dump(std::cout);
    try {
      auto info = NC::createInfo(*cfg);
      NC::dump(info);
    } catch ( NC::Error::BadInput& e ) {
      if (expectbadload) {
        std::cout<<"Got expected createInfo ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
        return;
      }
      throw;
    }
    nc_assert_always(!expectbadload);
  };
  auto testBadCfg = [&testOK](const std::string& a, const std::string& b="") { return testOK(a,b,true,false);  };
  auto testBadLoad = [&testOK](const std::string& a, const std::string& b="") { return testOK(a,b,false,true);  };
  testOK(hegasncmat + "@ATOMDB\nHe is He3\n");
  testBadLoad(hegasncmat + "@ATOMDB\n\n   \nnodefaults\n");
  testOK(hegasncmat + "@ATOMDB\n\n   \nnodefaults\nHe 4u 1fm 1b 1b\n");
  testBadLoad(hegasncmat + "@ATOMDB\n\n   \nnodefaults\nHe is He3\n");
  testBadCfg(hegasncmat,"atomdb=He is He3@nodefaults");
  testBadLoad(hegasncmat,"atomdb=nodefaults@He is He3@He3:4u:1fm:1b:1b");
  testOK(hegasncmat,"atomdb=nodefaults@He3:4u:1fm:1b:1b@He is He3");
  testBadLoad(hegasncmat + "@ATOMDB\n\n   \nnodefaults He\nHe 4u 1fm 1b 1b\n");
  testOK(hegasncmat+"@ATOMDB\n\n   \nnodefaults\nHe 4u 1fm 1b 1b\n","atomdb=nodefaults@He3:4u:1fm:1b:1b@He is He3");
  testOK(hegasncmat,"atomdb=nodefaults@He:4u:1fm:1b:1b");
  testBadLoad(hegasncmat,"atomdb=nodefaults");
}

void test_density_multiphase()
{
  std::cout<<"======= DENSITY_MULTIPHASE ====================================="<<std::endl;

  //TODO: in test_density also test comparision operators! Make large number, pass to std::sort, go through and verify ordering.

  auto testcfg_impl = [](const char* cfg, bool expectbad)
  {
    std::cout<<"\n=> Trying to load cfg \""<<cfg<<"\":\n";
    NC::Optional<NC::MatCfg> opt_cfg, opt_cfg2;
    try {
      opt_cfg.emplace( cfg );
      static int i = 0;
      if ( i++ %2 == 0 )
        opt_cfg2 = opt_cfg.value();//trigger cowpimpl refcounts>1 once in a while
      opt_cfg.value().checkConsistency();
    } catch ( NC::Error::BadInput& e ) {
      if (expectbad) {
        std::cout<<"=> Got expected MatCfg ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
        return;
      }
      throw;
    }
    if ( expectbad )
      NCRYSTAL_THROW(LogicError,"Did not fail as expected!!");
    nc_assert_always(opt_cfg.has_value());
    std::cout<<"=> Went OK! Resulting toStrCfg(): \""<<opt_cfg.value().toStrCfg()<<"\""<<std::endl;
  };

  auto testcfg     = [&testcfg_impl](const char* cfg) { return testcfg_impl(cfg,false); };
  auto testcfg_bad = [&testcfg_impl](const char* cfg) { return testcfg_impl(cfg,true); };

  testcfg("phases<0.3*Al_sg225.ncmat&0.7*Ge_sg227.ncmat>");
  testcfg_bad("phases<0.3*Al_sg225.ncmat&0.6*Ge_sg227.ncmat>");

  testcfg("phases<0.3*Al_sg225.ncmat&0.7*Ge_sg227.ncmat>;density=2gcm3");

  testcfg("phases<0.3*Al_sg225.ncmat;density=2gcm3&0.7*Ge_sg227.ncmat>");
  testcfg_bad("phases<0.3*Al_sg225.ncmat;density=2gcm3&0.7*Ge_sg227.ncmat>;packfact=1.4");
  testcfg("phases<0.3*Al_sg225.ncmat;density=2gcm3&0.7*Ge_sg227.ncmat>;density=1.4x");

  testcfg_bad( "Al_sg225.ncmat;ignorefilecfg");//ignorefilecfg obsolete

  testcfg( "Al_sg225.ncmat;sccutoff = 1e-5;density=1.2gcm3" );
  testcfg_bad( "phases<0.2*Al_sg225.ncmat;sccutoff = 1e-5;density=1.2gcm3&0.8*freegas::CO2" );
  testcfg( "phases<0.2*Al_sg225.ncmat;sccutoff = 1e-5;density=1.2gcm3&0.8*freegas::CO2/2kgm3;density=1.5x>" );
  testcfg( "phases<0.2*Al_sg225.ncmat;sccutoff = 1e-5;density=1.2gcm3&0.8*freegas::CO2/2kgm3;density=1.5x>;density=2x" );

  std::cout<<"\n======= EODENSITY_MULTIPHASE ====================================="<<std::endl;

}

void test_misc_leftovers()
{
  //leftover debug code which might as well be used in the unit test:

  std::cout<< NC::MatCfg( "Al_sg225.ncmat;atomdb= @ @ X :  :::: is Al  @  @@ @ Al is 0.99 Al 0.0100000000000000000000000000000000001  ::Cr"
                          "@Be10:0.1e-15u:1fm:1b:1e+5b").get_atomdb()<<std::endl;

  std::cout<<NC::MatCfg( "Al_sg225.ncmat;sccutoff = 1e-5;density=1.2gcm3" )<<std::endl;
  std::cout<<NC::MatCfg( "phases<0.2*Al_sg225.ncmat;sccutoff = 1e-5;density=1.2gcm3&0.8*freegas::CO2/2kgm3;density=1.5x>" )<<std::endl;
  std::cout<<NC::MatCfg( "phases<0.2*Al_sg225.ncmat;sccutoff = 1e-5;density=1.2gcm3&0.8*freegas::CO2/2kgm3;density=1.5x>;density=2x" )<<std::endl;
}

int main(int,char**) {

  std::cout<<"NCrystal::getBuildNameSpace() = " << NCrystal::getBuildNameSpace() << std::endl;
  std::cout<<"NCrystal::getVersion() = " << NCrystal::getVersion() << std::endl;

  NC::registerInMemoryStaticFileData("Al_fake.ncmat","#Just some fake file which exists and has no embedded config.\n");
  NC::registerInMemoryStaticFileData("Al_fake2.ncmat","#Another one.\n");
  NC::registerInMemoryStaticFileData("Al_fake_n_c_m_a_t.bla","NCMAT v1\n#ncmat data but file has weird extension.\n");

  auto td = NC::FactImpl::createTextData("Al_fake_n_c_m_a_t.bla");
  nc_assert_always(td->dataType()=="ncmat");

  {
    auto cfg = NC::MatCfg::createFromRawData( "NCMAT v2\nbla bla", "temp=600K;dcutoff=2");
    std::cout<<cfg<<std::endl;
    std::cout<<cfg.textData();
  }

  test();
  test_copyonwrite();
  test_filecfg();
  test_units();
  test_density();
  test_density_multiphase();
  test_atomdb();//do last-ish, as it also involves NCMAT loading.
  test_misc_leftovers();
  std::cout<<"All tests completed."<<std::endl;

  //For completely pure valgrind:
  fclose( stdin );
  fclose( stdout );
  fclose( stderr );
  return 0;
}
