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

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/factories/NCDataSources.hh"
#include "NCrystal/interfaces/NCSCOrientation.hh"
#include <iostream>

namespace NC = NCrystal;

void testcfg( const char * cfgstr )
{
  std::cout<<"MatCfg test-case:"<<std::endl;
  std::cout<<"  Input          : "<<cfgstr<<std::endl;
  NC::MatCfg cfg(cfgstr);
  std::cout<<"  Object ostream : "<<cfg<<std::endl;
  std::cout<<"  Object.toStrCfg: "<<cfg.toStrCfg()<<std::endl;

}

void test1()
{
  testcfg( "Al.ncmat;vdoslux=1;dcutoff=0.2");
  testcfg( "Al_400K.ncmat");
  // testcfg( "Al_400K.ncmat;ignorefilecfg");
  testcfg( "phases<0.9999999999999 * Al.ncmat;vdoslux=1>;dcutoff=0.2" );
  testcfg( "phases<0.1 * Al.ncmat;vdoslux=1  & 0.9 * water.ncmat>;dcutoff=0.2" );
  testcfg( "phases<0.1 * Al.ncmat & 0.9 * Al.ncmat  >;dcutoff=0.2" );
  testcfg( "phases<0.1 * Al.ncmat & 0.9 * Al.ncmat;vdoslux=1>;dcutoff=0.2 " );
  testcfg( "phases<1e-10*Al.ncmat&1.0*Al.ncmat\t>;dcutoff=0.2" );


  testcfg( "phases<0.1*Al.ncmat&0.3*Al_400K.ncmat&0.4*water.ncmat&0.2*Al_ec.ncmat>" );
  //  testcfg( "phases<0.1*Al.ncmat&0.3*Al_400K.ncmat;ignorefilecfg&0.4*water.ncmat&0.2*Al_ec.ncmat>" );
  testcfg( "phases<0.1*Al.ncmat;temp=400K&0.3*Al_400K.ncmat&0.4*water.ncmat;temp=400K&0.2*Al_ec.ncmat;temp=400>" );

  bool failed = false;
  auto badcfg = "phases<0.1*Al.ncmat&0.9*phases<1.0*Al_400K.ncmat>>";
  try {
    testcfg(badcfg);
  } catch ( NC::Error::BadInput& e ) {
    failed = true;
    std::cout<<"Failed as expected : "<<e.what()<<std::endl;
  }
  if (!failed)
    NCRYSTAL_THROW2(LogicError,"Cfg string did not fail as expected: "<<badcfg);
}

namespace {

  template <class TSetType,class TGetType,
            TGetType(NC::MatCfg::*get_xxx)() const,
            void(NC::MatCfg::*set_xxx)(TSetType)>
  inline void testAccess(const char* parname, const std::string& parstrval, TSetType valexample1, TSetType valexample2,bool isOrientPar = false) {
    std::cout<<"\nTesting get_"<<parname<<"() and set_"<<parname<<"():"<<std::endl;
    //Construct cfg object where parameter differs on some phases (i.e. get_xxx should fail):
    NC::MatCfg::PhaseList phaselist = { { 0.2, NC::MatCfg("Al.ncmat") },
                                        { 0.75,NC::MatCfg("water.ncmat") },
                                        { 0.05,NC::MatCfg("Be.ncmat") } };
    if ( isOrientPar ) {
      //must always set all or none of mos,dir1,dir2 (and those must be set when
      //dirtol is set). We also increase dirtol to avoid "spurious" errors.
      for (auto& ph : phaselist)
        ph.second.applyStrCfg("mos=0.123arcsec;dir1=@crys:0,0,1@lab:0,0,1;dir2=@crys:0,1,0@lab:0,1,0;dirtol=120deg");
    }

    ((&phaselist.front().second)->*set_xxx)( valexample1 );
    ((&phaselist.back().second)->*set_xxx)( valexample2 );

    NC::MatCfg cfg(std::move(phaselist));
    std::cout<<"  initial cfg: "<<cfg<<std::endl;
    cfg.checkConsistency();
    try {
      std::cout<<"  trying get_"<<parname<<"() with result: ";
      (void)((&cfg)->*get_xxx)();
      nc_assert_always(false);//should not get here
    } catch ( NC::Error::CalcError& e ) {
      std::cout<<"BadInputException["<<e.what()<<"]"<<std::endl;
    }
    ((&cfg)->*set_xxx)(valexample1);
    std::cout<<"  Calling set_"<<parname<<"() changing state to: "<<cfg<<std::endl;
    cfg.checkConsistency();
    std::cout<<"  Calling again get_"<<parname<<"() with result: "<<((&cfg)->*get_xxx)()<<std::endl;
    nc_assert_always(((&cfg)->*get_xxx)() == valexample1 );
    std::string strcfg(parname);
    strcfg+="=";
    strcfg+=parstrval;
    if ( parstrval.find("=") != std::string::npos )
      strcfg = parstrval;
    cfg.applyStrCfg(strcfg);
    std::cout<<"  Applying str cfg \""<<strcfg<<"\" changing state to: "<<cfg<<std::endl;
    cfg.checkConsistency();
    std::cout<<"  Calling again get_"<<parname<<"() with result: "<<((&cfg)->*get_xxx)()<<std::endl;
  }
}

void test2()
{
  std::cout<<"test2-----------------------"<<std::endl;
  NC::MatCfg cfg("phases<0.1*Al.ncmat;temp=400K&0.9*water.ncmat>;dcutoffup=99.0");
  std::cout<<"  cfg = "<<cfg<<std::endl;
  std::cout<<"  cfg.get_inelas() = "<<cfg.get_inelas()<<std::endl;
  std::cout<<"  cfg.get_dcutoff() = "<<cfg.get_dcutoff()<<std::endl;
  std::cout<<"  cfg.get_dcutoffup() = "<<cfg.get_dcutoffup()<<std::endl;
  std::cout<<"  cfg.get_temp() = ";
  try {
    cfg.get_temp();
    nc_assert_always(false);//should not get here
  } catch ( NC::Error::CalcError& e ) {
    std::cout<<"BadInputException["<<e.what()<<"]";
  }
  std::cout<<std::endl;
  std::cout<<"  Calling cfg.set_temp(550)."<<std::endl;
  cfg.set_temp(550);
  std::cout<<"  cfg.get_temp() = "<<cfg.get_temp()<<std::endl;
  testAccess<NC::Temperature,NC::Temperature,&NC::MatCfg::get_temp,&NC::MatCfg::set_temp>("temp","10C",NC::Temperature{210},NC::Temperature{430});
  testAccess<NC::MosaicityFWHM,NC::MosaicityFWHM,&NC::MatCfg::get_mos,&NC::MatCfg::set_mos>("mos","1arcmin",NC::MosaicityFWHM{0.1*NC::kDeg},NC::MosaicityFWHM{1.0*NC::kDeg},true);
  testAccess<double,double,&NC::MatCfg::get_dcutoff,&NC::MatCfg::set_dcutoff>("dcutoff","2Aa",0.9,1.2);
  testAccess<double,double,&NC::MatCfg::get_dcutoffup,&NC::MatCfg::set_dcutoffup>("dcutoffup","10Aa",0.9,1.2);
  //Obsolete: testAccess<double,double,&NC::MatCfg::get_packfact,&NC::MatCfg::set_packfact>("packfact","0.456",0.9,0.3);
  testAccess<double,double,&NC::MatCfg::get_mosprec,&NC::MatCfg::set_mosprec>("mosprec","0.00002",1e-4,1e-5);
  testAccess<double,double,&NC::MatCfg::get_sccutoff,&NC::MatCfg::set_sccutoff>("sccutoff","0.5Aa",0.1,1.0);
  testAccess<double,double,&NC::MatCfg::get_dirtol,&NC::MatCfg::set_dirtol>("dirtol","20deg",1*NC::kDeg,180*NC::kDeg,true);
  testAccess<bool,bool,&NC::MatCfg::get_coh_elas,&NC::MatCfg::set_coh_elas>("coh_elas","true",true,false);
  testAccess<bool,bool,&NC::MatCfg::get_incoh_elas,&NC::MatCfg::set_incoh_elas>("incoh_elas","0",true,false);
  testAccess<const std::string&,std::string,&NC::MatCfg::get_atomdb,&NC::MatCfg::set_atomdb>("atomdb","Al is B","Al:is:Xe","Be:is:Ne");
  testAccess<const std::string&,std::string,&NC::MatCfg::get_inelas,&NC::MatCfg::set_inelas>("inelas","vdosdebye","freegas","sterile");
  testAccess<const std::string&,std::string,&NC::MatCfg::get_infofactory,&NC::MatCfg::set_infofactory>("infofactory","ccc","aaa","bbb");
  testAccess<const std::string&,std::string,&NC::MatCfg::get_scatfactory,&NC::MatCfg::set_scatfactory>("scatfactory","ccc","aaa","bbb");
  testAccess<const std::string&,std::string,&NC::MatCfg::get_absnfactory,&NC::MatCfg::set_absnfactory>("absnfactory","ccc","aaa","bbb");
  testAccess<int,int,&NC::MatCfg::get_lcmode,&NC::MatCfg::set_lcmode>("lcmode","30",1000,-5);
  testAccess<int,int,&NC::MatCfg::get_vdoslux,&NC::MatCfg::set_vdoslux>("vdoslux","5",2,4);
  testAccess<const NC::LCAxis&,const NC::LCAxis&,&NC::MatCfg::get_lcaxis,&NC::MatCfg::set_lcaxis>("lcaxis","0,0.5,0.5",NC::LCAxis(0.0,1.0,0.0),NC::LCAxis(1.0,0.0,0.0));

  using ODir = NC::OrientDir;
  testAccess<const ODir&,ODir,&NC::MatCfg::get_dir1,&NC::MatCfg::set_dir1>("dir1","@crys:0.2,0,0.8@lab:0,0,1",
                                                                           ODir{NC::CrystalAxis{1,2,3},NC::LabAxis{4,5,6}},
                                                                           ODir{NC::HKLPoint{0,1,0},NC::LabAxis{1,0,0}},true);
  testAccess<const ODir&,ODir,&NC::MatCfg::get_dir2,&NC::MatCfg::set_dir2>("dir2","@crys:0.3,0.1,0.8@lab:0,1,1",
                                                                           ODir{NC::CrystalAxis{1,-2,3},NC::LabAxis{2,5,-6}},
                                                                           ODir{NC::HKLPoint{0,1,2},NC::LabAxis{1,-1,0}},true);

  //Obsolete code which was needed before get_dir1/get_dir2 worked independently
  // //Test dir1 and dir2 via SCOrientation interface:
  // NC::SCOrientation sco1, sco2;
  // sco1.setPrimaryDirection(NC::CrystalAxis{0.5,0.5,0.0},NC::LabAxis{1,0,0});
  // sco1.setSecondaryDirection(NC::HKLPoint{0,0,2},NC::LabAxis{1,1,0});
  // sco2.setPrimaryDirection(NC::HKLPoint{1,1,1},NC::LabAxis{0,1,0});
  // sco2.setSecondaryDirection(NC::CrystalAxis{0.5,0.5,0.5},NC::LabAxis{1,0,1});
  // testAccess<const NC::SCOrientation&,NC::SCOrientation,&NC::MatCfg::createSCOrientation,&NC::MatCfg::setOrientation>("SCORIENTDUMMY",
  //                                                                                                                 "dir1=@crys:0.2,0,0.8@lab:0,0,1;dir2=@crys:0,0.6,0.4@lab:0,1,0",
  //                                                                                                                 sco1,sco2,true);
}

int main(int,char**) {

  NC::DataSources::removeAllDataSources();
  NC::DataSources::registerInMemoryStaticFileData("Al.ncmat","#Just some fake file which exists and has no embedded config.\n");
  NC::DataSources::registerInMemoryStaticFileData("Be.ncmat","#Just some fake file which exists and has no embedded config.\n");
  NC::DataSources::registerInMemoryStaticFileData("Al_400K.ncmat","#has embedded cfg: NCRYSTALMATCFG[temp=400K].\n");
  NC::DataSources::registerInMemoryStaticFileData("Al_ec.ncmat","#has embedded cfg: NCRYSTALMATCFG[dcutoff=1.0].\n");
  NC::DataSources::registerInMemoryStaticFileData("water.ncmat","#Another one.\n");

  test1();
  test2();

  //For completely pure valgrind:
  fclose( stdin );
  fclose( stdout );
  fclose( stderr );
  return 0;
}
