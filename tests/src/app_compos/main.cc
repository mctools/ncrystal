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

#include "NCrystal/misc/NCCompositionUtils.hh"
#include "NCrystal/internal/atomdb/NCAtomDBExtender.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include <iostream>
namespace NC = NCrystal;
namespace NCCU = NCrystal::CompositionUtils;

void test(const std::vector<std::string>& atomdb_lines) {
  std::cout<<"\n===================================================================="<<std::endl;
  for (auto& e: atomdb_lines)
  std::cout<<"ATOMDB line: "<<e<<std::endl;
  std::cout<<"--------------------------------"<<std::endl;

  NC::AtomDBExtender atomdb;
  for (auto& e: atomdb_lines)
    atomdb.addData(e);
  //For convenience: To set up a fake Info::Composition object, we use the
  //AtomDBExtender to create X, and put the components of that X into the
  //Info::Composition object.
  auto x = atomdb.lookupAtomData("X");
  std::cout<<"AtomData description: "<<x->description(false)<<std::endl;
  NC::Info::Composition compos;
  if (!x->isComposite()) {
    compos.emplace_back(NC::Info::CompositionEntry{1.0,NC::IndexedAtomData{x,NC::AtomIndex{0}}});
  } else {
    for (unsigned i = 0; i < x->nComponents(); ++i) {
      auto& cc = x->getComponent(i);
      compos.emplace_back(NC::Info::CompositionEntry{cc.fraction,NC::IndexedAtomData{cc.data,NC::AtomIndex{i}}});
    }
  }

  NCCU::NaturalAbundanceProvider abp = [](unsigned Z)
  {
    //Silly database:
    std::vector<std::pair<unsigned,double>> natab = { {Z*2,0.5},{Z*2+1,0.5}  };
    return natab;
  };

  std::cout<<"Full breakdown:"<<std::endl;
  auto breakdown = NCCU::createFullBreakdown(compos,abp);
  for (const auto& e : breakdown ) {
    std::cout<<"    --->Z="<<e.first<<" "<<NC::elementZToName(e.first)<<std::endl;
    for ( const auto& A_frac : e.second )
      std::cout<<"        A="<<A_frac.first<<" "<<A_frac.second<<std::endl;
    //    continue;
    NCCU::ElementBreakdownLW lw(e);
    nc_assert_always(lw.Z() == e.first);
    nc_assert_always(lw.isNaturalElement() == (e.second.size()==1&&e.second.front().first==0));
    nc_assert_always( (lw.isNaturalElement()&&lw.nIsotopes()==0) || lw.nIsotopes() == e.second.size());
    double totfrac(0.0);
    for ( const auto& A_frac : e.second )
      totfrac += A_frac.second;
    for (unsigned i = 0; i<lw.nIsotopes();++i) {
      nc_assert_always( lw.A(i) == e.second.at(i).first );
      double expectedFrac = e.second.at(i).second/totfrac;
      nc_assert_always( NC::floateq(lw.fraction(i),expectedFrac,1e-15,1e-15) );
    }
  }

  std::cout<<"LWBreakdown as str: "<<NCCU::breakdownToStr(NCCU::createLWBreakdown(compos,abp))<<std::endl;
}

void testcfg(const std::string& cfgstr)
{
  std::cout<<"Loading cfgstr: \""<<cfgstr<<"\""<<std::endl;
  auto info = NC::FactImpl::createInfo(cfgstr);
  auto compos = info->getComposition();
  std::cout <<"  -> ncomponents: "<<compos.size()<<std::endl;
  for ( auto& e : compos )
    std::cout <<"  -> "<<NC::fmtg(e.fraction)<<" * "<<e.atom.data()<<std::endl;
}

int main() {
  test({"X2 is B",  "B is 0.5 B10 0.5 B11", "X2 is 0.1 B 0.9 Al", "X is 0.2 B10 0.5 Al27 0.3 X2"});
  test({"X is Al"});
  test({"X is 0.9 Al 0.1 Al"});
  test({"Al50 50u 1fm 1b 1b","X is 0.5 Al 0.5 Al50"});
  test({"X2 is 0.5 B10 0.5 Xe","X is 0.5 X2 0.5 B"});
  test({"X2 is 0.5 B 0.5 Xe","X is 0.5 X2 0.5 B"});
  test({"X is B10"});
  test({"X is D"});
  testcfg("C_sg227_Diamond.ncmat");
  testcfg("phases<0.7*C_sg227_Diamond.ncmat&0.28*solid::V/6.12gcm3&0.02*Polyethylene_CH2.ncmat;comp=inelas>;temp=380K");
  testcfg("phases<1.0*C_sg227_Diamond.ncmat>;temp=380K");
  testcfg("solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11");
  testcfg("phases<0.5*solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11&0.5*gasmix::BF3/2bar>");
}
