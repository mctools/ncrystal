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

#include "NCrystal/internal/atomdb/NCAtomDBExtender.hh"
#include <iostream>
namespace NC = NCrystal;

namespace test{
  void addData(NC::AtomDBExtender&db, const std::string& data)
  {
    std::cout<<"Adding data: \""<< data<<"\""<<std::endl;
    db.addData(data);
  }
  void addBadData(NC::AtomDBExtender&db, const std::string& data)
  {
    std::cout<<"Trying to add invalid data: \""<< data<<"\""<<std::endl;
    bool gotexpectederr(false);
    try {
      db.addData(data);
    } catch ( NC::Error::BadInput &e) {
      std::cout<<"    -> Exception caught: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
      gotexpectederr=true;
    }
    if (!gotexpectederr)
      throw std::runtime_error("ERROR: Did not trigger exception like it should!");
  }
  uint64_t lookup(NC::AtomDBExtender&db, const std::string& lbl)
  {
    auto res = db.lookupAtomData(lbl);
    std::cout << "Looking up \""<<lbl<<"\" gives: "<<*res << " [uid="<<res->getUniqueID().value<<"]"<<std::endl;
    return res->getUniqueID().value;
  }
  void lookupBad(NC::AtomDBExtender&db, const std::string& lbl)
  {
    std::cout<<"Trying to lookup invalid label: \""<< lbl <<"\""<<std::endl;
    bool gotexpectederr(false);
    try {
      db.lookupAtomData(lbl);
    } catch ( NC::Error::BadInput &e) {
      std::cout<<"    -> Exception caught: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
      gotexpectederr=true;
    }
    if (!gotexpectederr)
      throw std::runtime_error("ERROR: Did not trigger exception like it should!");
  }

}

int main (){
  uint64_t uid_c14(1234567);
  {
    NCrystal::AtomDBExtender db;
    test::lookup(db,"Al");
    test::addData(db,"Al 20u -2fm 5b 1b");
    test::lookup(db,"Al");
    test::lookup(db,"H2");
    test::lookup(db,"D");
    test::addData(db,"Og 20u -2fm 5b 1b");
    test::lookup(db,"Og");
    test::addBadData(db,"Xb 20u -2fm 5b 1b");
    test::lookup(db,"Og");
    test::lookupBad(db,"Xb");
    test::lookup(db,"H2");

    test::lookupBad(db,"C14");
    test::addData(db,"C14 30u 15.0e-3fm 5b 1b");//same as below
    uid_c14 = test::lookup(db,"C14");

    std::cout<<"------------------ Trying another DB instance (while original still in scope):"<<std::endl;
    NCrystal::AtomDBExtender db2;
    test::lookup(db2,"Al");
    test::lookup(db2,"H2");
    test::lookupBad(db2,"C14");
    test::addData(db2,"C14 30u 15.0e-3fm 5b 1b");//same as above
    auto uid_c14_b = test::lookup(db2,"C14");
    nc_assert_always(uid_c14==uid_c14_b);
    test::addData(db2,"C14 30u 15.0e-3fm 5000b 1b");//NOT same as above (both can exist in global database!!)
    auto uid_c14_c = test::lookup(db2,"C14");
    nc_assert_always(uid_c14!=uid_c14_c);
  }

  {
    std::cout<<"------------------ Trying another DB instance (original has gone out of scope):"<<std::endl;
    NCrystal::AtomDBExtender db3;
    test::lookupBad(db3,"C14");
    test::addData(db3,"C14 30u 15.0e-3fm 5b 1b");//same as original above
    auto uid_c14_d = test::lookup(db3,"C14");
    nc_assert_always(uid_c14==uid_c14_d);
    std::cout<<"------------------ Clear global cache (but keep local):"<<std::endl;
    NC::AtomDBExtender::clearGlobalCache();
    test::lookup(db3,"C14");
    std::cout<<"------------------ Clear local cache as well:"<<std::endl;
    db3 = NCrystal::AtomDBExtender(false);
    test::lookupBad(db3,"C14");
    test::addData(db3,"C14 30u 15.0e-3fm 5b 1b");//same as original above (but now both local and global caches are empty)
    test::lookup(db3,"C14");
  }

  {
    std::cout<<"------------------ Test mixtures:"<<std::endl;
    NCrystal::AtomDBExtender db;
    test::lookup(db,"Al");
    test::lookup(db,"Cr");
    test::addData(db,"Al is 0.99 Al 0.01 Cr");
    test::lookup(db,"Al");
    test::addData(db,"Al is 0.5 Al 0.5 Cr");
    test::lookup(db,"Al");
    test::addData(db,"B11 1u 2fm 3b 4b");
    test::addBadData(db,"B10 is 0.5 B 0.5 B11");
    test::addData(db,"B10 1u 2fm 3b 400b");
    test::addBadData(db,"B is 0.95 B10 0.06 B11");
    test::addData(db,"B is 0.95 B10 0.05 B11");
  }

}
