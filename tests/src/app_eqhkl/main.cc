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

#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include <string>
#include <vector>
#include <cstdlib>//std::abs

int main(int argc,char**argv) {

  int minhkl(-2);
  int maxhkl(3);
  bool show_all_equrefl(false);
  if (argc>=2&&std::string(argv[1])=="detailed") {
    minhkl=-10;
    maxhkl=10;
    show_all_equrefl = true;
  }
  std::vector<int> sgs = {1,3,16,75,89,143,149,168,177,195,207};

  for (unsigned isg = 0; isg<sgs.size(); ++isg) {
    int sg = sgs.at(isg);
    printf("spacegroup %i\n",sg);
    NCrystal::EqRefl eqrefcalc(sg);
    for (int h=minhkl;h<maxhkl;++h)
      for (int k=minhkl;k<maxhkl;++k)
        for (int l=minhkl;l<maxhkl;++l)
          {
            auto r_raw = eqrefcalc.getEquivalentReflections(h,k,l);
            std::set<NCrystal::HKL> r;
            for (auto it = r_raw.begin();it!=r_raw.end();++it) {
              r.insert(*it);
              r.insert(NCrystal::HKL(-it->h,-it->k,-it->l));
            }
            printf("   Neqrfl[%i,%i,%i] = %i\n",h,k,l,(int)r.size());
            bool show_equrefl = show_all_equrefl || (std::abs(h)+std::abs(k)+std::abs(l)<7);
            if (show_equrefl) {
              printf("     ");
              for (auto it = r.begin();it!=r.end();++it)
                printf(" (%i,%i,%i)",it->h,it->k,it->l);
            }
            printf("\n");
          }
  }
}

