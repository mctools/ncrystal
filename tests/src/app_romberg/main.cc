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

#include "NCrystal/internal/NCRomberg.hh"
#include "NCrystal/internal/NCMath.hh"
#include <iostream>

namespace NCRYSTAL_NAMESPACE {

  class TestFct : public Romberg {
  private:
  public:
    TestFct() : nevals(0), acc(1e-6) {}
    mutable int nevals;
    double acc;
    virtual bool accept(unsigned level, double prev_estimate, double estimate,double,double) const
    {
      (void)level;//unused in default implementation
      return ncabs(estimate-prev_estimate)<acc;
    }

    virtual double evalFunc(double x) const {
       //       std::cout<<"evalFunc("<<x<<")"<<std::endl;
      ++nevals;
      //sin(x)*cos^3(x)
      double c = std::cos(x);
      return std::sin(x) * c*c*c;
    }
    double integrateExact(double a, double b) const
    {
      //  [cos^4(a)-cos^4(b)]/4
      double ca(std::cos(a)), cb(std::cos(b));
      double ca2(ca*ca), cb2(cb*cb);
      return 0.25 * ( ca2*ca2 - cb2*cb2 );
    }
  };

}

void testint( double a, double b, double acc ) {
  std::cout<<"---- Testing integration over ["<<a<<","<<b<<"] with acc="<<acc<<std::endl;
  NCrystal::TestFct tf;
  tf.acc = acc;
  double res = tf.integrate(a,b);
  double ref = tf.integrateExact(a,b);
  bool bothvanish = NCrystal::ncabs(ref)<1e-16&&NCrystal::ncabs(res)<1e-16;//Avoid tiny fluctuations to spoil unit tests on some platforms (OSX)
  std::cout << "res: "; if (bothvanish) { std::cout<< "<almost zero>" << std::endl; } else { std::cout<< res << std::endl; }
  std::cout << "ref: "; if (bothvanish) { std::cout<< "<almost zero>" << std::endl; } else { std::cout<< ref << std::endl; }
  std::cout << "err: "; if (bothvanish) { std::cout<< "<almost zero>" << std::endl; } else { std::cout<< res-ref << std::endl; }
  std::cout << "relerr: "<< (bothvanish ? 0.0 : (res-ref)/ref ) << std::endl;
  std::cout << "nevals: "<< tf.nevals << std::endl;
}

int main() {
  testint(1,5.4,1);
  testint(1,5.4,1e-3);
  testint(1,2,1e-4);
  testint(0,NCrystal::kPi,1e-4);
  testint(0,3,1e-4);
  testint(0,1.5,1e-4);
  testint(0,1.5,1e-12);
  try {
    testint(0,1.5,-1);
  } catch (NCrystal::Error::CalcError&e) {
    std::cout<<"NCrystal::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
  }

  return 0;
}