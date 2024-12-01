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

#include <iostream>
#include "NCrystal/internal/utils/NCMatrix.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCVector.hh"
namespace NC = NCrystal;

int main(int , char**)
{
    //test the rref (Reduced Row Echelon Form) function
    double  data[12]  = {1,2,-1,-4,2,3,-1,-11,-2,0,-3,22};
    NC::Matrix matrix(3,4,data);
    std::cout<<matrix<<std::endl;
    matrix.rref();
    std::cout<<matrix<<std::endl;
    double  expected_data[12]  = {1,0,0,-8,0,1,0,1,0,0,1,-2};
    NC::Matrix  expected_matrix(3,4,expected_data);

    std::cout<<matrix<<std::endl;
    if( ! matrix.isEqual( expected_matrix ) )
      return 1;

    //vector transform
    double darray[] = {0.,-4.,3.,-1.,1.,0.,-1.,4.,-2.} ;
    NC::RotMatrix rot(darray);
    NC::Vector ori(1,2,3);
    NC::Vector newVec=rot*ori;
    NC::Vector expected_newvec(1.,1.,1.);
    if( newVec != expected_newvec )
      return 1;

    std::cout << rot*rot << std::endl;

    //sneak in a few vector tests here:
    if (! (  NC::Vector(0,0,1).isOrthogonal(NC::Vector(1,0,0))
             && NC::Vector(0,5e20,0).isOrthogonal(NC::Vector(-1,0,0))
             && !NC::Vector(0,0,0).isOrthogonal(NC::Vector(-1,0,0))
             && !NC::Vector(0,0,0).isOrthogonal(NC::Vector(0,0,0)) ) )
      return 1;

    return 0;
}
