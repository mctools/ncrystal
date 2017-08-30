////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCRotMatrix.hh"
#include "NCrystal/NCException.hh"
#include <iostream>

NCrystal::RotMatrix::~RotMatrix()
{
}

NCrystal::RotMatrix::RotMatrix()
  : Matrix()
{
  m_rowcount=3;
  m_colcount=3;
  m_data.resize(9);
}

NCrystal::RotMatrix::RotMatrix(const Matrix &m)
  : Matrix(m)
{
  if (m.nRows()!=3||m.nCols()!=3)
    NCRYSTAL_THROW(BadInput,"Can only convert 3x3 NCMatrix to NCRotMatrix");
}

NCrystal::RotMatrix::RotMatrix(const NCrystal::RotMatrix &m)
  : Matrix(m)
{
}

NCrystal::RotMatrix::RotMatrix(const double *data)
  : Matrix(3,3,data)
{
}

NCrystal::RotMatrix::RotMatrix(const std::vector<double>& data)
  : Matrix(3,3,data)
{
  if (data.size()!=9)
    NCRYSTAL_THROW(BadInput,"Input vector to NCRotMatrix must have size 3*3=9");
}
