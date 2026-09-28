#ifndef NCrystal_Romberg_FMA_hh
#define NCrystal_Romberg_FMA_hh

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

#include "NCrystal/internal/utils/NCFMADispatch.hh"
#include "NCrystal/internal/utils/NCRomberg.hh"

//The ordinary and MSVC/x86-only Windows-fast variants of NCRomberg.cc's
//detail_romberg_integrate -- see NCFMADispatch.hh for the macros used
//below and the full mechanism/rationale, and NCRomberg_WINFMA.cc for the
//other half.
//
//Uses NCRYSTAL_FMADISPATCH_DECLARATOR_C (not the plain DECLARATOR): this
//function is already extern "C" + NCRYSTAL_APPLY_C_NAMESPACE-wrapped in
//its ordinary form too, for the unrelated (Apple/Mach-O) rule-5 reason
//documented at its original definition site -- see
//docs/devel_fma_attribute.md. Calling the virtual self->evalFuncMany/
//evalFuncManySum/accept/convergenceError from the /arch:AVX2-compiled twin
//is safe: /arch: only changes what instructions *this* TU's own code may
//emit, not the calling convention used to reach externally-defined
//(possibly differently-compiled) functions, which is exactly what a
//virtual call already is.

namespace NCRYSTAL_NAMESPACE {

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, romberg_integrate,
    ( const NCrystal::Romberg* self, double a, double b )
  )

  //Every R(i,j) below combines two products (a trapezoidal-rule update, or a
  //Richardson-extrapolation step) into a single sum: a plain "a*b+c*d" is
  //exactly the shape a compiler may or may not silently fuse the *second*
  //product into (giving fma(c,d,a*b)), and since either fusion choice is a
  //valid reading of the same source expression, different platforms/flags
  //can legitimately pick different ones -- not just contract-or-not, but
  //*which* product gets fused. Made unambiguous throughout with an explicit
  //std::fma that always fuses the *first* product (the second is computed
  //separately first, then used as fma's addend), so every platform performs
  //the identical sequence of roundings. Every remaining plain expression in
  //this function is either a bare subtraction/exact-power-of-2-multiply (no
  //a*b+c shape for a "fma" clone to silently misfuse) or already explicit
  //std::fma, so the whole function is safe to decorate per
  //doc/devel_fma_attribute.md rule 1. evalFuncMany/evalFuncManySum are
  //virtual and so cannot be decorated themselves (rule: never on a virtual
  //member function); their own std::fma calls remain an ordinary
  //(undispatched) call:
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_romberg_integrate,romberg_integrate)
  ( const NCrystal::Romberg* self, double a, double b )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(romberg_integrate,(self,a,b));
    double h = (b-a);
    double fvals[17];//R(4,4) needs 17 equally spaced evaluations, we do them in one go:
    self->evalFuncMany(&fvals[0], 17, a, h*0.0625);

    //To reduce overhead, we unroll the calculations for R(n,k) up to R(5,5),
    //since they are anyway short enough to carry out before entering the main
    //loop:

    h *= 0.5;
    const double R00 = (fvals[0] + fvals[16])*h;
    const double R10 = std::fma( h, fvals[8], 0.5*R00 );
    const double R11 = std::fma( (4./3.), R10, (-1./3.)*R00 );
    h *= 0.5;
    const double R20 = std::fma( h, fvals[4]+fvals[12], 0.5*R10 );
    const double R21 = std::fma( (4./3.), R20, (-1./3.)*R10 );
    const double R22 = std::fma( (16./15.), R21, (-1./15.)*R11 );
    h *= 0.5;
    const double R30 = std::fma( h, (fvals[2]+fvals[6])+(fvals[10]+fvals[14]), 0.5*R20 );
    const double R31 = std::fma( (4./3.), R30, (-1./3.)*R20 );
    const double R32 = std::fma( (16./15.), R31, (-1./15.)*R21 );
    const double R33 = std::fma( (64./63.), R32, (-1./63.)*R22 );
    h *= 0.5;
    const double R40 = std::fma( h, ((fvals[1]+fvals[3])+(fvals[5]+fvals[7]))+((fvals[9]+fvals[11])+(fvals[13]+fvals[15])), 0.5*R30 );
    const double R41 = std::fma( (4./3.), R40, (-1./3.)*R30 );
    const double R42 = std::fma( (16./15.), R41, (-1./15.)*R31 );
    const double R43 = std::fma( (64./63.), R42, (-1./63.)*R32 );
    const double R44 = std::fma( (256./255.), R43, (-1./255.)*R33 );

    if (self->accept(4,R33,R44,a,b))
      return R44;

    //R(4,4) was not enough, try R(5,5):
    const double c5 = self->evalFuncManySum(16, a+h*0.5, h);
    h *= 0.5;
    const double R50 = std::fma( h, c5, 0.5*R40 );
    const double R51 = std::fma( (4./3.), R50, (-1./3.)*R40 );
    const double R52 = std::fma( (16./15.), R51, (-1./15.)*R41 );
    const double R53 = std::fma( (64./63.), R52, (-1./63.)*R42 );
    const double R54 = std::fma( (256./255.), R53, (-1./255.)*R43 );
    const double R55 = std::fma( (1024./1023.), R54, (-1./1023.)*R44 );

    if (self->accept(5,R44,R55,a,b))
      return R55;

    //Still not accepted. Use generic loop for R(6,6) or higher.

    //Set up cache arrays to keep row data of current and previous rows:
    const unsigned maxlevel = 16;
    double cache1[maxlevel], cache2[maxlevel];
    double *row_prev = &cache1[0], *row = &cache2[0];

    row_prev[0] = R50;
    row_prev[1] = R51;
    row_prev[2] = R52;
    row_prev[3] = R53;
    row_prev[4] = R54;
    row_prev[5] = R55;

    unsigned nj = 16;
    for(unsigned i = 6; i < maxlevel; ++i){
      double hh = h;
      h *= 0.5;
      nj *= 2;
      double c = self->evalFuncManySum(nj, a+h, hh);

      row[0] = std::fma( h, c, 0.5*row_prev[0] ); //R(i,0)

      double n_k = 1.;
      for(unsigned j = 0; j < i; ++j) {
        n_k *= 4.0;
        //extrapolate value for R(i,j):
        row[j+1] = std::fma( n_k, row[j], -row_prev[j] ) / (n_k-1.0);
      }

      if (self->accept(i,row_prev[i-1],row[i],a,b))
        return row[i];

      std::swap(row_prev,row);
    }

    //Did not converge:
    self->convergenceError(a,b);

    return row_prev[maxlevel-1];//convergenceError() did not throw or otherwise die, so return best estimate.
  }

}

#endif
