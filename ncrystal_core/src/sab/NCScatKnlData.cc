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

#include "NCrystal/internal/sab/NCScatKnlData.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

void NC::validateScatKnlData( const NC::ScatKnlDataView& data )
{
  auto xlabel = (data.knltype == ScatKnlData::KnlType::SQW?"Q":"alpha");
  auto ylabel = (data.knltype == ScatKnlData::KnlType::SQW?"omega":"beta");

  if ( data.knltype == ScatKnlData::KnlType::Unspecified )
    NCRYSTAL_THROW(BadInput,"Scatter kernel data has Unspecified type.");

  if ( !(data.temperature.get()>0.0) )
    NCRYSTAL_THROW(BadInput,"Scatter kernel data has invalid temperature");

  if ( !(data.elementMassAMU.get()>0.0) )
    NCRYSTAL_THROW(BadInput,"Scatter kernel data has invalid elementMass");

  if ( !(data.boundXS.get()>=0.0) ) {
    //NB: We allow ==0.0 since ppl might make species sterile with @ATOMDB
    //syntax (the stdscat factory checks that boundXS>0 before instantiating
    //SABScatter objects).
    NCRYSTAL_THROW(BadInput,"Scatter kernel data has invalid boundXS");
  }

  {
    auto al = std::make_pair(&data.alphaGrid,xlabel);
    auto bl = std::make_pair(&data.betaGrid,ylabel);
    for (const auto& e : {al,bl}) {
      if ( e.first->size() < 5 )
        NCRYSTAL_THROW2(BadInput,"Scatter kernel data has invalid "
                        <<e.second<<" grid (must have at least 5 entries)");
      if ( e.first->size() > 65534 )
        NCRYSTAL_THROW2(BadInput,"Scatter kernel data has invalid "
                        <<e.second<<" grid (must have at most 65534 entries)");
      if ( !nc_is_grid(*e.first) )
        NCRYSTAL_THROW2(BadInput,"Scatter kernel data has invalid "
                        <<e.second<<" grid (must consist of sorted, unique, regular numbers)");
    }
  }

  if (! ( data.alphaGrid.front() > 0.0 ) )
    NCRYSTAL_THROW2(BadInput,"Scatter kernel data has non-positive entries in "<<xlabel<<" grid");

  if ( data.knltype == ScatKnlData::KnlType::SCALED_SYM_SAB ) {
    if ( data.betaGrid.front() != 0.0 )
      NCRYSTAL_THROW2(BadInput,"Scatter kernel data "<<ylabel<<" grid must always start with 0.0 when specified as a symmetric table.");
  } else {
    if (! ( data.betaGrid.front() < 0.0 ) )
      NCRYSTAL_THROW2(BadInput,"Scatter kernel data "<<ylabel<<" grid must always start with a negative entry (if the table was"
                     " symmetric, it could start with 0.0)");
  }

  for (auto&e: data.sab)
    if (ncisnan(e) || ncisinf(e) || e<0.0 )
      NCRYSTAL_THROW(BadInput,"Scatter kernel data has negative or NaN/inf S-values");

  if ( data.alphaGrid.size()*data.betaGrid.size() != data.sab.size() )
    NCRYSTAL_THROW(BadInput,"Scatter kernel data has inconsistent array sizes (table size is not product of grid axis sizes)");

  if ( ! (data.suggestedEmax >= 0.0) )
    NCRYSTAL_THROW(BadInput,"Scatter kernel data has invalid suggestedEmax field (must be >=0.0)");

  if ( data.suggestedEmax > 0.0 ) {
    //The energy associated with the kinematic curve touching (alphamax,betamin)
    //gives us an approximate upper bound on the upper energy value of the table:
    const double bmin = data.betaGrid.front();
    const double amax = data.alphaGrid.back();
    const double EmaxUpperBound = constant_boltzmann*data.temperature.get()*(bmin-amax)*(bmin-amax)/(4*amax);
    if ( data.suggestedEmax > EmaxUpperBound*1.000001 )
      NCRYSTAL_THROW2(BadInput,"Scatter kernel data has suggestedEmax ("<<data.suggestedEmax
                      <<" eV) which is clearly too high (grid ranges implies Emax must be less than "
                      <<EmaxUpperBound<<" eV)");
  }

  data.temperature.validate();
  data.boundXS.validate();
  data.elementMassAMU.validate();
}
