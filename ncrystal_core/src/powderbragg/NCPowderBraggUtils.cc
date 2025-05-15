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

#include "NCrystal/internal/powderbragg/NCPowderBraggUtils.hh"

namespace NC = NCrystal;
namespace NCPB = NCrystal::PowderBraggUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace PowderBraggUtils {
    namespace {

      double extractV0MultNAtom( const StructureInfo& si )
      {
        if (!(si.volume>0) || !(si.n_atoms>=1) )
          NCRYSTAL_THROW(BadInput,"Passed structure info object has"
                         " invalid volume or n_atoms fields.");
        nc_assert_always(si.n_atoms>0);
        nc_assert_always(si.volume>0);
        double res = si.volume * si.n_atoms;
        nc_assert_always( std::isfinite(res) );
        return res;
      }
    }
  }
}

void NCPB::sortData( PreparedPowderInputData& data )
{
  auto& v = data.d_fm_list;
  std::sort(v.begin(),v.end(),std::greater<PairDD>());
}

void NCPB::checkData( const PreparedPowderInputData& data )
{
  if ( !std::isfinite(data.v0_times_natoms)
       || !(data.v0_times_natoms>0) )
    NCRYSTAL_THROW2(BadInput,"v0_times_natoms is not a finite positive"
                    " number: "<<fmt(data.v0_times_natoms));
  auto& v = data.d_fm_list;
  if ( !std::is_sorted(std::begin(v), std::end(v), std::greater<PairDD>()))
    NCRYSTAL_THROW2(BadInput,"Received unsorted list of (d,f*m) values");
  for ( auto& e : v ) {
    if ( !std::isfinite(e.first) || !(e.first>0.0) ) {
      NCRYSTAL_THROW2(BadInput,"d-spacing not finite and positive"
                      " in list of (d,f*m) values");
    }
    if ( !std::isfinite(e.second) || !(e.second>0.0) ) {
      NCRYSTAL_THROW2(BadInput,"fm not finite and positive"
                      " in list of (d,f*m) values");
    }
  }
}

NC::PreparedPowderInputData NCPB::prepareData( double v0_times_natoms,
                                               VectDFM&& data )
{
  PreparedPowderInputData res;
  res.v0_times_natoms = v0_times_natoms;
  res.d_fm_list = std::move(data);
  sortData(res);
  checkData(res);
  return res;
}

NC::PreparedPowderInputData NCPB::prepareData( no_init_t )
{
  PreparedPowderInputData res;
  res.v0_times_natoms = 1.0;//dummy value of 1 Aa^3 for empty d_fm_list
  return  res;
}

NC::PreparedPowderInputData NCPB::prepareData( const StructureInfo& si,
                                               VectDFM&& data )
{
  return prepareData( extractV0MultNAtom(si), std::move(data) );
}

NC::PreparedPowderInputData NCPB::prepareData( const Info& info )
{
  if (!info.hasHKLInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks HKL information.");
  if (!info.hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,
                   "Passed Info object lacks Structure information.");
  const auto& hklList = info.hklList();
  VectDFM data;
  data.reserve(hklList.size());

  for ( const auto& hkl : hklList ) {
    double f = hkl.fsquared * hkl.multiplicity;
    if ( f == 0.0 )
      continue;
    if (f<0)
      NCRYSTAL_THROW(CalcError,
                     "Inconsistent data implies negative |F|^2*multiplicity.");
    if (data.empty()||data.back().first!=hkl.dspacing) {
      data.emplace_back(hkl.dspacing,f);
    } else {
      data.back().second += f;
    }
  }
  return prepareData(info.getStructureInfo(),std::move(data));
}
