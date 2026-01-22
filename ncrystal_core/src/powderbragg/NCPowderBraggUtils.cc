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

#include "NCrystal/internal/powderbragg/NCPowderBraggUtils.hh"

namespace NC = NCrystal;
namespace NCPB = NCrystal::PowderBraggUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace PowderBraggUtils {
    namespace {


      bool cellVolNatomsOK( double volume, unsigned n_atoms )
      {
        return ( volume > 0
                 && n_atoms >= 1
                 && std::isfinite(volume)
                 && std::isfinite(volume*n_atoms) );
      }

      void checkInfoOK( const Info& info )
      {
        if (!info.hasHKLInfo())
          NCRYSTAL_THROW(MissingInfo,
                         "Passed Info object lacks HKL information.");
        if (!info.hasStructureInfo())
          NCRYSTAL_THROW(MissingInfo,
                         "Passed Info object lacks Structure information.");
      }

      void updatePlaneList( PIData::PlaneList& data,
                            const HKLInfo& hkl )
      {
        if ( data.empty()
             || data.back().dsp != hkl.dspacing
             || data.back().fsq != hkl.fsquared ) {
          data.emplace_back();
          data.back().dsp = hkl.dspacing;
          data.back().fsq = hkl.fsquared;
          data.back().mult = hkl.multiplicity;
        } else {
          data.back().mult += hkl.multiplicity;
        }
      }

      void updatePlaneList( PIMergedData::PlaneList& data,
                            const HKLInfo& hkl )
      {
        const double v = hkl.fsquared * hkl.multiplicity;
        nc_assert_always( std::isfinite(v) );
        if ( data.empty() || data.back().dsp != hkl.dspacing ) {
          data.emplace_back();
          data.back().dsp = hkl.dspacing;
          data.back().fsqmult = v;
        } else {
          data.back().fsqmult += v;
        }
      }

      template<class TData>
      TData prepareDataImpl( const HKLList& hklList, const CellData& cell )
      {
        TData res;
        res.cell = cell;
        auto& data = res.planes;
        data.reserve(hklList.size());

        for ( const auto& hkl : hklList ) {
          double f = hkl.fsquared;
          if ( !(f>0.0) ) {
            if (f<0)
              NCRYSTAL_THROW(CalcError,"Inconsistent data implies"
                             " negative |F|^2*multiplicity.");
            continue;
          }
          updatePlaneList( data, hkl );
        }
        sortData( res );
        return res;
      }

      void checkPlaneData( const PIMergedData::plane_t& p )
      {
        if ( !std::isfinite(p.dsp) || !(p.dsp>0.0) ) {
          NCRYSTAL_THROW2(BadInput,"d-spacing not finite and positive"
                          " in list of (dsp,fsq*mult) values");
        }
        if ( !std::isfinite(p.fsqmult) || !(p.fsqmult>0.0) ) {
          NCRYSTAL_THROW2(BadInput,"fm not finite and positive"
                          " in list of (dsp,fsq*mult) values");
        }
      }

      void checkPlaneData( const PIData::plane_t& p )
      {
        if ( !std::isfinite(p.dsp) || !(p.dsp>0.0) ) {
          NCRYSTAL_THROW2(BadInput,"d-spacing not finite and positive"
                          " in list of (dsp,fsq,mult) values");
        }
        if ( !std::isfinite(p.fsq) || !(p.fsq>0.0) ) {
          NCRYSTAL_THROW2(BadInput,"fsq not finite and positive"
                          " in list of (dsp,fsq,mult) values");
        }
        if ( !std::isfinite(p.mult) || !(p.mult>=1.0) ) {
          NCRYSTAL_THROW2(BadInput,"fsq not finite and >=1"
                          " in list of (dsp,fsq,mult) values");
        }
        if ( std::trunc(p.mult) != p.mult ) {
          NCRYSTAL_THROW2(BadInput,"mult not an integral value "
                          "in list of (dsp,fsq,mult) values");
        }
      }

      template<class TData>
      void checkDataImpl( const TData& data )
      {
        checkData( data.cell );
        for ( auto& e : data.planes )
          checkPlaneData( e );
      }

      bool cmpMergedPlanes(const PIMergedData::plane_t& a,
                           const PIMergedData::plane_t& b)
      {
        nc_assert( std::isfinite(a.dsp) );
        nc_assert( std::isfinite(b.dsp) );
        nc_assert( std::isfinite(a.fsqmult) );
        nc_assert( std::isfinite(b.fsqmult) );
        //larger d-spacing first, then larger fsqmult first
        if ( a.dsp != b.dsp )
          return a.dsp > b.dsp;//larger d-spacing first
        return a.fsqmult > b.fsqmult;
      }

      bool cmpPlanes(const PIData::plane_t& a, const PIData::plane_t& b)
      {
        nc_assert( std::isfinite(a.dsp) );
        nc_assert( std::isfinite(b.dsp) );
        nc_assert( std::isfinite(a.fsq) );
        nc_assert( std::isfinite(b.fsq) );
        nc_assert( std::isfinite(a.mult) );
        nc_assert( std::isfinite(b.mult) );
        if ( a.dsp != b.dsp )
          return a.dsp > b.dsp;//larger d-spacing first
        if ( a.fsq != b.fsq )
          return a.fsq > b.fsq;//then larger fsq first
        return a.mult > b.mult;//then larger mult first
      }
    }
  }
}

NCPB::PIMergedData NCPB::prepareMergedData( const Info& info )
{
  checkInfoOK(info);
  auto& si = info.getStructureInfo();
  return prepareDataImpl<PIMergedData>( info.hklList(), prepareCellData( si ) );
}

NCPB::PIData NCPB::prepareData( const Info& info )
{
  checkInfoOK(info);
  return prepareDataImpl<PIData>( info.hklList(),
                                  prepareCellData( info.getStructureInfo() ) );
}

NCPB::CellData NCPB::prepareCellData( const StructureInfo& si )
{
  if ( !cellVolNatomsOK(si.volume,si.n_atoms) )
    NCRYSTAL_THROW(BadInput,"Passed structure info object has"
                   " invalid volume or n_atoms fields.");
  nc_assert_always(si.n_atoms>0);
  nc_assert_always(si.volume>0);
  CellData res;
  res.volume = si.volume;
  res.n_atoms = si.n_atoms;
  checkData(res);
  return res;
}

NCPB::CellData NCPB::prepareCellData( double volume,
                                      unsigned n_atoms )
{
  if ( !cellVolNatomsOK(volume,n_atoms) )
    NCRYSTAL_THROW(BadInput,"Passed cell volume or n_atoms data are invalid.");
  CellData res;
  res.volume = volume;
  res.n_atoms = n_atoms;
  checkData(res);
  return res;
}

void NCPB::checkData( const CellData& cd )
{
  if ( !cellVolNatomsOK(cd.volume,cd.n_atoms) )
    NCRYSTAL_THROW(BadInput,"CellData object has"
                   " invalid volume or n_atoms fields.");
}

void NCPB::checkData( const PIMergedData& data )
{
  checkDataImpl( data );
  auto & v = data.planes;
  if ( !std::is_sorted(std::begin(v), std::end(v), cmpMergedPlanes))
    NCRYSTAL_THROW2(BadInput,"Received unsorted list of (dsp,fsq*mult) values");
}

void NCPB::checkData( const PIData& data )
{
  checkDataImpl( data );
  auto & v = data.planes;
  if ( !std::is_sorted(std::begin(v), std::end(v), cmpPlanes))
    NCRYSTAL_THROW2(BadInput,"Received unsorted list of (dsp,fsq,mult) values");
}

void NCPB::sortData( PIMergedData::PlaneList& v )
{
  std::sort(v.begin(),v.end(),cmpMergedPlanes);
}

void NCPB::sortData( PIData::PlaneList& v )
{
  std::sort(v.begin(),v.end(),cmpPlanes);
}

void NCPB::sortData( PIMergedData& data )
{
  sortData(data.planes);
}

void NCPB::sortData( PIData& data )
{
  sortData(data.planes);
}

NCPB::PIMergedData NCPB::prepareMergedData( const PIData& input )
{
  //fixme: unit test
  checkData(input);
  NCPB::PIMergedData res;
  res.cell = input.cell;
  auto& tgt = res.planes;
  tgt.reserve( input.planes.size() );
  for ( auto& e : input.planes ) {
    if ( tgt.empty() || tgt.back().dsp != e.dsp ) {
      tgt.emplace_back();
      tgt.back().dsp = e.dsp;
      tgt.back().fsqmult = e.fsq * e.mult;
    } else {
      tgt.back().fsqmult += e.fsq * e.mult;
    }
  }
  checkData(res);
  return res;
}

NCPB::PIMergedData NCPB::prepareMergedData( const StructureInfo& si,
                                            PIMergedData::PlaneList&& v )
{
  return prepareMergedData( prepareCellData( si ), std::move(v) );
}

NCPB::PIMergedData NCPB::prepareMergedData( const CellData& c,
                                            PIMergedData::PlaneList&& v )
{
  sortData(v);
  PIMergedData res;
  res.cell = c;
  res.planes = std::move(v);
  checkData(res);
  return res;
}

NCPB::PIData NCPB::prepareData( const StructureInfo& si,
                                PIData::PlaneList&& v )
{
  return prepareData( prepareCellData( si ), std::move(v) );
}

NCPB::PIData NCPB::prepareData( const CellData& c,
                                PIData::PlaneList&& v )
{
  sortData(v);
  PIData res;
  res.cell = c;
  res.planes = std::move(v);
  checkData(res);
  return res;
}

