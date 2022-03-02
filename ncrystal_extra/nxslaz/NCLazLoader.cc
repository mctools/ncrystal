////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2022 NCrystal developers                                   //
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

#include "NCLazLoader.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include <sstream>
#include <iostream>

namespace NC = NCrystal;

//NB: lower priority version without sginfo dep could go to ncrystal core code? But is it worth it?

namespace NCrystal {
  double str2dbl_laz(const std::string& s) { return str2dbl(s,"Invalid number in .laz/.lau data"); }
  int str2int_laz(const std::string& s) { return str2int(s,"Invalid integer in .laz/.lau data"); }
}

bool NC::LazLoader::setupSgInfo(unsigned spaceGroupNbr, nxs::T_SgInfo& sgInfo) const
{
  nc_assert_always(!nxs::SgError);
  sgInfo.MaxList = 1024;
  sgInfo.ListSeitzMx = (nxs::T_RTMx*)malloc( sgInfo.MaxList * sizeof(*sgInfo.ListSeitzMx) );
  /* no list info needed here */
  sgInfo.ListRotMxInfo = NULL;
  const nxs::T_TabSgName *tsgn = NULL;
  std::ostringstream s;
  s<<spaceGroupNbr;
  char spaceGroup[1024];
  strncpy(spaceGroup,s.str().c_str(),1023);
  if( isdigit(spaceGroup[0]) ) {
    tsgn = nxs::FindTabSgNameEntry(spaceGroup, 'A');
    if (!tsgn||nxs::SgError)
      return false;
    strncpy(spaceGroup,tsgn->HallSymbol,1023);
  }
  if (nxs::SgError)
    return false;
  InitSgInfo( &sgInfo );
  if (nxs::SgError)
    return false;
  sgInfo.TabSgName = tsgn;
  if ( tsgn )
    sgInfo.GenOption = 1;
  ParseHallSymbol( spaceGroup, &sgInfo );
  if (nxs::SgError)
    return false;
  CompleteSgInfo( &sgInfo );
  if (nxs::SgError)
    return false;
  return true;
}

NC::LazLoader::LazLoader(const TextData& data,double dcutlow, double dcutup, Temperature temp)
  : m_inputDescription(data.dataSourceName()),
    m_dcutlow(dcutlow),
    m_dcutup(dcutup),
    m_temp(temp)
{
  nc_assert_always( dcutlow==-1 || ( dcutlow>=0 && dcutlow < dcutup ) );
  preParse(data);//copies all data - pretty obsolete but I don't feel like
                 //rewriting the Laz loader without some sort of reason...
}

void NC::LazLoader::preParse(const TextData& data)
{
  //copy ascii raw data into memory as a vector of strings
  for ( const auto& line : data )    {
    if (line.empty())
      continue;
    auto parts = split2( line );
    if ( parts.empty() )
      continue;
    if (parts.front().at(0)=='#')
      m_raw_header.push_back(std::move(parts));
    else
      m_raw_data.push_back(std::move(parts));
  }
}

NC::InfoBuilder::SinglePhaseBuilder NC::LazLoader::read()
{
  NC::InfoBuilder::SinglePhaseBuilder builder;

  //set density
  double density=0.;
  if(search_parameter("density", density ))
    builder.density = Density{ density };

  //set structure info
  StructureInfo structure_info;
  if(!search_spacegroup(structure_info.spacegroup)||structure_info.spacegroup==0)
    NCRYSTAL_THROW2(DataLoadError,"Can not find space group definition in the input data \""<<m_inputDescription<<"\"");

  if(!search_multiplicity(structure_info.n_atoms) || structure_info.n_atoms==0)
    NCRYSTAL_THROW2(DataLoadError,"Can not find multiplicity in the input data \""<<m_inputDescription<<"\"");
  double n = 0.;
  search_parameter("multiplicity", n );
  structure_info.n_atoms *= n;

  if(!search_parameter("lattice_a", structure_info.lattice_a ))
    NCRYSTAL_THROW2(DataLoadError,"The unit cell lattice_a is not defined in the input data \""<<m_inputDescription<<"\"");

  if(!search_parameter("lattice_b", structure_info.lattice_b ))
    structure_info.lattice_b = 0;

  if(!search_parameter("lattice_c", structure_info.lattice_c ))
    structure_info.lattice_c = 0;

  if(!search_parameter("lattice_aa", structure_info.alpha ))
    structure_info.alpha = 90;

  if(!search_parameter("lattice_bb", structure_info.beta ))
    structure_info.beta = 90;

  if(!search_parameter("lattice_cc", structure_info.gamma ))
    structure_info.gamma = 90;

  if(!search_parameter("Vc", structure_info.volume ))
    structure_info.volume = 0.0;//Ok: InfoBuilder will simply calculate.

  //Need to determine composition. First we look for a special line like:
  //
  //# formula Al2O3

  std::string composition_chemform_str;
  for ( auto& l : m_raw_header ) {
    auto it = l.begin();
    auto itE = l.end();
    while ( it!=itE && startswith(*it,"#") )
      ++it;
    if ( it!=itE && std::next(it)!=itE && *it == "formula" ) {
      composition_chemform_str = *std::next(it);
      break;
    }
  }
  if ( composition_chemform_str.empty() ) {
    //Try to fix monoatomic files by looking for ATOM lines.
    SmallVector<std::string,4> atom_strs;
    for ( auto& l : m_raw_header ) {
      auto it = l.begin();
      auto itE = l.end();
      while ( it!=itE && startswith(*it,"#") )
        ++it;
      if ( it!=itE && std::next(it)!=itE && *it == "ATOM" ) {
        if ( atom_strs.empty() || atom_strs.front() != *std::next(it) )
          atom_strs.push_back(*std::next(it));
      }
    }
    if ( atom_strs.size() == 1 && ( atom_strs.front().size() == 1 || atom_strs.front().size() == 2 ) ) {
      std::string s = atom_strs.front();
      //Fix casing, e.g. AL -> Al, al->Al:
      if ( s.at(0)>='a' && s.at(0)<='z' )
        s.at(0) += ('A'-'a');
      if ( s.size()==2 && s.at(1)>='A' && s.at(1)<='Z' )
        s.at(1) -= ('A'-'a');
      composition_chemform_str = s;
    }
  }
  if ( composition_chemform_str.empty() ) {
    NCRYSTAL_THROW(BadInput,"Could not determine atomic composition from provided .laz/.lau data. Try to"
                   " add a line in the file's header specifying the formula. For instance (using Al2O3 as an example): # formula Al2O3");
  }
  builder.composition = InfoBuilder::buildCompositionFromChemForm( composition_chemform_str );

  //Delayed until after sanity check below: builder.hklPlanes = ...;

  //simply pass along requested temperature (there is anyway no material
  //temperature info given in the lazy data):
  builder.temperature = m_temp;

  //structure factors
  unsigned d_index = 0;
  if(!search_index("column_d", d_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for d-spacing in the table is not defined in the input data \""<<m_inputDescription<<"\"");

  unsigned mul_index = 0;
  if(!search_index("column_j", mul_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for multiplicity in the table is not defined in the input data \""<<m_inputDescription<<"\"");

  unsigned f_index = 0;
  unsigned f2_index = 0;
  if(! (search_index("column_F", f_index) || search_index("column_F2", f2_index) ))
    NCRYSTAL_THROW2(DataLoadError,"The index for structure factor in the table is not defined in the input data \""<<m_inputDescription<<"\"");

  unsigned h_index = 0;
  if(!search_index("column_h", h_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for h in the table is not defined in the input data \""<<m_inputDescription<<"\"");

  unsigned k_index = 0;
  if(!search_index("column_k", k_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for k in the table is not defined in the input data \""<<m_inputDescription<<"\"");

  unsigned l_index = 0;
  if(!search_index("column_l", l_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for l in the table is not defined in the input data \""<<m_inputDescription<<"\"");

  if (m_raw_data.empty())
    NCRYSTAL_THROW2(DataLoadError,"No data lines found in input data \""<<m_inputDescription<<"\"");

  // double lower_d = str2dbl_laz(m_raw_data.back().at(d_index-1));
  // double upper_d = str2dbl_laz(m_raw_data.begin()->at(d_index-1));

  checkAndCompleteLattice( structure_info.spacegroup, structure_info.lattice_a, structure_info.lattice_b, structure_info.lattice_c );

  const bool enable_hkl(m_dcutlow!=-1);
  HKLList hklList;
  bool structure_info_is_sane(true);
  if (enable_hkl) {
    double dlow(0.0);
    unsigned cache_mult = 0;
    double cache_d =-1.;
    double cache_f =-1.;

    for (RawItr it = m_raw_data.begin() ; it != m_raw_data.end(); ++it)
      {
        double dspacing(str2dbl_laz((*it).at(d_index-1)));
        if (dspacing<m_dcutlow||dspacing>m_dcutup)
          continue;
        if (!dlow||dspacing<dlow) dlow = dspacing;
        HKLInfo info;
        info.h = str2int_laz((*it).at(h_index-1));
        info.k = str2int_laz((*it).at(k_index-1));
        info.l = str2int_laz((*it).at(l_index-1));
        info.dspacing = dspacing;
        info.multiplicity = str2int_laz((*it).at(mul_index-1));
        if(f_index) {
          info.fsquared = pow(str2dbl_laz((*it).at(f_index-1)),2) ;
        } else {
          info.fsquared = str2dbl_laz((*it).at(f2_index-1)) * 0.01;//0.01 is fm^2/barn
        }

        if(info.fsquared < 1e-20)//NB: Hardcoded to lower value as in
                                 //.nxs/.ncmat factory, since it was actually
                                 //specified directly in a file rather than
                                 //calculated.
          continue;

        //TODO: .lau files should get multiplicity by counting number of lines
        //(but warn if does not correspond to number of encountered lines!). This
        //is inspired by mcstas's C_graphite.lau which has incorrect
        //multiplicity (1) specified for 0,0,2/0,0,-2 but mcstas's
        //single_crystal.comp simply ignores the multiplicity parameter.
        //
        //Probably the one true way to fix this loading would be to use our
        //EqRefl calculator to check all lines, group them properly, and check
        //if all multiplicities are correct (and allow all group members to
        //specificy either 1 or the actual group multiplicity).
        bool repeated_line(false);
        if (f2_index&&ncabs(cache_d-info.dspacing)<1.0e-4&&ncabs(cache_f-info.fsquared)<1.0e-4&&cache_mult==info.multiplicity)
          repeated_line = true;//remove repeated planes in a lau file.
        //TODO: Check globally against repeated lines.
#if 0
        //debug what to do with repeated lines in .lau
        repeated_line=false;
#endif
        if(!repeated_line)
          {
            cache_d = info.dspacing;
            cache_f = info.fsquared;
            cache_mult = info.multiplicity;
#if 0
            //debug BeO
            if (f2_index)
              info.fsquared *= 100.0*6;
            else
              info.fsquared *= 15.0;
#endif
            hklList.push_back( std::move(info) );
          }
      }

    if ( m_dcutlow == 0 ) {
      //m_dcutlow=0 means auto-select, which in laz/lau files we implement as put
      //to the smallest dspacing in the file (i.e. use all lines in file).
      m_dcutlow = dlow>0 ? dlow : 0.9*m_dcutup;
      if (std::getenv("NCRYSTAL_DEBUGINFO"))
        std::cout<<"NCrystal::NCLAZFactory::automatically selected dcutoff level "<< m_dcutlow << " Aa"<<std::endl;
    }

    {
      //Sanity-check .laz/.lau header, but verifying that no forbidden hkl planes were found in the list.
      nxs::T_SgInfo sgInfo;
      const char * old_SgError = nxs::SgError;
      nxs::SgError = 0;
      if (!setupSgInfo(structure_info.spacegroup,sgInfo)||nxs::SgError) {
        printf("NCrystal::loadLazCrystal WARNING: Problems setting up space-group info for sanity checking data.\n");
        if (nxs::SgError)
          printf("NCrystal::loadLazCrystal WARNING: Problem was \"%s\".\n",nxs::SgError);
        nxs::SgError = old_SgError;
      } else {
        nxs::SgError = old_SgError;
        //Do stuff
        int warn(0);
        for ( auto& hkl : hklList ) {
          int dummy(0);
          if (IsSysAbsent_hkl(&sgInfo,hkl.h,hkl.k,hkl.l,&dummy)!=0) {
            printf("NCrystal::loadLazCrystal WARNING: Forbidden hkl=(%i,%i,%i) for spacegroup found in data.\n",hkl.h,hkl.k,hkl.l);
            ++warn;
            if (warn==5) {
              printf("NCrystal::loadLazCrystal WARNING: Suppressing further warnings about forbidden hkl for this data.\n");
              break;
            }
            structure_info_is_sane=false;
          }
        }
      }
      free(sgInfo.ListSeitzMx);
    }
  }

  if (structure_info_is_sane) {
    builder.unitcell.emplace();
    builder.unitcell.value().structinfo = std::move(structure_info);
  } else {
    std::cout<< "NCrystal::loadLazCrystal WARNING: Apparently forbidden hkl values"
      " found indicating malformed header in .laz data \""<<m_inputDescription<<'"'<<std::endl;
  }

  if (enable_hkl) {
    builder.hklPlanes.emplace();//Set up uninitialised HKLPlanes struct
    builder.hklPlanes.value().dspacingRange = { m_dcutlow,m_dcutup };
    builder.hklPlanes.value().source = std::move(hklList);
  }

  return builder;
}

bool NC::LazLoader::search_parameter(const std::string& attr, double &result ) const
{
  for ( const auto& hdrline : m_raw_header) {
    const auto n = hdrline.size();
    for ( auto i : ncrange( std::min<decltype(n)>(2,n) ) ) {
      if ( hdrline.at(i) == attr ) {
        if ( i+1 >= n )
          NCRYSTAL_THROW2(DataLoadError,"Missing value after keyword: "<<attr);
        result = str2dbl_laz( hdrline.at(i+1) );
        return true;
      }
    }
  }
  return false;
}

bool NC::LazLoader::search_index(const std::string& attr, unsigned &result) const
{
  for (RawItr it = m_raw_header.begin() ; it != m_raw_header.end(); ++it)
    {
      for(StrVecItr it_str = (*it).begin() ; it_str != (*it).end(); ++it_str )
        {
          if(attr==(*it_str))
            {
              if ( std::next(it_str) == it->end() )
                NCRYSTAL_THROW2(DataLoadError,"Missing value after keyword: "<<attr);
              result = str2int_laz( (*++it_str));
              return true;
            }
        }
    }
  return false;
}

unsigned NC::LazLoader::countAtom(const std::string& formula) const
{
  const char * sp = formula.c_str();
  unsigned nCap =0;
  unsigned nN = 0;
  uint64_t sum = 0;
  for(unsigned i=0;i<formula.size();i++)
    {
      if(sp[i]>='A' && sp[i]<='Z')
        nCap++;
    }


  for(unsigned i=0;i<formula.length();i++)
    {
      if(sp[i]>='0' && sp[i]<='9')
        {
          std::string str="";
          while( i <= formula.size() && sp[i]>='0' && sp[i]<='9')
            {
              str += sp[i++];
            }
          sum += str2int_laz(str);
          nN++;

        }
    }

  return sum+(nCap-nN);
}


bool NC::LazLoader::search_spacegroup(unsigned &result) const
{
  for (RawItr it = m_raw_header.begin() ; it != m_raw_header.end(); ++it)
    {
      for(StrVecItr it_str = (*it).begin() ; it_str != (*it).end(); ++it_str )
        {
          if((*it_str)=="SPCGRP")
            {

              if ( std::next(it_str) == it->end() )
                NCRYSTAL_THROW2(DataLoadError,"Missing value after keyword: SPCGRP");
              std::string sg_name = *(it_str+1);
              unsigned pos = 2;
              while ( (it_str+pos) != (*it).end() && (*(it_str+pos)).size() <=4 )
                {
                  sg_name += "_";
                  sg_name += *(it_str+pos);
                  pos++;
                }
              const char * old_SgError = nxs::SgError;
              nxs::SgError = 0;
              const nxs::T_TabSgName * tsgn = nxs::FindTabSgNameEntry(sg_name.c_str(), 'A');
              if (nxs::SgError) {
                NCRYSTAL_THROW2(BadInput,"loadLazCrystal ERROR: Problems using SgInfo to decode"
                                " space group symbol. Problem was: \""<<nxs::SgError<<"\".");
                nxs::SgError = old_SgError;
                return false;
              }
              nxs::SgError = old_SgError;
              result= static_cast<unsigned> (tsgn->SgNumber);
              return true;
            }
        }
    }
  return false;
}

bool NC::LazLoader::search_multiplicity(unsigned &result ) const
{
  //Look for the "SiO2" in a line like:
  //# multiplicity 3    in [SiO2/unit cell]
  //and calculate the number of atoms (1+2=3 in case of SiO2).
  //If it says:
  //# multiplicity 240  in [atoms/unit cell]
  //Then we should return 1.
  result = 0;
  for (RawItr it = m_raw_header.begin() ; it != m_raw_header.end(); ++it)
    {
      nc_assert_always(!it->empty());
      for(StrVecItr it_str = (*it).begin() ; it_str != (*it).end(); ++it_str )
        {
          if((*it_str)=="multiplicity")
            {
              auto it_fs = it_str;
              for ( auto jj : ncrange(3) ) {
                (void)jj;
                if ( ++it_fs == it->end() )
                  return false;
              }
              const std::string& full_string = *it_fs;
              std::string mul_str;
              for(unsigned i=1;i<full_string.size();++i)
                {
                  if (full_string[i]=='/')
                    break;
                  mul_str += full_string[i];
                }
              if(mul_str=="atoms")
                result=1;
              else
                result = countAtom(mul_str);

              return true;
            }
        }
    }
  return false;
}
