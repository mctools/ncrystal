////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCDefs.hh"
#include "NCString.hh"
#include "NCMath.hh"
#include "NCLatticeUtils.hh"
#include <sstream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iostream>

//TODO for NC2: lower priority version without sginfo dep could go to dist/

namespace NCrystal {
  double str2dbl_laz(const std::string& s) { return str2dbl(s,"Invalid number in .laz/.lau file"); }
  double str2int_laz(const std::string& s) { return str2int(s,"Invalid integer in .laz/.lau file"); }
}

NCrystal::LazLoader::LazLoader(std::string laz_file,double dcutlow, double dcutup, double temp)
  : m_full_path(laz_file), m_dcutlow(dcutlow), m_dcutup(dcutup), m_temp(temp)
{
  nc_assert_always( dcutlow==-1 || ( dcutlow>=0 && dcutlow < dcutup ) );
  m_cinfo = new NCrystal::Info();
  m_cinfo->ref();
}

bool NCrystal::LazLoader::setupSgInfo(unsigned spaceGroupNbr, nxs::T_SgInfo& sgInfo) {
  nc_assert_always(!nxs::SgError);
  sgInfo.MaxList = 1024;
  sgInfo.ListSeitzMx = (nxs::T_RTMx*)malloc( sgInfo.MaxList * sizeof(*sgInfo.ListSeitzMx) );
  /* no list info needed here */
  sgInfo.ListRotMxInfo = NULL;
  const nxs::T_TabSgName *tsgn = NULL;
  std::stringstream s;
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

NCrystal::LazLoader::~LazLoader()
{
  m_cinfo->unref();
}

NCrystal::Info* NCrystal::LazLoader::getCrystalInfo()
{
  return m_cinfo;
}

void NCrystal::LazLoader::read()
{
  std::ifstream infile;
  std::string line="";
  std::string line2;
  infile.open(m_full_path.c_str());
  if (!infile.good())
    NCRYSTAL_THROW2(FileNotFound,"Could not find and open input file \""<<m_full_path<<"\"");

  //copy ascii raw data into memory as a vector of strings
  std::vector<std::string> strVec;
  while (std::getline(infile, line))
    {
      if (line.empty())
        continue;
      split(strVec,line);
      if (strVec.empty())
        continue;
      if (line[0]=='#')
        m_raw_header.push_back(strVec);
      else
        m_raw_data.push_back(strVec);
    }

  infile.close();


  //set density
  double density=0.;
  if(search_parameter("density", density ))
    m_cinfo->setDensity(density);

  //set structure info
  NCrystal::StructureInfo structure_info;
  if(!search_spacegroup(structure_info.spacegroup))
    NCRYSTAL_THROW2(DataLoadError,"Can not find space group definition in the input file \""<<m_full_path<<"\"");

  if(!search_multiplicity(structure_info.n_atoms))
    NCRYSTAL_THROW2(DataLoadError,"Can not find multiplicity in the input file \""<<m_full_path<<"\"");
  double n = 0.;
  search_parameter("multiplicity", n );
  structure_info.n_atoms *= n;

  if(!search_parameter("lattice_a", structure_info.lattice_a ))
    NCRYSTAL_THROW2(DataLoadError,"The unit cell lattice_a is not defined in the input file \""<<m_full_path<<"\"");

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
    NCRYSTAL_THROW2(DataLoadError,"The unit cell volume is not defined in the input file \""<<m_full_path<<"\"");//TODO for NC2: just calculate (after completing structure info below)


  //Delayed until after sanity check below: m_cinfo->setStructInfo(structure_info);

  //simply pass along requested temperature (there is anyway no material
  //temperature info given in the lazy file):
  m_cinfo->setTemperature(m_temp);

  //set cross section info
  double sigma_abs(-1.0);
  if(search_parameter("sigma_abs", sigma_abs ))
    m_cinfo->setXSectAbsorption(sigma_abs);

  //structure factors
  unsigned d_index = 0;
  if(!search_index("column_d", d_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for d-spacing in the table is not defined in the input file \""<<m_full_path<<"\"");

  unsigned mul_index = 0;
  if(!search_index("column_j", mul_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for multiplicity in the table is not defined in the input file \""<<m_full_path<<"\"");

  unsigned f_index = 0;
  unsigned f2_index = 0;
  if(! (search_index("column_F", f_index) || search_index("column_F2", f2_index) ))
    NCRYSTAL_THROW2(DataLoadError,"The index for structure factor in the table is not defined in the input file \""<<m_full_path<<"\"");

  unsigned h_index = 0;
  if(!search_index("column_h", h_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for h in the table is not defined in the input file \""<<m_full_path<<"\"");

  unsigned k_index = 0;
  if(!search_index("column_k", k_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for k in the table is not defined in the input file \""<<m_full_path<<"\"");

  unsigned l_index = 0;
  if(!search_index("column_l", l_index) )
    NCRYSTAL_THROW2(DataLoadError,"The index for l in the table is not defined in the input file \""<<m_full_path<<"\"");

  if (m_raw_data.empty())
    NCRYSTAL_THROW2(DataLoadError,"No data lines found in input file \""<<m_full_path<<"\"");

  // double lower_d = str2dbl_laz(m_raw_data.back().at(d_index-1));
  // double upper_d = str2dbl_laz(m_raw_data.begin()->at(d_index-1));

  checkAndCompleteLattice( structure_info.spacegroup, structure_info.lattice_a, structure_info.lattice_b, structure_info.lattice_c );

  const bool enable_hkl(m_dcutlow!=-1);
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
        //TODO for NC2: Check globally against repeated lines.
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
            m_cinfo->addHKL(info);
          }
      }

    if ( m_dcutlow == 0 ) {
      //m_dcutlow=0 means auto-select, which in laz/lau files we implement as put
      //to the smallest dspacing in the file (i.e. use all lines in file).
      m_dcutlow = dlow>0 ? dlow : 0.9*m_dcutup;
      if (std::getenv("NCRYSTAL_DEBUGINFO"))
        std::cout<<"NCrystal::NCLAZFactory::automatically selected dcutoff level "<< m_dcutlow << " Aa"<<std::endl;
    }
    m_cinfo->enableHKLInfo(m_dcutlow,m_dcutup);
    if (m_cinfo->hasHKLInfo()) {
      {
        //Sanity-check .laz/.lau header, but verifying that no forbidden hkl planes were found in the list.
        nxs::T_SgInfo sgInfo;
        const char * old_SgError = nxs::SgError;
        nxs::SgError = 0;
        if (!setupSgInfo(structure_info.spacegroup,sgInfo)||nxs::SgError) {
          printf("NCrystal::loadLazCrystal WARNING: Problems setting up space-group info for sanity checking file.\n");
          if (nxs::SgError)
            printf("NCrystal::loadLazCrystal WARNING: Problem was \"%s\".\n",nxs::SgError);
          nxs::SgError = old_SgError;
        } else {
          nxs::SgError = old_SgError;
          //Do stuff
          HKLList::const_iterator it = m_cinfo->hklBegin();
          HKLList::const_iterator itE = m_cinfo->hklEnd();
          int warn(0);
          for (;it!=itE;++it) {
            int dummy(0);
            if (IsSysAbsent_hkl(&sgInfo,it->h,it->k,it->l,&dummy)!=0) {
              printf("NCrystal::loadLazCrystal WARNING: Forbidden hkl=(%i,%i,%i) for spacegroup found in file.\n",it->h,it->k,it->l);
              ++warn;
              if (warn==5) {
                printf("NCrystal::loadLazCrystal WARNING: Suppressing further warnings about forbidden hkl for this file.\n");
                break;
              }
              structure_info_is_sane=false;
            }
          }
        }
        free(sgInfo.ListSeitzMx);
      }
    }
  }
  if (structure_info_is_sane)
    m_cinfo->setStructInfo(structure_info);
  else
    printf("NCrystal::loadLazCrystal WARNING: Apparently forbidden hkl values found indicating malformed header in .laz file \"%s\"\n",m_full_path.c_str());

  m_cinfo->objectDone();
}

bool NCrystal::LazLoader::search_parameter(std::string attr, double &result )
{
  for (RawItr it = m_raw_header.begin() ; it != m_raw_header.end(); ++it)
    {
      for(StrVecItr it_str = (*it).begin() ; it_str != (*it).end(); ++it_str )
        {
          if(attr==(*it_str))
            {
              result = str2dbl_laz( (*++it_str));
              return true;
            }
        }
    }
  return false;
}

bool NCrystal::LazLoader::search_index(std::string attr, unsigned &result)
{
  for (RawItr it = m_raw_header.begin() ; it != m_raw_header.end(); ++it)
    {
      for(StrVecItr it_str = (*it).begin() ; it_str != (*it).end(); ++it_str )
        {
          if(attr==(*it_str))
            {
              result = str2int_laz( (*++it_str));
              return true;
            }
        }
    }
  return false;
}

unsigned NCrystal::LazLoader::countAtom(std::string formula)
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


bool NCrystal::LazLoader::search_spacegroup(unsigned &result)
{
  for (RawItr it = m_raw_header.begin() ; it != m_raw_header.end(); ++it)
    {
      for(StrVecItr it_str = (*it).begin() ; it_str != (*it).end(); ++it_str )
        {
          if((*it_str)=="SPCGRP")
            {
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
                printf("NCrystal::loadLazCrystal ERROR: Problems using SgInfo to decode space group symbol\n");
                printf("NCrystal::loadLazCrystal ERROR: Problem was \"%s\".\n",nxs::SgError);
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

bool NCrystal::LazLoader::search_multiplicity(unsigned &result )
{
  for (RawItr it = m_raw_header.begin() ; it != m_raw_header.end(); ++it)
    {
      nc_assert_always(!it->empty());
      for(StrVecItr it_str = (*it).begin() ; it_str != (*it).end(); ++it_str )
        {
          if((*it_str)=="multiplicity")
            {
              if (it->size()<3)
                return false;
              std::string full_string = it->at(it->size()-2);
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

