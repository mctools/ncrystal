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

#include "NCParseNCMAT.hh"

#include "NCString.hh"
#include "NCFile.hh"
#include "NCrystal/NCException.hh"

NCrystal::NCMATParser::NCMATParser(const char* fn)
  : m_debye_temp(-1.), m_sg(0)
{
  nc_assert(NCrystal::file_exists(fn));
  std::ifstream infile;
  infile.open(fn);
  if (!infile.good())
    NCRYSTAL_THROW2(DataLoadError,"Problems opening file: "<<fn);

  std::string tmpstr;
  getline(infile, tmpstr);
  if (tmpstr !="NCMAT v1")//Todo for NC2: update check here (still support v1, but don't allow new fields in v1)
    NCRYSTAL_THROW(DataLoadError,"File is in an unsupported NCMAT version" );
  m_version.swap(tmpstr);

  if(ignoreCharNTimes(infile,1,'@').eof())
    NCRYSTAL_THROW2(DataLoadError,"data file does not contain any @" << fn );

  infile.seekg(-1, infile.cur);

  while(infile >> tmpstr)
  {
    if(tmpstr=="@CELL")
    {
      ignoreCharNTimes(infile,1,'\n');

      double a=0., b=0.,c=0.,alpha=0.,beta=0.,gamma=0. ;
      bool hasLen=false; bool hasAng=false;
      std::streampos oldpos;
      while(infile >> tmpstr)
      {
        //#TODO for NC2. The code was using atof instead of the safer str2dbl. After the migration, it now barfs on input data like:
        //@CELL
        //  lengths 4.06869 4.06869 4.06869#Some comment here
        //When "4.06869#Some" is fed to str2dbl. We need to rewrite this whole class for much more sane parsing!

        if (tmpstr=="lengths")
        {
          infile >> tmpstr;
          a =  str2dbl(tmpstr.c_str());
          m_a=a;
          infile >> tmpstr;
          b =  str2dbl(tmpstr.c_str());
          m_b=b;
          infile >> tmpstr;
          c =  str2dbl(tmpstr.c_str());
          m_c=c;
          hasLen=true;
        }
        else if (tmpstr=="angles")
        {
          infile >> tmpstr;
          alpha =  str2dbl(tmpstr.c_str());
          m_alpha=alpha;
          infile >> tmpstr;
          beta =  str2dbl(tmpstr.c_str());
          m_beta=beta;
          infile >> tmpstr;
          gamma =  str2dbl(tmpstr.c_str());
          m_gamma=gamma;
          hasAng=true;
        }
        else if(tmpstr[0]=='@')
        {
          infile.seekg(-tmpstr.size(), infile.cur);
          break;
        }
      }

      if(!hasLen)
        NCRYSTAL_THROW2(DataLoadError,"NCMATParser data file does not contain unitcell lengths " << fn );

      if(!hasAng)
        NCRYSTAL_THROW2(DataLoadError,"NCMATParser data file does not contain unitcell angles " << fn );

    }
    else if (tmpstr=="@ATOMPOSITIONS")
    {
      m_atomnum=0;
      ignoreCharNTimes(infile,1,'\n');
      while(infile >> tmpstr)
      {
        if(tmpstr[0]=='@')
        {
          infile.seekg(-tmpstr.size(), infile.cur);
          break;
        }
        else if (tmpstr[0]!='#')
        {
          std::string ele = tmpstr;
          infile >> tmpstr;
          double px =  str2dbl(tmpstr.c_str());
          infile >> tmpstr;
          double py =  str2dbl(tmpstr.c_str());
          infile >> tmpstr;
          double pz =  str2dbl(tmpstr.c_str());

          Vector pos(px,py,pz);
          //pos.print();
          std::map<std::string, std::vector<Vector> >::iterator it = m_atomic_pos.find(ele);
          if ( it == m_atomic_pos.end() )
          {
            std::vector<Vector> vv;
            vv.push_back(pos);
            m_atomic_pos[ele]=vv;
          }
          else
          {
            it->second.push_back(pos);
          }
          m_atomnum++;
          ignoreCharNTimes(infile,1,'\n');

        }
      }
    }
    else if (tmpstr=="@SPACEGROUP")
    {
      ignoreCharNTimes(infile,1,'\n');
      while(infile >> tmpstr)
      {
        if(tmpstr[0]=='@')
        {
          infile.seekg(-tmpstr.size(), infile.cur);
          break;
        }
        else if (tmpstr[0]!='#')
        {
          m_sg = str2int(tmpstr.c_str());
        }
      }
    }
    else if (tmpstr=="@DEBYETEMPERATURE")
    {
      ignoreCharNTimes(infile,1,'\n');
      while(infile >> tmpstr)
      {
        if(tmpstr[0]=='@')
        {
          infile.seekg(-tmpstr.size(), infile.cur);
          break;
        }
        else if (tmpstr[0]=='#')  {
          ignoreCharNTimes(infile,1,'\n');
        }
        else
        {
          if (std::isdigit (tmpstr[0]) )
            m_debye_temp = str2dbl(tmpstr.c_str());
          else
          {
            std::string element = tmpstr;
            infile >> tmpstr ;
            m_debye_map[element] =  str2dbl(tmpstr.c_str());
          }
        }

      }
    }
  }
  infile.close();


  if(m_atomic_pos.size()==0)
    NCRYSTAL_THROW2(DataLoadError,"failed to load atomic positions info from file " << fn);

  if(!m_sg)
    NCRYSTAL_THROW2(DataLoadError,"unit cell space group is not provided in the file " << fn );

  //TODO for NC2: Debye temperature will become optional when other inelastic models become available
  if(m_debye_temp==-1. && m_debye_map.empty())
    NCRYSTAL_THROW2(DataLoadError,"Debye temperature is not provided in the file " << fn );

  if (!m_debye_map.empty())
    {
      std::map<std::string, std::vector<Vector> >::const_iterator itAtom, itAtomE(m_atomic_pos.end());
    for (itAtom=m_atomic_pos.begin();itAtom!=itAtomE;++itAtom) {
      if (m_debye_map.find(itAtom->first)==m_debye_map.end())
        NCRYSTAL_THROW2(DataLoadError,"Some per-elment Debye temperatures are specified but one is missing for \""<<itAtom->first<<"\" in the file: "<<fn);
    }
    if (m_debye_map.size()!=m_atomic_pos.size())
      NCRYSTAL_THROW2(DataLoadError,"Per-elment Debye temperatures for elements not in the lattice are specified in the file: "<<fn);
  }


}

NCrystal::NCMATParser::~NCMATParser()
{
}



std::map<std::string, unsigned> NCrystal::NCMATParser::getAtomMap()
{
  std::map<std::string, unsigned> an;
  std::map<std::string, std::vector<Vector> >::const_iterator it = m_atomic_pos.begin();
  for( ; it!=m_atomic_pos.end() ; ++it )
    an[it->first] = it->second.size();
  return an;
}
