////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCrystal/NCLoadNCMAT.hh"
#include "NCParseNCMAT.hh"
#include "NCrystal/NCDefs.hh"
#include "NCFillHKL.hh"
#include "NCFile.hh"
#include "NCRotMatrix.hh"
#include "NCNeutronSCL.hh"
#include "NCLatticeUtils.hh"
#include "NCDebyeMSD.hh"

#include <vector>
#include <iostream>
#include <cstdlib>

const NCrystal::Info * NCrystal::loadNCMAT( const char * ncmat_file,
                                            double temperature,
                                            double dcutoff,
                                            double dcutoffup,
                                            bool expandhkl )
{
  const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
  if (verbose)
    std::cout<<"NCrystal::loadNCMAT called with ("<< ncmat_file
             <<", temp="<<temperature
             <<", dcutoff="<<dcutoff
             <<", dcutoffup="<<dcutoffup
             <<", expandhkl="<<expandhkl<<")"<<std::endl;

  NCMATParser parser(ncmat_file);
  nc_assert_always(parser.getAtomPerCell()>0);

  ////////////////////////
  // Create Info object //
  ////////////////////////

  Info * crystal = new Info();
  const NeutronSCL *nscl = NeutronSCL::instance();

  //1. structure info
  StructureInfo si;
  si.spacegroup = parser.getSpacegroupNum();
  parser.getLatticeParameters( si.lattice_a, si.lattice_b, si.lattice_c,
                               si.alpha, si.beta, si.gamma );
  const RotMatrix cell = getLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                        si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
  si.volume = cell.colX().cross(cell.colY()).dot(cell.colZ());
  si.n_atoms = parser.getAtomPerCell();
  crystal->setStructInfo(si);

  //2. temperature info
  if (temperature>0.0)
    crystal->setTemperature(temperature);
  if (parser.getDebyeTemp()>0.0)
    crystal->setDebyeTemperature(parser.getDebyeTemp());

  //3. atom info, cross sections and density:
  double cell_mass_amu = 0.;
  double totalFreeScaXS = 0.;
  double totalAbsXS = 0.;

  const std::map<std::string, double>& debye_map = parser.getDebyeMap();

  for(std::map<std::string, std::vector<Vector> >::const_iterator it = parser.getAtomicPosMap().begin();
      it != parser.getAtomicPosMap().end();
      ++it )
    {
      const std::string& elem_name = it->first;
      double elem_mass_amu = nscl->getAtomicMass(elem_name);
      size_t npos = it->second.size();
      //if per-element debye temps exists, they are preferable and should be
      //available for all elements. Otherwise a global debye temp should be
      //available (asserting here as a LogicError, since the parser should already
      //have checked this):
      double debye_temp = ( debye_map.empty() ? parser.getDebyeTemp() : debye_map.at(elem_name) );
      nc_assert_always( debye_temp > 0.0 );

      double msd = debyeIsotropicMSD(debye_temp, temperature, elem_mass_amu );


      //cross sections and mass:
      unsigned atnum = npos;
      totalAbsXS += atnum* nscl->getCaptureXS(elem_name);
      totalFreeScaXS += atnum* nscl->getFreeXS(elem_name);
      cell_mass_amu += elem_mass_amu * npos;

      //atominfo
      AtomInfo ai;
      ai.atomic_number = nscl->getAtomicNumber(elem_name);
      ai.number_per_unit_cell = atnum;
      ai.mean_square_displacement = msd;
      if (!debye_map.empty()) {
        nc_assert_always(debye_map.find(elem_name)!=debye_map.end());//should have been checked by parser
        ai.debye_temp = debye_map.find(elem_name)->second;
      } else {
        ai.debye_temp = 0;
      }
      ai.positions.reserve(npos);
      for (std::vector<Vector>::const_iterator itPos = it->second.begin(); itPos!=it->second.end(); ++itPos) {
        ai.positions.push_back(AtomInfo::Pos(itPos->x(),itPos->y(),itPos->z()));
      }
      crystal->addAtom(ai);
    }

  crystal->setXSectAbsorption( totalAbsXS / parser.getAtomPerCell() );
  crystal->setXSectFree( totalFreeScaXS / parser.getAtomPerCell() );
  crystal->setDensity( 1e27 * cell_mass_amu * constant_dalton2kg / si.volume );//1e27 converts kg/Aa^3 to g/cm^3

  //4. HKL info
  if(dcutoff==0) {
    //Very simple heuristics here for now to select appropriate dcutoff value
    //(specifically we needed to raise the value to 0.4Aa for expensive Y2O3/SiLu2O5
    //with ~80/65 atoms/cell):
    dcutoff = ( parser.getAtomPerCell()>40 ? 0.25 : 0.1 ) ;
    std::string cmt;
    if (dcutoff>=dcutoffup) {
      //automatically selected conflicts with value of dcutoffup.
      cmt = " (lower than usual due to value of dcutoffup)";
      dcutoff = 0.5*dcutoffup;
    }
    if (verbose)
      std::cout<<"NCrystal::NCMATFactory::automatically selected dcutoff level "<< dcutoff << " Aa"<<cmt<<std::endl;
  }

  if ( dcutoff != -1 ) {
    const double fsquare_cut = 1e-5;//NB: Hardcoded to same value as in .nxs factory
    const double merge_tolerance = 1e-6;
    fillHKL(*crystal,  dcutoff , dcutoffup, expandhkl, fsquare_cut, merge_tolerance);
  }


  ///////////
  // Done! //
  ///////////

  crystal->objectDone();
  return crystal;

}

