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

#include "NCrystal/NCInfo.hh"
#include "NCNCMatLoader.hh"
#include "NCPhononDebye.hh"
#include "NCrystal/NCException.hh"
#include "NCFile.hh"
#include "NCNeutronSCL.hh"
#include "NCReflections.hh"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>

namespace NCrystal {

  const Info * loadNCMATCrystal( const char * ncmat_file, double temperature_kelvin,
                                 double dcutoff_lower_aa, double dcutoff_upper_aa,
                                 bool expandhkl)
{
  const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
  if (verbose)
    std::cout<<"NCrystal::NCMATFactory::invoked loadNCMATCrystal("<< ncmat_file
             <<", temp="<<temperature_kelvin
             <<", dcutoff="<<dcutoff_lower_aa
             <<", dcutoffup="<<dcutoff_upper_aa
             <<", expandhkl="<<expandhkl<<")"<<std::endl;

  NCMatLoader loader(ncmat_file);
  nc_assert_always(loader.getAtomPerCell()>0);

  ////////////////////////
  // Create Info object //
  ////////////////////////

  Info * crystal = new Info();
  const NeutronSCL *nscl = NeutronSCL::instance();

  //1. structure info
  StructureInfo si;
  si.spacegroup = loader.getSpacegroupNum();
  loader.getLatticeParameters(si.lattice_a,si.lattice_b,si.lattice_c,si.alpha,si.beta,si.gamma);
  double volume = loader.getLattice().cell.colX().cross(loader.getLattice().cell.colY()).dot(loader.getLattice().cell.colZ());
  si.volume = volume;
  si.n_atoms = loader.getAtomPerCell();
  crystal->setStructInfo(si);

  //2. temperature info
  if (temperature_kelvin>0.0)
    crystal->setTemperature(temperature_kelvin);
  if (loader.getDebyeTemp()>0.0)
    crystal->setDebyeTemperature(loader.getDebyeTemp());

  //3. cross sections and atom info
  double sigma_abs=0;

  std::map<std::string, double> msd_map;
  std::vector<double> xs;

  double cellsPerGram = 1e24/volume;//1e24 is the conversion factor from Aa**3 to cm**3
  double density =0.;
  double totalFreeScaXS = 0;

  const std::map<std::string, double>& debye_map = loader.getDebyeMap();

  for(std::map<std::string, std::vector<Vector> >::const_iterator it = loader.getLattice().atomic_pos.begin();
      it != loader.getLattice().atomic_pos.end();
      ++it )
  {
    //if per-element debye temps exists, they are preferable and should be
    //available for all elements. Otherwise a global debye temp should be
    //available (asserting here as a LogicError, since the loader should already
    //have checked this):
    nc_assert_always( debye_map.empty() ? (loader.getDebyeTemp()>0.0) : (debye_map.find(it->first)!=debye_map.end()));

    //phonon inelastic (create with nphonon=1 since we just need the MSD map - we are never going to invoke inel.doit here).
    PhononDebye inel( ( debye_map.empty() ? loader.getDebyeTemp()*constant_boltzmann : debye_map.at(it->first)*constant_boltzmann ),
                      temperature_kelvin*constant_boltzmann, it->first, 1 );
    msd_map[it->first] = inel.getMSD();

    //capture, bound incoherent and coherent
    unsigned atnum = it->second.size();
    sigma_abs += atnum* nscl->getCaptureXS(it->first);

    //free sacttering xs
    totalFreeScaXS += atnum* nscl->getFreeXS(it->first);
    //density
    density += it->second.size()*nscl->getNeutronWeightedMass(it->first)*const_neutron_mass;

    //atominfo
    AtomInfo ai;
    ai.atomic_number = nscl->getAtomicNumber(it->first);
    ai.number_per_unit_cell = atnum;
    if (!debye_map.empty()) {
      nc_assert_always(debye_map.find(it->first)!=debye_map.end());//should have been checked by loader
      ai.debye_temp = debye_map.find(it->first)->second;
    } else {
      ai.debye_temp = 0;
    }
    crystal->addAtom(ai);
  }
  density*=cellsPerGram;
  sigma_abs /= loader.getAtomPerCell();
  double sigma_free = totalFreeScaXS/ loader.getAtomPerCell();
  crystal->setXSectAbsorption(sigma_abs);
  crystal->setXSectFree(sigma_free);

  //4. HKL info
  if ( dcutoff_lower_aa != -1 )
  {
    if(dcutoff_lower_aa==0) {
    //Very simple heuristics here for now, can likely be improved (specifically
    //we needed 0.4Aa for expensive Y2O3 with ~80 atoms/cell):
      if (loader.getAtomPerCell()>40)
        dcutoff_lower_aa = 0.4;
      else
        dcutoff_lower_aa = 0.15;
      std::string cmt;
      if (dcutoff_lower_aa>=dcutoff_upper_aa) {
        cmt = " (lower than usual due to value of dcutoffup)";
        dcutoff_lower_aa = 0.5*dcutoff_upper_aa;
      }
      if (verbose)
        std::cout<<"NCrystal::NCMATFactory::automatically selected dcutoff level "<< dcutoff_lower_aa << " Aa"<<cmt<<std::endl;
    }
    crystal->enableHKLInfo(dcutoff_lower_aa,dcutoff_upper_aa);
    loader.fillHKL(*crystal,  msd_map , dcutoff_lower_aa , dcutoff_upper_aa, expandhkl);
  }

  //5. density
  crystal->setDensity(density);


  ///////////
  // Done! //
  ///////////

  crystal->objectDone();
  return crystal;

}
}
///////////////////////////////////////////////////////////////////
//Finally, the code which can be used to enable our NCMAT loader //
///////////////////////////////////////////////////////////////////

#include "NCrystal/NCFactoryRegistry.hh"
#include "NCrystal/NCMatCfg.hh"

namespace NCrystal {

  class NCMATFactory : public FactoryBase {
  public:
    const char * getName() const { return "NCrystalNCMATFactory"; }

    virtual int canCreateInfo( const MatCfg& cfg ) const {
      return cfg.getDataFileExtension()=="ncmat" ? 100 : 0;
    }
    virtual const Info * createInfo( const MatCfg& cfg ) const
    {
      nc_assert_always(canCreateInfo(cfg));
      return loadNCMATCrystal( cfg.getDataFile().c_str(), cfg.get_temp(),
                               cfg.get_dcutoff(), cfg.get_dcutoffup(),
                               cfg.get_expandhkl() );
    }
  };

}


//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_ncmat_factory()
{
  if (!NCrystal::hasFactory("NCrystalNCMATFactory"))
    registerFactory(new NCrystal::NCMATFactory);
}
