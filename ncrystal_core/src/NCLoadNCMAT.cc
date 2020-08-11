////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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
#include "NCrystal/NCParseNCMAT.hh"
#include "NCrystal/NCNCMATData.hh"
#include "NCrystal/NCDefs.hh"
#include "NCFillHKL.hh"
#include "NCrystal/NCFile.hh"
#include "NCRotMatrix.hh"
#include "NCNeutronSCL.hh"
#include "NCLatticeUtils.hh"
#include "NCDebyeMSD.hh"
#include "NCSABUtils.hh"
#include "NCScatKnlData.hh"
#include "NCVDOSEval.hh"

#include <algorithm>
#include <iostream>
#include <cstdlib>

namespace NC = NCrystal;

namespace NCrystal {

  class DI_ScatKnlImpl final : public DI_ScatKnlDirect {
  public:
    virtual ~DI_ScatKnlImpl(){}

    DI_ScatKnlImpl( double fraction,
                    const std::string& elementName,
                    VectD&& egrid,
                    ScatKnlData&& data )
      : DI_ScatKnlDirect(fraction,elementName,data.temperature),
        m_inputdata(std::make_unique<ScatKnlData>(std::move(data)))
    {
      if (!egrid.empty())
        m_egrid = std::make_shared<const VectD>(std::move(egrid));
    }

    std::shared_ptr<const VectD> energyGrid() const final {return m_egrid;}


  protected:
    virtual std::shared_ptr<const SABData> buildSAB() const override final
    {
      //Someone actually requested the SAB data! Validate and (if needed) expand
      //to complete asymmetric description (thus we save a bit of memory for the
      //use-cases that don't actually need to access the SAB - like when running
      //in inelas=0 mode).
      //
      //NB: Invocation of this method is protected by object-specific mutex lock.
      nc_assert_always(!!m_inputdata);
      SABData data = SABUtils::transformKernelToStdFormat(std::move(*m_inputdata));
      return std::make_shared<const SABData>(std::move(data));
    }
  private:
    mutable std::unique_ptr<ScatKnlData> m_inputdata;
    std::shared_ptr<const VectD> m_egrid;
  };

  class DI_VDOSImpl final : public DI_VDOS {
  public:
    virtual ~DI_VDOSImpl() = default;
    DI_VDOSImpl( double fraction,
                 const std::string& elementName,
                 double temperature,
                 VectD&& egrid,
                 VDOSData&& data,
                 VectD&& orig_vdos_egrid,
                 VectD&& orig_vdos_density)
      : DI_VDOS(fraction,elementName,temperature),
        m_vdosdata(std::move(data)),
        m_vdosOrigEgrid(orig_vdos_egrid),
        m_vdosOrigDensity(orig_vdos_density)
    {
      if (!egrid.empty())
        m_egrid = std::make_shared<const VectD>(std::move(egrid));
    }

    std::shared_ptr<const VectD> energyGrid() const final {return m_egrid;}
    const VDOSData& vdosData() const final { return m_vdosdata; }

    const VectD& vdosOrigEgrid() const final { return m_vdosOrigEgrid; }
    const VectD& vdosOrigDensity() const final { return m_vdosOrigDensity; }

  private:
    VDOSData m_vdosdata;
    std::shared_ptr<const VectD> m_egrid;
    VectD m_vdosOrigEgrid;
    VectD m_vdosOrigDensity;
  };

}

const NC::Info * NC::loadNCMAT( const char * ncmat_file,
                                double temperature,
                                double dcutoff,
                                double dcutoffup,
                                bool expandhkl )
{
  nc_assert_always(ncmat_file);
  return loadNCMAT(std::string(ncmat_file),temperature,dcutoff,dcutoffup,expandhkl);
}

const NC::Info * NC::loadNCMAT( const std::string& ncmat_file,
                                double temperature,
                                double dcutoff,
                                double dcutoffup,
                                bool expandhkl )
{
  UniquePtr<TextInputStream> inputstream;
  createTextInputStream( ncmat_file, inputstream );
  NCMATParser parser(inputstream);
  NCMATData data;
  parser.getData(data);
  return loadNCMAT(std::move(data),temperature,dcutoff,dcutoffup,expandhkl);
}

const NC::Info * NC::loadNCMAT( NCMATData&& data,
                                double temperature,
                                double dcutoff,
                                double dcutoffup,
                                bool expandhkl )
{
  const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
  if (verbose)
    std::cout<<"NCrystal::loadNCMAT called with ("
             << data.sourceFullDescr
             <<", temp="<<temperature
             <<", dcutoff="<<dcutoff
             <<", dcutoffup="<<dcutoffup
             <<", expandhkl="<<expandhkl<<")"<<std::endl;

  //(re)validate data here, to be defensive:
  data.validate();

  ////////////////////////
  // Handle temperature //
  ////////////////////////

  //Check for any temperature indicated in the input data itself (only possible
  //in ScatKnl sections of dyninfos):
  double input_temperature = -1.0;
  if (data.hasDynInfo()) {
    for (auto& e : data.dyninfos) {
      if ( e.dyninfo_type != NCMATData::DynInfo::ScatKnl )
        continue;
      double dit = e.fields.at("temperature").at(0);
      if ( input_temperature == -1.0 ) {
        input_temperature = dit;
      } else {
        nc_assert_always( floateq(input_temperature, dit ) );
      }
    }
  }
  if ( temperature == -1.0 ) {
    temperature = ( input_temperature==-1.0 ? 293.15 : input_temperature );
  } else {
    if ( input_temperature != -1.0 && !floateq(input_temperature, temperature) )
      NCRYSTAL_THROW2(BadInput,data.sourceFullDescr <<" specified temperature ("<<temperature<<"K)"
                      " is incompatible with temperature ("<<input_temperature<<"K) at which input data is valid.");
  }
  nc_assert_always( temperature > 0.0 && ( input_temperature==-1.0 || floateq(input_temperature,temperature) ) );


  ////////////////////////
  // Create Info object //
  ////////////////////////

  Info * info = new Info();
  const NeutronSCL *nscl = NeutronSCL::instance();

  //Cache all the hasXXX results, before we start messing them up by
  //std::move'ing stuff out of the data object:
  const bool data_hasDynInfo = data.hasDynInfo();
  const bool data_hasUnitCell = data.hasUnitCell();
  const bool data_hasCell = data.hasCell();
  const bool data_hasAtomPos = data.hasAtomPos();
  const bool data_hasSpaceGroup = data.hasSpaceGroup();
  const bool data_hasDebyeTemperature = data.hasDebyeTemperature();
  const bool data_hasDensity = data.hasDensity();

  //Cache Debye temperatures (if any):
  const double data_debyetemp_global = data.debyetemp_global;
  std::map<std::string, double> perelemdebye_map;
  data.fillPerElementDebyeTempMap(perelemdebye_map);
  auto elementName2DebyeTemp = [data_debyetemp_global,&perelemdebye_map](std::string en)
                               {
                                 if (data_debyetemp_global)
                                   return data_debyetemp_global;
                                 auto it = perelemdebye_map.find(en);
                                 nc_assert_always(it!=perelemdebye_map.end());
                                 return it->second;
                               };


  //==> Dynamics (deal with them first as they might provide MSD info)
  std::map<std::string,double> elem2msd;
  std::map<std::string,double> elem2frac;
  std::vector<std::unique_ptr<DynamicInfo>> dyninfolist;
  auto getEgrid = [](NCMATData::DynInfo::FieldMapT& fields)
                  {
                    VectD egrid;
                    if (fields.count("egrid"))
                      egrid = std::move(fields.at("egrid"));
                    if ( egrid.size() == 1 )
                      egrid = { 0.0, egrid.front(), 0.0};//convert egrid={<emax>} to egrid={<emin>,<emax>,<npts> format
                    return egrid;
                  };

  if (data_hasDynInfo) {
    for (auto& e : data.dyninfos) {
      //nscl->getAtomicNumber(e.element_name) throws BadInput in case element
      //name is invalid, so call just to trigger that - for now:
      nscl->getAtomicNumber(e.element_name);
      nc_assert_always(e.fraction>0.0&&e.fraction<=1.0);
      std::unique_ptr<DynamicInfo> di;
      elem2frac[e.element_name] = e.fraction;
      switch (e.dyninfo_type) {
      case NCMATData::DynInfo::Sterile:
        di = std::make_unique<DI_Sterile>(e.fraction, e.element_name, temperature);
        break;
      case NCMATData::DynInfo::FreeGas:
        di = std::make_unique<DI_FreeGas>(e.fraction, e.element_name, temperature);
        break;
      case NCMATData::DynInfo::VDOSDebye:
        nc_assert_always(data_hasDebyeTemperature);
        di = std::make_unique<DI_VDOSDebye>(e.fraction, e.element_name, temperature,
                                            elementName2DebyeTemp(e.element_name));
        break;
      case NCMATData::DynInfo::VDOS:
        {
          VectD egrid = getEgrid(e.fields);
          auto vdos_egrid_orig = std::move(e.fields.at("vdos_egrid"));
          auto vdos_density_orig = std::move(e.fields.at("vdos_density"));

          nc_assert_always( vdos_egrid_orig.size()==2 || vdos_egrid_orig.size()==vdos_density_orig.size());
          nc_assert_always( vdos_density_orig.size() >= 5 );

          VectD vdos_egrid_reg, vdos_density_reg;
          std::tie(vdos_egrid_reg, vdos_density_reg) = regulariseVDOSGrid( vdos_egrid_orig, vdos_density_orig);

          nc_assert_always(vdos_egrid_reg.size()==2);
          PairDD vdos_egrid_pair(vdos_egrid_reg.front(),vdos_egrid_reg.back());
          SigmaBound boundXS = SigmaBound{NeutronSCL::instance()->getBoundXS(e.element_name)};
          double elementMassAMU = NeutronSCL::instance()->getAtomicMass(e.element_name);
          di = std::make_unique<DI_VDOSImpl>( e.fraction, e.element_name, temperature,
                                              std::move(egrid),
                                              VDOSData(vdos_egrid_pair,std::move(vdos_density_reg),
                                                       temperature,boundXS,elementMassAMU),
                                              std::move(vdos_egrid_orig),
                                              std::move(vdos_density_orig) );
        }
        break;
      case NCMATData::DynInfo::ScatKnl:
        {
          nc_assert( floateq( temperature, e.fields.at("temperature").at(0) ) );
          //Prepare scatter kernel object:
          ScatKnlData knldata;
          knldata.temperature = temperature;
          knldata.boundXS = SigmaBound{NeutronSCL::instance()->getBoundXS(e.element_name)};
          knldata.elementMassAMU = NeutronSCL::instance()->getAtomicMass(e.element_name);
          //Move acquire expensive fields:

          if (e.fields.count("sab")) {
            knldata.alphaGrid = std::move(e.fields.at("alphagrid"));
            knldata.betaGrid = std::move(e.fields.at("betagrid"));
            knldata.sab = std::move(e.fields.at("sab"));
            knldata.knltype = ScatKnlData::KnlType::SAB;
          } else if (e.fields.count("sab_scaled")) {
            knldata.alphaGrid = std::move(e.fields.at("alphagrid"));
            knldata.betaGrid = std::move(e.fields.at("betagrid"));
            knldata.sab = std::move(e.fields.at("sab_scaled"));
            nc_assert_always(knldata.betaGrid.size()>0);
            if (knldata.betaGrid.front()==0.0)
              knldata.knltype = ScatKnlData::KnlType::SCALED_SYM_SAB;
            else
              knldata.knltype = ScatKnlData::KnlType::SCALED_SAB;
          } else if (e.fields.count("sqw")) {
            knldata.alphaGrid = std::move(e.fields.at("qgrid"));
            knldata.betaGrid = std::move(e.fields.at("omegagrid"));
            knldata.sab = std::move(e.fields.at("sqw"));
            knldata.knltype = ScatKnlData::KnlType::SQW;
          } else {
            NCRYSTAL_THROW(LogicError,"Unexpected SAB type in input data");//logic-error, since we should have caught this earlier.
          }

          //Egrid:
          VectD egrid = getEgrid(e.fields);
          di = std::make_unique<DI_ScatKnlImpl>(e.fraction, e.element_name,
                                                std::move(egrid),
                                                std::move(knldata));
        }
        break;
      default:
        nc_assert_always(false);//data.validate() should have caught this.
      };
      //Collect here for now (we might adjust the fractions before adding to NCInfo):
      dyninfolist.push_back(std::move(di));
    }
  }

  //==> Temperature:
  info->setTemperature(temperature);
  double cell_volume = 0.0;
  std::size_t natoms_per_cell = 0;

  if ( data_hasUnitCell ) {
    //==> Cell layout (structure info)
    if ( data_hasCell ) {
      nc_assert_always(data_hasAtomPos);
      StructureInfo si;
      si.spacegroup = data_hasSpaceGroup ? data.spacegroup : 0;
      si.lattice_a = data.cell.lengths[0];
      si.lattice_b = data.cell.lengths[1];
      si.lattice_c = data.cell.lengths[2];
      si.alpha = data.cell.angles[0];
      si.beta = data.cell.angles[1];
      si.gamma = data.cell.angles[2];
      const RotMatrix cell = getLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                            si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
      si.volume = cell_volume = cell.colX().cross(cell.colY()).dot(cell.colZ());
      natoms_per_cell += ( si.n_atoms = data.atompos.size() );
      info->setStructInfo(si);
    }
    //==> Get sorted list of atom positions (so adjacant entries have same element name):
    decltype(data.atompos) atompos = std::move(data.atompos);
    nc_assert_always(!atompos.empty());
    std::stable_sort(atompos.begin(),atompos.end());
    //Convert to elementname->positions map:
    std::map<std::string,std::vector<AtomInfo::Pos>> elem2pos;
    std::vector<std::pair<std::string,NCMATData::Pos>>::const_iterator itAP(atompos.begin()), itAPE(atompos.end());
    for( ; itAP != itAPE; ) {
      const std::string& elem_name = itAP->first;
      std::vector<AtomInfo::Pos> positions;
      for ( ; itAP!=itAPE && itAP->first == elem_name; ++itAP )
        positions.emplace_back(itAP->second[0],itAP->second[1],itAP->second[2]);
      positions.shrink_to_fit();
      elem2pos[elem_name] = std::move(positions);
    }

    //==> Calculate element fractions (sanity check consistency with dyninfo
    //fractions, which should have already been validated by data.validate()):
    for (auto& ep : elem2pos) {
      auto& elem_name = ep.first;
      double fr_calc = double(ep.second.size()) / atompos.size();
      //Whether or not the number already exists, we update with these higher precision fractions (e.g. 0.666667
      //specified in file's dyninfo section might be consistent with 2./3. but it is not as precise):
      elem2frac[elem_name] = fr_calc;
    }

    if ( dyninfolist.empty() && data_hasDebyeTemperature ) {
      //Special case: dyninfo not specified, but we know elemental fractions
      //from atom positions, and we have debye temperatures available. We must
      //add DI_VDOSDebye entries for all elements in this case:
      for ( auto& ef: elem2frac )
        dyninfolist.emplace_back(std::make_unique<DI_VDOSDebye>(ef.second, ef.first, temperature,
                                                                elementName2DebyeTemp(ef.first)));
    }

    //Check and update fractions on dyninfo objects with numbers in elem2frac,
    //for consistency and increased precision in the scenario where numbers in
    //dyninfos have less precision than those calculated from atomic
    //compositions:
    for (auto& di: dyninfolist) {
      double precise_frac = elem2frac.at(di->elementName());
      if ( precise_frac != di->fraction() ) {
        nc_assert_always( floateq( di->fraction(), precise_frac ) );//check consistency
        di->changeFraction(precise_frac);//ensure highest precision
      }
    }

    nc_assert_always(elem2pos.size()==elem2frac.size());
    //==> Prepare to add Info::AtomInfo by figuring out MSD's for all atoms in the cell:
    nc_assert_always(data_hasDebyeTemperature);//todo: relax condition, when we have other ways of finding MSDs.
    for ( auto& ef: elem2frac ) {
      auto& elem_name = ef.first;
      if ( elem2msd.count(elem_name) ) {
        //we already calculated msd from dyninfo, just sanity check that we
        //didn't do this when we were not supposed to and proceed:
        nc_assert(!perelemdebye_map.count(ef.first));
        continue;
      }
      //Calculate MSDs from Debye temp:
      double debye_temp = ( perelemdebye_map.count(elem_name) ? perelemdebye_map.at(elem_name) : data.debyetemp_global );
      nc_assert_always(debye_temp>0.0);
      elem2msd[elem_name] = debyeIsotropicMSD( debye_temp, temperature, nscl->getAtomicMass(elem_name) );
    }

    //==> Fill Info::AtomInfo
    for ( auto& ef : elem2frac ) {
      auto& elem_name = ef.first;
      AtomInfo ai;
      ai.element_name = elem_name;
      ai.atomic_number = nscl->getAtomicNumber(elem_name);
      ai.debye_temp = ( perelemdebye_map.count(elem_name) ? perelemdebye_map.at(elem_name) : 0.0 );
      ai.mean_square_displacement = elem2msd.at(elem_name);
      ai.positions = std::move(elem2pos.at(elem_name));
      ai.number_per_unit_cell = ai.positions.size();//kind of redundant
      info->addAtom(std::move(ai));
    }

    //==> Fill global debye temp
    if (data_hasDebyeTemperature && data.debyetemp_global > 0.0)
      info->setGlobalDebyeTemperature(data.debyetemp_global);

  }//end of unit-cell specific stuff

  //==> Actually register the dyninfo objects now, since fractions have final values:
  for (auto& di: dyninfolist)
    info->addDynInfo(std::move(di));

  //==> Figure out densities and cross sections:

  double xs_abs = 0.;
  double xs_freescat = 0.;
  double avr_elem_mass_amu = 0.;
  for (auto& ef : elem2frac) {
    xs_abs            += nscl->getCaptureXS( ef.first )  * ef.second;
    xs_freescat       += nscl->getFreeXS( ef.first )     * ef.second;
    avr_elem_mass_amu += nscl->getAtomicMass( ef.first ) * ef.second;
  }

  info->setXSectAbsorption( xs_abs );
  info->setXSectFree( xs_freescat );

  double density (0.0);
  double numberdensity (0.0);

  //1e27 in next line converts kg/Aa^3 to g/cm^3:
  const double numberdensity_2_density = 1e27 * avr_elem_mass_amu * constant_dalton2kg;

  if ( cell_volume ) {
    //Calculate from unit cell
    nc_assert_always( natoms_per_cell>0 && !data_hasDensity );
    numberdensity = natoms_per_cell / cell_volume;
    density = numberdensity * numberdensity_2_density;
  } else {
    //No unit cell, user must have specified:
    nc_assert_always( !natoms_per_cell && data_hasDensity );
    if ( data.density_unit == NCMATData::KG_PER_M3 ) {
      density = data.density * 0.001;
      numberdensity = density / numberdensity_2_density;
    } else {
      nc_assert( data.density_unit == NCMATData::ATOMS_PER_AA3 );
      numberdensity = data.density;
      density = numberdensity * numberdensity_2_density;
    }
  }

  info->setDensity( density );
  info->setNumberDensity( numberdensity );

  //==> Finally populate HKL list if appropriate:
  if ( data_hasUnitCell ) {
    if(dcutoff==0) {
      //Very simple heuristics here for now to select appropriate dcutoff value
      //(specifically we needed to raise the value for expensive Y2O3/SiLu2O5
      //with ~80/65 atoms/cell):
      dcutoff = ( natoms_per_cell>40 ? 0.25 : 0.1 ) ;
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
      fillHKL(*info,  dcutoff , dcutoffup, expandhkl, fsquare_cut, merge_tolerance);
    }
  }

  ///////////
  // Done! //
  ///////////

  info->objectDone();
  return info;

}
