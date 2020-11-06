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
#include "NCrystal/internal/NCAtomDBExtender.hh"
#include "NCrystal/internal/NCFillHKL.hh"
#include "NCrystal/NCFile.hh"
#include "NCrystal/internal/NCRotMatrix.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCDebyeMSD.hh"
#include "NCrystal/internal/NCIter.hh"
#include "NCrystal/internal/NCSABUtils.hh"
#include "NCrystal/internal/NCScatKnlData.hh"
#include "NCrystal/internal/NCVDOSEval.hh"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <atomic>

namespace NC = NCrystal;

namespace NCrystal {

  class DI_ScatKnlImpl final : public DI_ScatKnlDirect {
  public:
    virtual ~DI_ScatKnlImpl(){}

    DI_ScatKnlImpl( double fraction,
                    IndexedAtomData atom,
                    VectD&& egrid,
                    ScatKnlData&& data )
      : DI_ScatKnlDirect(fraction,std::move(atom),data.temperature),
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
                 IndexedAtomData atom,
                 double temperature,
                 VectD&& egrid,
                 VDOSData&& data,
                 VectD&& orig_vdos_egrid,
                 VectD&& orig_vdos_density)
      : DI_VDOS(fraction,std::move(atom),temperature),
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

namespace NCrystal {

  static std::atomic<bool> s_NCMATWarnOnCustomSections(!getenv("NCRYSTAL_NCMAT_NOWARNFORCUSTOM"));

}

bool NC::getNCMATWarnOnCustomSections()
{
  return s_NCMATWarnOnCustomSections;
}

void NC::setNCMATWarnOnCustomSections(bool bb)
{
  s_NCMATWarnOnCustomSections = bb;
}

const NC::Info * NC::loadNCMAT( const char * ncmat_file,
                                NC::NCMATCfgVars&& cfgvars )
{
  nc_assert_always(ncmat_file);
  return loadNCMAT( std::string(ncmat_file), std::move(cfgvars) );
}

const NC::Info * NC::loadNCMAT( const std::string& ncmat_file,
                                NC::NCMATCfgVars&& cfgvars )
{
  auto inputstream = createTextInputStream( ncmat_file );
  const bool doFinalValidation = false;
  //don't validate at end of the parseNCMATData call, since the loadNCMAT call
  //anyway validates.
  NCMATData data = parseNCMATData(std::move(inputstream),doFinalValidation);
  return loadNCMAT( std::move(data), std::move(cfgvars) );
}

const NC::Info * NC::loadNCMAT( NCMATData&& data,
                                NC::NCMATCfgVars&& cfgvars )
{
  const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
  if (verbose) {
    std::cout<<"NCrystal::loadNCMAT called with ("
             << data.sourceFullDescr
             <<", temp="<<cfgvars.temp
             <<", dcutoff="<<cfgvars.dcutoff
             <<", dcutoffup="<<cfgvars.dcutoffup
             <<", expandhkl="<<cfgvars.expandhkl
             <<", atomdb=";
    if (cfgvars.atomdb.empty()) {
      std::cout<<"<none>";
    } else {
      for (unsigned i = 0; i < cfgvars.atomdb.size(); ++i) {
        if (i>0)
          std::cout<<"@";
        std::cout<<joinstr(cfgvars.atomdb.at(i),":");
      }
    }
    std::cout<<")"<<std::endl;
  }

  /////////////////////////
  // Map isotope aliases //
  /////////////////////////

  //The specification says that D and T will be treated exactly as if the file
  //had instead contained H2 and H3 in those positions. The easiest way to
  //handle that is to perform the mapping here, before the data is parsed further.

  data.unaliasElementNames();

  //////////////
  // Validate //
  //////////////

  //Be defensive and (re)validate here. This should be done *after* the call to
  //unaliasElementNames above.
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
  if ( cfgvars.temp == -1.0 ) {
    cfgvars.temp = ( input_temperature==-1.0 ? 293.15 : input_temperature );
  } else {
    if ( input_temperature != -1.0 && !floateq(input_temperature, cfgvars.temp) )
      NCRYSTAL_THROW2(BadInput,data.sourceFullDescr <<" specified temperature ("<<cfgvars.temp<<"K)"
                      " is incompatible with temperature ("<<input_temperature<<"K) at which input data is valid.");
  }
  nc_assert_always( cfgvars.temp > 0.0 && ( input_temperature==-1.0 || floateq(input_temperature,cfgvars.temp) ) );

  /////////////////////////////////////////////////////////
  // Setup AtomDB, possibly extended via @ATOMDB section //
  /////////////////////////////////////////////////////////

  //Add AtomDB "lines", first from the NCMAT data itself and secondly from the
  //configuration variable named atomdb. Both of these sources are allowed to
  //start with an initial "line" containing the word "nodefaults". If either
  //does, access to the in-built database of isotopes and natural elements is
  //disabled (thus all atom data used in the file must be provided directly in
  //the other AtomDB "lines").
  const bool ncmat_atomdb_nodefs = (data.hasAtomDB() && data.atomDBLines.at(0).size()==1 && data.atomDBLines.at(0).at(0)=="nodefaults");
  const bool cfgvar_atomdb_nodefs = ( !cfgvars.atomdb.empty() && cfgvars.atomdb.at(0).size()==1 && cfgvars.atomdb.at(0).at(0)=="nodefaults");
  const bool allowInbuilt = !( ncmat_atomdb_nodefs || cfgvar_atomdb_nodefs );

  AtomDBExtender atomdb(allowInbuilt);//<-- the database

  if (data.hasAtomDB()) {
    for (unsigned i = (ncmat_atomdb_nodefs?1:0); i<data.atomDBLines.size(); ++i)
      atomdb.addData(data.atomDBLines.at(i),data.version);
  }
  if ( !cfgvars.atomdb.empty() ) {
    for (unsigned i = (cfgvar_atomdb_nodefs?1:0); i<cfgvars.atomdb.size(); ++i)
      atomdb.addData(cfgvars.atomdb.at(i),NCMATData::latest_version);
  }

  static_assert(AtomDBExtender::latest_version==NCMATData::latest_version,"inconsistent latest_version");//consistency check (putting in this file as it is convenient)

  ////////////////////////////////////
  // Create IndexedAtomData objects //
  ////////////////////////////////////

  //Create AtomData objects and associated indices. To do this, we simply find
  //all the "element" names seen in the @ATOMPOSITIONS, @DEBYETEMPERATURE and
  //@DYNINFO sections. We do not look in the ATOMDB section, since it can
  //contain entries not used in the other sections. While doing this, we
  //construct IndexedAtomData instances (with our local atomdb) and retain them
  //for usage below.

  std::map<std::string,IndexedAtomData> elementname_2_indexedatomdata;
  {
    std::set<std::string> allNames;
    for ( const auto& e : data.atompos )
      allNames.insert(e.first);
    for ( const auto& e : data.debyetemp_perelement )
      allNames.insert(e.first);
    for ( const auto& e : data.dyninfos )
      allNames.insert(e.element_name);
    std::vector<std::pair<AtomDataSP,std::string> > v;
    for ( auto& name : allNames ) {
      v.emplace_back(atomdb.lookupAtomData(name),name);
    }
    std::sort( v.begin(), v.end(),
               [](const std::pair<AtomDataSP,std::string>& a,
                  const std::pair<AtomDataSP,std::string>& b) {
                 return ( a.first->getUniqueID() == b.first->getUniqueID() ? a.second < b.second : *a.first < *b.first );
               } );
    nc_assert_always( (uint64_t)v.size() < (uint64_t)std::numeric_limits<unsigned>::max() );
    for ( auto&& e : enumerate(v) )
      elementname_2_indexedatomdata[ e.val.second ] = IndexedAtomData{std::move(e.val.first),{static_cast<unsigned>(e.idx)}};
  }
  std::vector<const IndexedAtomData*> index2iad;
  index2iad.resize( elementname_2_indexedatomdata.size(), nullptr );
  for( const auto& e : elementname_2_indexedatomdata ) {
    nc_assert( e.second.index.value < elementname_2_indexedatomdata.size() );
    index2iad.at(e.second.index.value) = &e.second;
  }

  NCRYSTAL_DEBUGONLY(for (auto& e: index2iad) { nc_assert_always(e!=nullptr); });

  ////////////////////////
  // Create Info object //
  ////////////////////////

  Info * info = new Info();

  //Cache all the hasXXX results, before we start messing them up by
  //std::move'ing stuff out of the data object:
  const bool data_hasDynInfo = data.hasDynInfo();
  const bool data_hasUnitCell = data.hasUnitCell();
  const bool data_hasCell = data.hasCell();
  const bool data_hasAtomPos = data.hasAtomPos();
  const bool data_hasSpaceGroup = data.hasSpaceGroup();
  const bool data_hasDebyeTemperature = data.hasDebyeTemperature();
  const bool data_hasDensity = data.hasDensity();

  //Find Debye temperatures of elements:
  const double data_debyetemp_global = data.debyetemp_global;
  std::map<AtomIndex,double> perelemdebye_map;
  for( const auto& e : data.debyetemp_perelement ) {
    const auto& iad = elementname_2_indexedatomdata[e.first];
    nc_assert(iad.atomDataSP!=nullptr);
    perelemdebye_map[iad.index] = e.second;
  }
  auto element2DebyeTemp = [data_debyetemp_global,&perelemdebye_map](const AtomIndex& idx)
  {
    if (data_debyetemp_global)
      return data_debyetemp_global;
    auto it = perelemdebye_map.find(idx);
    nc_assert(it!=perelemdebye_map.end());
    return it->second;
  };

  //==> Dynamics (deal with them first as they might one day provide MSD info)
  std::map<AtomIndex,double> elem2frac;
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
      const auto& iad = elementname_2_indexedatomdata[e.element_name];
      nc_assert_always(e.fraction>0.0&&e.fraction<=1.0);
      std::unique_ptr<DynamicInfo> di;
      elem2frac[iad.index] = e.fraction;
      switch (e.dyninfo_type) {
      case NCMATData::DynInfo::Sterile:
        di = std::make_unique<DI_Sterile>(e.fraction, iad, cfgvars.temp);
        break;
      case NCMATData::DynInfo::FreeGas:
        di = std::make_unique<DI_FreeGas>(e.fraction, iad, cfgvars.temp);
        break;
      case NCMATData::DynInfo::VDOSDebye:
        nc_assert_always(data_hasDebyeTemperature);
        di = std::make_unique<DI_VDOSDebye>(e.fraction, iad, cfgvars.temp,
                                            element2DebyeTemp(iad.index));
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
          di = std::make_unique<DI_VDOSImpl>( e.fraction, iad, cfgvars.temp,
                                              std::move(egrid),
                                              VDOSData(vdos_egrid_pair,
                                                       std::move(vdos_density_reg),
                                                       cfgvars.temp,
                                                       iad.data().scatteringXS(),
                                                       iad.data().averageMassAMU()),
                                              std::move(vdos_egrid_orig),
                                              std::move(vdos_density_orig) );
        }
        break;
      case NCMATData::DynInfo::ScatKnl:
        {
          nc_assert( floateq( cfgvars.temp, e.fields.at("temperature").at(0) ) );
          //Prepare scatter kernel object:
          ScatKnlData knldata;
          knldata.temperature = cfgvars.temp;
          knldata.boundXS = iad.data().scatteringXS();
          knldata.elementMassAMU = iad.data().averageMassAMU();
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
          di = std::make_unique<DI_ScatKnlImpl>(e.fraction, iad,
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
  info->setTemperature(cfgvars.temp);
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
    std::map<AtomIndex,std::vector<AtomInfo::Pos>> elem2pos;
    auto itAPE = atompos.end();
    for (auto itAP = atompos.begin() ; itAP != itAPE; ) {
      const std::string current_elementname = itAP->first;
      std::vector<AtomInfo::Pos> positions;
      for ( ; itAP!=itAPE && itAP->first == current_elementname; ++itAP )
        positions.emplace_back(itAP->second[0],itAP->second[1],itAP->second[2]);
      positions.shrink_to_fit();
      const auto& iad = elementname_2_indexedatomdata[current_elementname];
      elem2pos[iad.index] = std::move(positions);
    }

    //==> Calculate element fractions (sanity check consistency with dyninfo
    //fractions, which should have already been validated by data.validate()):
    for (auto& ep : elem2pos) {
      double fr_calc = double(ep.second.size()) / atompos.size();
      //Whether or not the number already exists, we update with these higher precision fractions (e.g. 0.666667
      //specified in file's dyninfo section might be consistent with 2./3. but it is not as precise):
      elem2frac[ep.first] = fr_calc;
    }

    if ( dyninfolist.empty() && data_hasDebyeTemperature ) {
      //Special case: dyninfo not specified, but we know elemental fractions
      //from atom positions, and we have debye temperatures available. We must
      //add DI_VDOSDebye entries for all elements in this case:
      for ( auto& ef: elem2frac )
        dyninfolist.emplace_back(std::make_unique<DI_VDOSDebye>(ef.second,
                                                                *index2iad.at(ef.first.value),
                                                                cfgvars.temp,
                                                                element2DebyeTemp(ef.first)));
    }

    //Check and update fractions on dyninfo objects with numbers in elem2frac,
    //for consistency and increased precision in the scenario where numbers in
    //dyninfos have less precision than those calculated from atomic
    //compositions:
    for (auto& di: dyninfolist) {
      double precise_frac = elem2frac.at(di->atom().index);
      if ( precise_frac != di->fraction() ) {
        nc_assert_always( floateq( di->fraction(), precise_frac ) );//check consistency
        di->changeFraction(precise_frac);//ensure highest precision
      }
    }

    nc_assert_always(elem2pos.size()==elem2frac.size());

    //==> Fill Info::AtomInfo
    for ( auto& ef : elem2frac ) {
      IndexedAtomData iad = *index2iad.at(ef.first.value);
      nc_assert(iad.index==ef.first);
      //Calculate MSDs from Debye temp:
      nc_assert_always(data_hasDebyeTemperature);
      const double debye_temp = element2DebyeTemp(ef.first);
      nc_assert_always(debye_temp>0.0);
      const double msd = debyeIsotropicMSD( debye_temp, cfgvars.temp, iad.data().averageMassAMU() );
      AtomInfo ai;
      ai.debye_temp = data_debyetemp_global ? 0.0 : debye_temp;//TODO: put the value here, even for global DT? Would simplify *a lot*!!!
      ai.mean_square_displacement = msd;
      ai.positions = std::move(elem2pos.at(ef.first));
      ai.number_per_unit_cell = ai.positions.size();//kind of redundant
      ai.atom = std::move(iad);
      info->addAtom( std::move(ai) );
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
    //TODO: Make the info object add it by itself from composition (and only allow it to be specified when composition is not, i.e. .laz files...)
    const AtomData & ad = *index2iad.at(ef.first.value)->atomDataSP;
    xs_abs            += ad.captureXS()  * ef.second;
    xs_freescat       += ad.freeScatteringXS().val * ef.second;
    avr_elem_mass_amu += ad.averageMassAMU() * ef.second;
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
    if(cfgvars.dcutoff==0) {
      //Very simple heuristics here for now to select appropriate dcutoff value
      //(specifically we needed to raise the value for expensive Y2O3/SiLu2O5
      //with ~80/65 atoms/cell):
      cfgvars.dcutoff = ( natoms_per_cell>40 ? 0.25 : 0.1 ) ;
      std::string cmt;
      if (cfgvars.dcutoff>=cfgvars.dcutoffup) {
        //automatically selected conflicts with value of dcutoffup.
        cmt = " (lower than usual due to value of dcutoffup)";
        cfgvars.dcutoff = 0.5*cfgvars.dcutoffup;
      }
      if (verbose)
        std::cout<<"NCrystal::NCMATFactory::automatically selected dcutoff level "<< cfgvars.dcutoff << " Aa"<<cmt<<std::endl;
    }
    if ( cfgvars.dcutoff != -1 ) {
      const double fsquare_cut = 1e-5;//NB: Hardcoded to same value as in .nxs factory
      const double merge_tolerance = 1e-6;
      fillHKL(*info,  cfgvars.dcutoff , cfgvars.dcutoffup, cfgvars.expandhkl, fsquare_cut, merge_tolerance);
    }
  }

  //==> Transfer any custom sections:
  if (!data.customSections.empty()) {
    if (s_NCMATWarnOnCustomSections)
      std::cout<<"NCrystal::NCMATFactory WARNING: Loading NCMAT data which has @CUSTOM_ section(s). This is OK if intended."<<std::endl;
    else if (verbose)
      std::cout<<"NCrystal::NCMATFactory:: Loaded NCMAT data has @CUSTOM_ section(s). This is OK if intended."<<std::endl;
    info->setCustomData(std::move(data.customSections));
  }

  ///////////
  // Done! //
  ///////////

  info->objectDone();
  return info;

}
