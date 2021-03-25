////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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
#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCAtomDBExtender.hh"
#include "NCrystal/internal/NCFillHKL.hh"
#include "NCrystal/internal/NCRotMatrix.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCDebyeMSD.hh"
#include "NCrystal/internal/NCIter.hh"
#include "NCrystal/internal/NCSABUtils.hh"
#include "NCrystal/internal/NCScatKnlData.hh"
#include "NCrystal/internal/NCVDOSEval.hh"
#include "NCrystal/internal/NCDynInfoUtils.hh"
#include <iostream>
#include <cstdlib>

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
                 Temperature temperature,
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

  static std::atomic<bool> s_NCMATWarnOnCustomSections(!ncgetenv_bool("NCMAT_NOWARNFORCUSTOM"));

}

bool NC::getNCMATWarnOnCustomSections()
{
  return s_NCMATWarnOnCustomSections;
}

void NC::setNCMATWarnOnCustomSections(bool bb)
{
  s_NCMATWarnOnCustomSections = bb;
}

NC::Info NC::loadNCMAT( const char * ncmat_file,
                        NC::NCMATCfgVars&& cfgvars )
{
  nc_assert_always(ncmat_file);
  return loadNCMAT( std::string(ncmat_file), std::move(cfgvars) );
}

NC::Info NC::loadNCMAT( const std::string& ncmat_file,
                           NC::NCMATCfgVars&& cfgvars )
{
  return loadNCMAT( FactImpl::createTextData( ncmat_file ),
                    std::move(cfgvars) );
}

NC::Info NC::loadNCMAT( const TextData& inputText,
                        NC::NCMATCfgVars&& cfgvars )
{
  const bool doFinalValidation = false;
  //don't validate at end of the parseNCMATData call, since the loadNCMAT call
  //anyway validates.
  NCMATData data = parseNCMATData( inputText, doFinalValidation);
  return loadNCMAT( std::move(data), std::move(cfgvars) );
}


NC::Info NC::loadNCMAT( const MatInfoCfg& cfg )
{
  cfg.infofactopt_validate({"expandhkl"});//only this infofactopt is supported
  NCMATCfgVars ncmatcfgvars;
  ncmatcfgvars.temp      = cfg.get_temp();
  ncmatcfgvars.dcutoff   = cfg.get_dcutoff();
  ncmatcfgvars.dcutoffup = cfg.get_dcutoffup();
  ncmatcfgvars.expandhkl = cfg.get_infofactopt_flag("expandhkl");
  ncmatcfgvars.atomdb    = cfg.get_atomdb_parsed();
  return loadNCMAT( cfg.textData(), std::move(ncmatcfgvars) );
}

NC::Info NC::loadNCMAT( const MatCfg& cfg )
{
  return loadNCMAT( cfg.createInfoCfg() );
}

NC::Info NC::loadNCMAT( NCMATData&& data,
                        NC::NCMATCfgVars&& cfgvars )
{
  const bool verbose = ncgetenv_bool("DEBUGINFO");

  if (verbose) {
    std::cout<<"NCrystal::loadNCMAT called with ("
             << data.sourceDescription
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
  Temperature input_temperature{-1.0};
  if (data.hasDynInfo()) {
    for (auto& e : data.dyninfos) {
      if ( e.dyninfo_type != NCMATData::DynInfo::ScatKnl )
        continue;
      Temperature dit{ e.fields.at("temperature").at(0) };
      if ( input_temperature.get() == -1.0 ) {
        input_temperature = dit;
      } else {
        nc_assert_always( floateq(input_temperature.get(), dit.get() ) );
      }
    }
  }
  if ( cfgvars.temp.get() == -1.0 ) {
    cfgvars.temp = ( input_temperature.get()==-1.0 ? Temperature{293.15} : input_temperature );
  } else {
    if ( input_temperature.get() != -1.0 && !floateq(input_temperature.get(), cfgvars.temp.get()) )
      NCRYSTAL_THROW2(BadInput,data.sourceDescription <<" specified temperature ("<<cfgvars.temp<<"K)"
                      " is incompatible with temperature ("<<input_temperature<<"K) at which input data is valid.");
  }
  nc_assert_always( cfgvars.temp.get() > 0.0 && ( input_temperature.get()==-1.0 || floateq(input_temperature.get(),cfgvars.temp.get()) ) );

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

  std::map<std::string,IndexedAtomData> elementname_2_indexedatomdata_map;
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
      nc_map_force_emplace( elementname_2_indexedatomdata_map, e.val.second, IndexedAtomData{ std::move(e.val.first), AtomIndex{static_cast<unsigned>(e.idx)} } );
  }
  std::vector<const IndexedAtomData*> index2iad;
  index2iad.resize( elementname_2_indexedatomdata_map.size(), nullptr );
  for( const auto& e : elementname_2_indexedatomdata_map ) {
    nc_assert( e.second.index.get() < elementname_2_indexedatomdata_map.size() );
    index2iad.at(e.second.index.get()) = &e.second;
  }

  auto elementname_2_indexedatomdata = [&elementname_2_indexedatomdata_map](const std::string& s) -> const IndexedAtomData&
  {
    auto it = elementname_2_indexedatomdata_map.find(s);
    if ( it == elementname_2_indexedatomdata_map.end() )
      NCRYSTAL_THROW(LogicError,"NCrystal::LoadNCMAT inconsistency detected in name2index map");
    return it->second;
  };

  NCRYSTAL_DEBUGONLY(for (auto& e: index2iad) { nc_assert_always(e!=nullptr); });

  ////////////////////////
  // Create Info object //
  ////////////////////////

  Info info;

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

  std::map<AtomIndex,DebyeTemperature> perelemdebye_map;
  if ( data.debyetemp_global.has_value() ) {
    //Transfer global value to all entries:
    for ( const auto& e : elementname_2_indexedatomdata_map )
      perelemdebye_map[e.second.index] = data.debyetemp_global.value();
  } else {
    //Set the per-element entries:
    for( const auto& e : data.debyetemp_perelement ) {
      const auto& iad = elementname_2_indexedatomdata(e.first);
      perelemdebye_map[iad.index] = e.second;
    }
  }

  auto tryElement2DebyeTemp = [&perelemdebye_map](const AtomIndex& idx)
  {
    Optional<DebyeTemperature> res;
    auto it = perelemdebye_map.find(idx);
    if ( it != perelemdebye_map.end() )
      res = DebyeTemperature{ it->second };
    return res;
  };

  auto element2DebyeTemp = [&perelemdebye_map](const AtomIndex& idx)
  {
    //only call if there!
    auto it = perelemdebye_map.find(idx);
    nc_assert_always( it != perelemdebye_map.end() );
    return DebyeTemperature{ it->second };
  };

  //==> Dynamics (deal with them first as they might provide MSD info)
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
      const auto& iad = elementname_2_indexedatomdata(e.element_name);
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
                                                       iad.data().scatteringXS(),//(full xs, incoherent approximation)
                                                       iad.data().averageMassAMU()),
                                              std::move(vdos_egrid_orig),
                                              std::move(vdos_density_orig) );
        }
        break;
      case NCMATData::DynInfo::ScatKnl:
        {
          nc_assert( floateq( cfgvars.temp.get(), e.fields.at("temperature").at(0) ) );
          //Prepare scatter kernel object:
          ScatKnlData knldata;
          knldata.temperature = cfgvars.temp;
          knldata.boundXS = iad.data().scatteringXS();//(full xs, incoherent approximation)
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
  info.setTemperature(cfgvars.temp);
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
      info.setStructInfo(si);
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
      const auto& iad = elementname_2_indexedatomdata(current_elementname);
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
                                                                *index2iad.at(ef.first.get()),
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

    //==> Prepare MSD numbers from VDOS or Debye temperature:

    //NB: msd from vdos only allowed in NCMAT v4 or later. If both msd-from-vdos
    //and msd-from-debyetemp is possible, the latter will be preferred (as per
    //the NCMAT doc), but we will emit a warning at the end.

    const bool do_allow_msd_from_vdos = ( data.version >= 4);
    bool warn_msd_from_debye_but_vdos_avail = false;

    std::map<AtomIndex,Optional<double>> elem2msd;

    for (auto& di: dyninfolist) {
      auto& msd = elem2msd[di->atom().index];
      //Possible sources are Debye temp or VDOS:
      auto debyeTemp = tryElement2DebyeTemp(di->atom().index);
      auto di_vdos = do_allow_msd_from_vdos ? dynamic_cast<const DI_VDOS*>(di.get()) : nullptr;
      if ( debyeTemp.has_value() ) {
        if ( di_vdos)
          warn_msd_from_debye_but_vdos_avail = true;
        //Estimate via Debye temperature (as described in sec. 2.5
        //of the first NCrystal paper [https://doi.org/10.1016/j.cpc.2019.07.015]):
        const auto mass = di->atomData().averageMassAMU();
        msd = debyeIsotropicMSD( debyeTemp.value(), cfgvars.temp, mass );
      } else if ( di_vdos ) {
        //Estimate MSD by VDOS integral:
        msd = VDOSEval( di_vdos->vdosData() ).getMSD();
      } else {
        NCRYSTAL_THROW(LogicError,"No MSD source available (this should not"
                       " happen - other code should have caught this earlier)");
      }
    }

    if (warn_msd_from_debye_but_vdos_avail)//Fixme: Make it possible to disable this warning
      std::cout<<"NCrystal::NCMATFactory WARNING: Loading NCMAT data which has Debye temperatures for elements with VDOS curves available (this might give sub-optimal MSD values)."<<std::endl;

    //==> Fill Info::AtomInfo
    for ( auto& ef : elem2frac ) {
      IndexedAtomData iad = *index2iad.at(ef.first.get());
      nc_assert(iad.index==ef.first);
      const auto& msd = elem2msd[iad.index];
      nc_assert_always( msd.has_value() );
      Optional<DebyeTemperature> ai_debye_temp;
      if ( data_hasDebyeTemperature )
        ai_debye_temp = tryElement2DebyeTemp(ef.first);
      if ( !ai_debye_temp.has_value() ) {
        //No direct debye temp was provided in input, but we can estimate one via MSD:
        const auto mass = iad.atomDataSP->averageMassAMU();
        ai_debye_temp = debyeTempFromIsotropicMSD( msd.value(), cfgvars.temp, mass );
      }
      info.addAtom( AtomInfo( std::move(iad), std::move(std::move(elem2pos.at(ef.first))), ai_debye_temp, msd ) );
    }

  }//end of unit-cell specific stuff

  //==> Actually register the dyninfo objects now, since fractions have final values:
  for (auto& di: dyninfolist)
    info.addDynInfo(std::move(di));

  //==> Figure out densities and cross sections:

  SigmaAbsorption xs_abs;
  SigmaFree xs_freescat;
  double avr_elem_mass_amu = 0.;
  for (auto& ef : elem2frac) {
    //TODO: Make the info object add it by itself from composition (and only allow it to be specified when composition is not, i.e. .laz files...)
    const AtomData& ad = index2iad.at(ef.first.get())->atomDataSP;
    xs_abs.dbl()      += ad.captureXS()  * ef.second;
    xs_freescat.dbl() += ad.freeScatteringXS().get() * ef.second;
    avr_elem_mass_amu += ad.averageMassAMU().get() * ef.second;
  }

  info.setXSectAbsorption( xs_abs );
  info.setXSectFree( xs_freescat );

  Density density;
  NumberDensity numberdensity;

  //1e27 in next line converts kg/Aa^3 to g/cm^3:
  const double numberdensity_2_density = 1e27 * avr_elem_mass_amu * constant_dalton2kg;

  if ( cell_volume ) {
    //Calculate from unit cell
    nc_assert_always( natoms_per_cell>0 && !data_hasDensity );
    numberdensity.dbl() = natoms_per_cell / cell_volume;
    density.dbl() = numberdensity.dbl() * numberdensity_2_density;
  } else {
    //No unit cell, user must have specified:
    nc_assert_always( !natoms_per_cell && data_hasDensity );
    if ( data.density_unit == NCMATData::KG_PER_M3 ) {
      density.dbl() = data.density * 0.001;
      numberdensity.dbl() = density.dbl() / numberdensity_2_density;
    } else {
      nc_assert( data.density_unit == NCMATData::ATOMS_PER_AA3 );
      numberdensity.dbl() = data.density;
      density.dbl() = numberdensity.dbl() * numberdensity_2_density;
    }
  }

  info.setDensity( density );
  info.setNumberDensity( numberdensity );

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

      FillHKLCfg hklcfg;
      hklcfg.dcutoff = cfgvars.dcutoff;
      hklcfg.dcutoffup = cfgvars.dcutoffup;
      hklcfg.expandhkl = cfgvars.expandhkl;

      fillHKL( info,  hklcfg );
    }
  }

  //==> Transfer any custom sections:
  if (!data.customSections.empty()) {
    if (s_NCMATWarnOnCustomSections)
      std::cout<<"NCrystal::NCMATFactory WARNING: Loading NCMAT data which has @CUSTOM_ section(s). This is OK if intended."<<std::endl;
    else if (verbose)
      std::cout<<"NCrystal::NCMATFactory:: Loaded NCMAT data has @CUSTOM_ section(s). This is OK if intended."<<std::endl;
    info.setCustomData(std::move(data.customSections));
  }

  ///////////
  // Done! //
  ///////////

  info.objectDone();
  return info;

}
