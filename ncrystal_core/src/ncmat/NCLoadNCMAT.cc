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

#include "NCrystal/internal/ncmat/NCLoadNCMAT.hh"
#include "NCrystal/internal/ncmat/NCParseNCMAT.hh"
#include "NCrystal/internal/ncmat/NCNCMATData.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/infobld/NCInfoBuilder.hh"
#include "NCrystal/internal/atomdb/NCAtomDBExtender.hh"
#include "NCrystal/internal/extd_utils/NCFillHKL.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/sab/NCScatKnlData.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

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

namespace NCRYSTAL_NAMESPACE {

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

NC::Info NC::loadNCMAT( const FactImpl::InfoRequest& cfg )
{
  NCMATCfgVars ncmatcfgvars;
  ncmatcfgvars.temp      = cfg.get_temp();
  ncmatcfgvars.dcutoff   = cfg.get_dcutoff();
  ncmatcfgvars.dcutoffup = cfg.get_dcutoffup();
  ncmatcfgvars.atomdb    = cfg.get_atomdb_parsed();
  ncmatcfgvars.dataSourceName = cfg.dataSourceName();
  ncmatcfgvars.originalInfoRequest = &cfg;
  return loadNCMAT( cfg.textData(), std::move(ncmatcfgvars) );
}

NC::Info NC::loadNCMAT( NCMATData&& data,
                        NC::NCMATCfgVars&& cfgvars )
{
  const bool verbose = ncgetenv_bool("DEBUGINFO");

  if (verbose) {
    std::ostringstream ss;
    ss<<"loadNCMAT called with ("
      << data.sourceDescription//usually (always?) the same as cfgvars.dataSourceName()
      <<", temp="<<cfgvars.temp
      <<", dcutoff="<<cfgvars.dcutoff
      <<", dcutoffup="<<cfgvars.dcutoffup
      <<", atomdb=";
    if (cfgvars.atomdb.empty()) {
      ss<<"<none>";
    } else {
      for (unsigned i = 0; i < cfgvars.atomdb.size(); ++i) {
        if (i>0)
          ss<<"@";
        ss<<joinstr(cfgvars.atomdb.at(i),":");
      }
    }
    ss<<", dataSourceName="
      << cfgvars.dataSourceName
      <<")";
    Msg::outputMsg(ss.str());
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
  //in ScatKnl sections of dyninfos or @TEMPERATURE sections in NCMAT v7+):
  Temperature dyninfo_temperature{-1.0};
  if (data.hasDynInfo()) {
    for (auto& e : data.dyninfos) {
      if ( e.dyninfo_type != NCMATData::DynInfo::ScatKnl )
        continue;
      Temperature dit{ e.fields.at("temperature").at(0) };
      if ( dyninfo_temperature.get() == -1.0 ) {
        dyninfo_temperature = dit;
      } else {
        if ( !floateq(dyninfo_temperature.get(), dit.get() ) )
          NCRYSTAL_THROW2(BadInput,data.sourceDescription << " incompatible "
                          "temperature values in different @DYNINFO sections");
      }
    }
  }
  Optional<std::pair<Temperature,NCMATData::TemperatureType>> input_temp_request;
  if ( dyninfo_temperature.dbl() > -1.0 ) {
    input_temp_request.emplace( dyninfo_temperature, NCMATData::TemperatureType::Fixed );
  }

  if ( data.temperature.has_value() ) {
    if ( input_temp_request.has_value() ) {
      if ( ! (input_temp_request.value().first == data.temperature.value().first) )
        NCRYSTAL_THROW2(BadInput,data.sourceDescription << " incompatible"
                        " temperature values in @TEMPERATURE and @DYNINFO"
                        " sections.");
    } else {
      input_temp_request = data.temperature.value();
    }
  }

  if ( !input_temp_request.has_value() ) {
    //implicitly defaults to 293.15K
    input_temp_request.emplace( Temperature{293.15}, NCMATData::TemperatureType::Default );
  }

  if ( cfgvars.temp.get() == -1.0 ) {
    //Nothing indicated in cfg strings, go for default/fixed value found in data:
    cfgvars.temp = input_temp_request.value().first;
  } else if ( input_temp_request.value().second == NCMATData::TemperatureType::Fixed
              && !floateq( cfgvars.temp.dbl(), input_temp_request.value().first.dbl() ) ) {
    NCRYSTAL_THROW2(BadInput,data.sourceDescription <<" requested temperature ("<<cfgvars.temp<<")"
                    " is incompatible with temperature ("<<input_temp_request.value().first<<") at which input data is valid.");
  }
  nc_assert_always( cfgvars.temp.dbl() > 0.0
                    && ( input_temp_request.value().second == NCMATData::TemperatureType::Default
                         || floateq( input_temp_request.value().first.dbl(),cfgvars.temp.dbl() ) ) );

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

  InfoBuilder::SinglePhaseBuilder builder;

  //Cache all the hasXXX results, before we start messing them up by
  //std::move'ing stuff out of the data object:
  const bool data_hasDynInfo = data.hasDynInfo();
  const bool data_hasUnitCell = data.hasUnitCell();
  const bool data_hasSpaceGroup = data.hasSpaceGroup();
  const bool data_hasDebyeTemperature = data.hasDebyeTemperature();//@DEBYETEMPERATURE section (must also check @DYNINFO section)
  const bool data_hasDensity = data.hasDensity();

  nc_assert_always( data_hasUnitCell == data.hasAtomPos() && data_hasUnitCell == data.hasCell() );

  //Find Debye temperatures of elements from @DEBYETEMPERATURE section:
  std::map<AtomIndex,DebyeTemperature> perelemdebye_map;
  std::map<AtomIndex,double> elem2msd;

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

  //==> Dynamics (deal with them first as they might provide MSD/DebyeTemp info)
  std::map<AtomIndex,double> elem2frac;
  builder.dynamics = DynamicInfoList();
  //builder.dynamics.value().reserve( elementname_2_indexedatomdata_map.size() );
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
    FactoryJobs jobs_msdcalc;
    bool warn_msd_from_debye_but_vdos_avail = false;
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
        {
          Optional<DebyeTemperature> debye_temp;
          auto itDT = e.fields.find("debye_temp");//can only be there in NCMAT v5+ (and
                                                  //then we are guaranteed that there
                                                  //is no @DEBYETEMPERATURE section)
          if ( itDT != e.fields.end() ) {
            nc_assert(data.version>=5);
            nc_assert(!data_hasDebyeTemperature);
            nc_assert(itDT->second.size()==1);
            debye_temp = DebyeTemperature{ itDT->second.at(0) };
            perelemdebye_map[iad.index] = debye_temp.value();
          }
          if ( !debye_temp.has_value() )
            debye_temp = element2DebyeTemp(iad.index);
          di = std::make_unique<DI_VDOSDebye>(e.fraction, iad, cfgvars.temp, debye_temp.value());
        }
        break;
      case NCMATData::DynInfo::VDOS:
        {
          VectD egrid = getEgrid(e.fields);
          auto vdos_egrid_orig = std::move(e.fields.at("vdos_egrid"));
          auto vdos_density_orig = std::move(e.fields.at("vdos_density"));

          //NB: The following code should remain similar to the equivalent code
          //in ncrystal.cc in ncc::createVDOSDataFromRaw:

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
          //If missing, estimate MSD/DebyeTemp from VDOS (in v4 or later):
          if ( data.version >= 4 ) {
            if ( tryElement2DebyeTemp(iad.index).has_value() ) {
              warn_msd_from_debye_but_vdos_avail = true;//@DEBYETEMPERATURE section takes precedence, but emit warning
            } else {
#if 0
              //Directly:
              double msd = VDOSEval( static_cast<const DI_VDOS*>(di.get())->vdosData() ).getMSD();
              auto debye_temp = debyeTempFromIsotropicMSD( msd, cfgvars.temp, di->atomData().averageMassAMU() );
              perelemdebye_map[iad.index] = debye_temp;
              elem2msd[iad.index] = msd;//record so we don't have to convert back from Debye temp again below
#else
              //Concurrently:
              auto di_raw_ptr = static_cast<const DI_VDOS*>(di.get());
              perelemdebye_map[iad.index] = DebyeTemperature(100.0);//dummy
              elem2msd[iad.index] = -1.0;//dummy
              DebyeTemperature * res_debyetemp_ptr = &perelemdebye_map.find(iad.index)->second;
              double * res_msd_ptr = &elem2msd.find(iad.index)->second;
              Temperature temp = cfgvars.temp;
              jobs_msdcalc.queue([di_raw_ptr,
                                  res_debyetemp_ptr,
                                  res_msd_ptr,
                                  temp]()
              {
                double msd = VDOSEval( di_raw_ptr->vdosData() ).getMSD();//NB: getMSD->calcGamma0 is significant (~10ms) work.
                auto debye_temp = debyeTempFromIsotropicMSD( msd, temp, di_raw_ptr->atomData().averageMassAMU() );
                *res_debyetemp_ptr = debye_temp;
                *res_msd_ptr = msd;
              });
#endif
            }
          }
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
          //TODO: Also try to estimate Debye temperature from SAB?
        }
        break;
      default:
        nc_assert_always(false);//data.validate() should have caught this.
      };
      //Collect here for now (we might adjust the fractions before adding to NCInfo):
      builder.dynamics.value().push_back(std::move(di));
    }
    if (warn_msd_from_debye_but_vdos_avail)
      NCRYSTAL_WARN("NCMATLoader: Loading NCMAT data which has Debye"
                    " temperatures for elements with VDOS curves available"
                    " (this might give sub-optimal MSD values).");
    jobs_msdcalc.waitAll();//make sure msd/debyetemp info is available.
    builder.dynamics.value().shrink_to_fit();
  }

  //==> Temperature:
  builder.temperature = cfgvars.temp;

  if ( data_hasUnitCell ) {
    builder.unitcell = InfoBuilder::UnitCell{};

    //==> Cell layout (structure info)
    StructureInfo& si = builder.unitcell.value().structinfo;
    si.spacegroup = data_hasSpaceGroup ? data.spacegroup : 0;
    si.lattice_a = data.cell.lengths[0];
    si.lattice_b = data.cell.lengths[1];
    si.lattice_c = data.cell.lengths[2];
    si.alpha = data.cell.angles[0];
    si.beta = data.cell.angles[1];
    si.gamma = data.cell.angles[2];
    //Leave si.volume with zero value (InfoBuilder will calculate it).
    si.n_atoms = data.atompos.size();

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

    if ( builder.dynamics.value().empty() && data_hasDebyeTemperature ) {
      //Special case: dyninfo not specified, but we know elemental fractions
      //from atom positions, and we have debye temperatures available. We must
      //add DI_VDOSDebye entries for all elements in this case:
      for ( auto& ef: elem2frac )
        builder.dynamics.value().emplace_back(std::make_unique<DI_VDOSDebye>(ef.second,
                                                                             *index2iad.at(ef.first.get()),
                                                                             cfgvars.temp,
                                                                             element2DebyeTemp(ef.first)));
    }

    nc_assert_always(elem2pos.size()==elem2frac.size());

    //==> Prepare MSD numbers from Debye temperature (some might have been already filled from VDOS's above):

    //Convert Debye temperatures to MSD values (as described in sec. 2.5 of the
    //first NCrystal paper [https://doi.org/10.1016/j.cpc.2019.07.015]):
    for ( auto& di: builder.dynamics.value() ) {
      if ( elem2msd.count(di->atom().index) )
        continue;//Already filled
      auto debyeTemp = tryElement2DebyeTemp(di->atom().index);
      if (!debyeTemp.has_value()) {
        NCRYSTAL_THROW(LogicError,"No MSD source available for element in crystal (this"
                       "  should not happen - other code should have caught this earlier)");
      }
      elem2msd[di->atom().index] = debyeIsotropicMSD( debyeTemp.value(), cfgvars.temp, di->atomData().averageMassAMU() );
    }

    //==> Fill AtomInfo list:
    builder.unitcell.value().atomlist.emplace();
    auto& bldr_atomlist = builder.unitcell.value().atomlist.value();
    //bldr_atomlist.reserve( elem2frac.size() );
    for ( auto& ef : elem2frac ) {
      IndexedAtomData iad = *index2iad.at(ef.first.get());
      nc_assert(iad.index==ef.first);
      nc_assert( elem2msd.count(iad.index) );
      double msd = elem2msd[iad.index];
      Optional<DebyeTemperature> ai_debye_temp;
      if ( !perelemdebye_map.empty() )
        ai_debye_temp = tryElement2DebyeTemp(ef.first);
      if ( !ai_debye_temp.has_value() ) {
        //No direct debye temp was provided in input, but we can estimate one via MSD:
        const auto mass = iad.atomDataSP->averageMassAMU();
        ai_debye_temp = debyeTempFromIsotropicMSD( msd, cfgvars.temp, mass );
      }
      bldr_atomlist.emplace_back( std::move(iad), std::move(std::move(elem2pos.at(ef.first))), ai_debye_temp, msd );
    }

  }//end of unit-cell specific stuff


  if ( data_hasUnitCell ) {
    //Density will be automatically calculated based on unit cell, no need to
    //set explicitly.
  } else {
    //No unit cell, user must have specified density:
    nc_assert_always( data_hasDensity );
    if ( data.density_unit == NCMATData::KG_PER_M3 ) {
      builder.density = Density{ data.density * 0.001 };
    } else {
      nc_assert( data.density_unit == NCMATData::ATOMS_PER_AA3 );
      builder.numberDensity = NumberDensity{ data.density };
    }
  }

  //==> Finally populate HKL list if appropriate:
  if ( data_hasUnitCell ) {
    if(cfgvars.dcutoff==0) {
      //Very simple heuristics here for now to select appropriate dcutoff value
      //(specifically we needed to raise the value for expensive Y2O3/SiLu2O5
      //with ~80/65 atoms/cell):
      cfgvars.dcutoff = ( builder.unitcell.value().structinfo.n_atoms > 40 ? 0.2 : 0.1 ) ;
      std::string cmt;
      if (cfgvars.dcutoff>=cfgvars.dcutoffup) {
        //automatically selected conflicts with value of dcutoffup.
        cmt = " (lower than usual due to value of dcutoffup)";
        cfgvars.dcutoff = 0.5 * cfgvars.dcutoffup;
      }
      if (verbose)
        NCRYSTAL_MSG("NCMATLoader automatically selected dcutoff level "
                     << cfgvars.dcutoff << " Aa"<<cmt);
    }
    if ( cfgvars.dcutoff != -1 ) {
      //HKL planes. Allow on-demand/delayed initialisation.
      builder.hklPlanes.emplace();//Set up uninitialised HKLPlanes struct
      builder.hklPlanes.value().dspacingRange = { cfgvars.dcutoff, cfgvars.dcutoffup };
      using HKLListGenFct = InfoBuilder::SinglePhaseBuilder::HKLPlanes::HKLListGenFct;
      HKLListGenFct genfct = [](const StructureInfo* si,
                                             const AtomInfoList* ai,
                                             PairDD dspacingRangeRequest) -> HKLList
      {
        nc_assert_always(si!=nullptr);
        nc_assert_always(ai!=nullptr);
        FillHKLCfg hklcfg;
        hklcfg.dcutoff = dspacingRangeRequest.first;
        hklcfg.dcutoffup = dspacingRangeRequest.second;
        return calculateHKLPlanes( *si, *ai, std::move(hklcfg) );
      };
      builder.hklPlanes.value().source = genfct;
    }
  }

  //==> Transfer any custom sections:
  if (!data.customSections.empty()) {
    Msg::outputMsg("Loading NCMAT data which has @CUSTOM_ section(s). This is OK if intended.",
                   (s_NCMATWarnOnCustomSections?MsgType::Warning:MsgType::Info));
    builder.customData = std::move(data.customSections);
  }

  //==> State of matter.

  if ( data.stateOfMatter.has_value() ) {
    switch( data.stateOfMatter.value() ) {
    case NCMATData::StateOfMatter::Solid:
      builder.stateOfMatter = Info::StateOfMatter::Solid;
      break;
    case NCMATData::StateOfMatter::Gas:
      builder.stateOfMatter = Info::StateOfMatter::Gas;
      break;
    case NCMATData::StateOfMatter::Liquid:
      builder.stateOfMatter = Info::StateOfMatter::Liquid;
      break;
    default:
      nc_assert_always(false);//should not happen
    }
  }

  //////////////////////////////////////
  // Multiphase support if requested. //
  //////////////////////////////////////

  if ( data.hasOtherPhases() ) {

    using Cfg::CfgManip;
    auto createInfoFromSecondaryCfgStr = []( const Cfg::CfgData& top_request_cfgdata, const std::string& infile_cfgstr ) -> InfoPtr
    {
      nc_assert_always( CfgManip::empty( top_request_cfgdata,
                                         [](Cfg::detail::VarId varid){ return Cfg::varGroup(varid) != Cfg::VarGroupId::Info; } ) );
      //Essentially we want to call the top-level createInfo on the
      //infile_cfgstr, but after applying any Info-variables from the
      //top_request:
      MatCfg cfg(infile_cfgstr);
      cfg.apply(top_request_cfgdata);
      return FactImpl::createInfo(cfg);
    };

    //Figure out fraction of the primary phase and the cfg vars of the top
    //request (which needs to be applied to all phases):
    StableSum sumfrac;
    for ( auto& e : data.otherPhases )
      sumfrac.add(e.first);
    const double frac_this = 1.0 - sumfrac.sum();
    nc_assert_always( frac_this > 0.0 && frac_this < 1.0 );

    Cfg::CfgData top_request_cfgdata_dummy;
    const Cfg::CfgData * top_request_cfgdata(nullptr);
    if ( cfgvars.originalInfoRequest == nullptr ) {
      //No original InfoRequest, setup cfgdata based on cfgvars values:
      CfgManip::set_temp( top_request_cfgdata_dummy, cfgvars.temp );
      CfgManip::set_dcutoff( top_request_cfgdata_dummy, cfgvars.dcutoff );
      CfgManip::set_dcutoffup( top_request_cfgdata_dummy, cfgvars.dcutoffup );
      CfgManip::set_atomdb_parsed( top_request_cfgdata_dummy, cfgvars.atomdb );
      top_request_cfgdata = &top_request_cfgdata_dummy;
    } else {
      top_request_cfgdata = &cfgvars.originalInfoRequest->rawCfgData();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // Create other phases Info objects and finish up (multi phase material with OTHERPHASES)! //
    /////////////////////////////////////////////////////////////////////////////////////////////

    //The phase defined directly in our data becomes accessible as "phasechoice=0":
    builder.dataSourceName = cfgvars.dataSourceName.str() + ";phasechoice=0";

    InfoBuilder::MultiPhaseBuilder mp_builder;
    mp_builder.phases.reserve( data.otherPhases.size() + 1 );
    mp_builder.phases.emplace_back( frac_this, InfoBuilder::buildInfoPtr(std::move(builder)) );
    for ( auto& e : data.otherPhases )
      mp_builder.phases.emplace_back( e.first, createInfoFromSecondaryCfgStr( *top_request_cfgdata, e.second ) );
    return InfoBuilder::buildInfo( std::move(mp_builder) );

  }


  ///////////////////////////////////
  // Done (single phase material)! //
  ///////////////////////////////////

  builder.dataSourceName = std::move(cfgvars.dataSourceName);
  return InfoBuilder::buildInfo(std::move(builder));
}
