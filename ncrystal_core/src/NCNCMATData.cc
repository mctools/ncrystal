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

#include "NCrystal/NCNCMATData.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCAtomUtils.hh"
#include "NCrystal/internal/NCIter.hh"
namespace NC = NCrystal;

void NC::NCMATData::DynInfo::validate() const
{
  //Check that required fields were present:
  if ( element_name.empty() )
    NCRYSTAL_THROW(BadInput,"Element name missing");//NB: NCMATData::validate() validates element_name more carefully.
  if ( fraction==-1.0 )
    NCRYSTAL_THROW(BadInput,"Fraction not set");
  if ( dyninfo_type==Undefined )
    NCRYSTAL_THROW(BadInput,"Type of dynamic info is unspecified");
  std::set<std::string> requiredfields,optionalfields;
  if ( dyninfo_type==VDOS ) {
    requiredfields.insert("vdos_egrid");
    requiredfields.insert("vdos_density");
    optionalfields.insert("egrid");
  } else if ( dyninfo_type==ScatKnl ) {
    requiredfields.insert("temperature");
    optionalfields.insert("egrid");
    //Must be exactly one sqw/sab/sab_scaled field:
    bool sk_sab(fields.count("sab")), sk_sab_scaled(fields.count("sab_scaled")), sk_sqw(fields.count("sqw"));
    auto nsk = (sk_sab?1:0) + (sk_sab_scaled?1:0) + (sk_sqw?1:0);
    if ( nsk > 1 )
      NCRYSTAL_THROW(BadInput,"Can not specify more than one of the following fields: \"sab\", \"sab_scaled\", and \"sqw\"");
    if ( nsk != 1 )
      NCRYSTAL_THROW(BadInput,"Must always specify exactly one of fields: \"sab\", \"sab_scaled\", and \"sqw\"");

    const auto sqw_grids = VectS{ "qgrid", "omegagrid"};
    const auto sab_grids = VectS{ "alphagrid", "betagrid"};
    const auto& grids = sk_sqw ? sqw_grids : sab_grids;
    for ( const auto& fn : grids ) {
      auto it = fields.find(fn);
      if ( it == fields.end() )
        NCRYSTAL_THROW2(BadInput,"missing field \""<<fn<<"\" (always required for \""<<(sk_sqw?"sqw":(sk_sab?"sab":"sab_scaled"))<<"\"-scatter kernels)");
      auto nn = it->second.size();
      if ( nn < 5 || nn > 65534 )
        NCRYSTAL_THROW2(BadInput,"Invalid number of entries ("<<nn<<") for field \""<<fn<<"\". Must be at least 5 and at most 65534.");
    }

    //already checked presence, but avoid spurious error further below:
    for ( auto& e : grids )
      requiredfields.insert(e);
    requiredfields.insert(sk_sab?"sab":(sk_sab_scaled?"sab_scaled":"sqw"));

  }
  {
    std::set<std::string>::const_iterator it(requiredfields.begin()),itE(requiredfields.end());
    for (;it!=itE;++it)
      if (!fields.count(*it))
        NCRYSTAL_THROW2(BadInput,"missing field \""<<*it<<"\" (always required for this type of dynamic info)");
  }

  //Check that no unexpected fields were present:
  {
    FieldMapT::const_iterator it(fields.begin()), itE(fields.end());
    for (;it!=itE;++it)
      if ( !requiredfields.count(it->first) && !optionalfields.count(it->first))
        NCRYSTAL_THROW2(BadInput,"Invalid (at least for this type of dynamic info) field \""<<it->first<<"\" specified");
  }


  //Validate specific entries:
  auto valvector = [](const std::string& name, const VectD& v, bool no_negative ) {
    for (auto&e:v) {
      if ( ncisnan(e) || ncisinf(e) || ( no_negative && e<0.0 ) )
        NCRYSTAL_THROW2(BadInput,"invalid entry in "<<name<<" array : "<<e);
    }
  };

  if ( dyninfo_type==VDOS || dyninfo_type==ScatKnl ) {

    //fields common for type=vdos and type=scatknl:
    if ( fields.count("egrid") ) {
      auto& v_e = fields.at("egrid");
      valvector("egrid",v_e,true);
      if ( v_e.size()<10 ) {
        if ( v_e.size() == 3 ) {
          if ( v_e.at(2) != static_cast<unsigned long>(v_e.at(2)) )
            NCRYSTAL_THROW(BadInput,"egrid with 3 entries must have integral number as third entry");
        } else if ( v_e.size() != 1 ) {
          NCRYSTAL_THROW(BadInput,"egrid must have either more than 10 entries or exactly 3 or 1 entries");
        }
      }
    }
  }

 if ( dyninfo_type==VDOS ) {

    //fields specific for type=vdos:
    auto& v_dos = fields.at("vdos_density");
    if ( v_dos.size() < 5 )
      NCRYSTAL_THROW(BadInput,"too few vdos_density parameters");
    valvector("vdos",v_dos, true);

    auto& v_er = fields.at("vdos_egrid");
    if ( v_er.size()!=2 && v_er.size()!=v_dos.size() )
      NCRYSTAL_THROW(BadInput,"vdos_egrid keyword must be followed by exactly two parameters"
                     " or by the same number of parameters as specified for the vdos_density field");
    if ( !(v_er.front()>=1e-5) )
      NCRYSTAL_THROW(BadInput,"invalid vdos_egrid parameters : first point must be at least 1e-5 (0.01 meV).");
    if ( !(v_er.back()<=1e6) )
      NCRYSTAL_THROW(BadInput,"invalid vdos_egrid parameters : last point is excessively large (larger than 1MeV !!).");
    if ( !nc_is_grid(v_er) )
      NCRYSTAL_THROW(BadInput,"invalid vdos_egrid parameters : points do not constitute a grid of increasing values");

  } else if ( dyninfo_type==ScatKnl ) {

    //fields specific for type=scatknl:
    auto& v_t = fields.at("temperature");
    if ( v_t.size()!=1 )
      NCRYSTAL_THROW(BadInput,"temperature keyword not followed by exactly one parameter");
    if ( !(v_t.at(0)>0.0) || !(v_t.at(0)<1e6) )
      NCRYSTAL_THROW(BadInput,"invalid temperature value");
    std::string sa,sb,ssab;
    if (fields.count("alphagrid")) {
      sa = "alphagrid";
      sb = "betagrid";
      ssab = fields.count("sab") ? "sab" : "sab_scaled";
    } else {
      sa = "qgrid";
      sb = "omegagrid";
      ssab = "sqw";
    }

    auto& v_a = fields.at(sa);
    auto& v_b = fields.at(sb);
    auto& v_sab = fields.at(ssab);
    if ( v_a.size() < 5 )
      NCRYSTAL_THROW2(BadInput,"too few "<<sa<<" parameters");
    if ( v_b.size() < 5 )
      NCRYSTAL_THROW2(BadInput,"too few "<<sb<<" parameters");
    if ( v_sab.size()!=v_a.size()*v_b.size() )
      NCRYSTAL_THROW2(BadInput,"number of "<<ssab<<" entries is not (size of "<<sa<<")*(size of "<<sb<<")");
    valvector(sa,v_a, true);
    valvector(sb,v_b, false);
    valvector(ssab,v_sab, false);

  }

}

void NC::NCMATData::unaliasElementNames()
{
  //Unalias D and T markers in NCMAT versions >=3.
  if (version<3)
    return;

  auto doUnalias = [](std::string& name)
  {
    if (name.size()==1) {
      auto c = name[0];
      if ( c == 'D' )
        name = "H2";
      else if ( c=='T' )
        name = "H3";
    }
  };
  for ( auto& e : atompos )
    doUnalias( e.first );
  for ( auto& e : debyetemp_perelement )
    doUnalias( e.first );
  for ( auto& e : dyninfos )
    doUnalias( e.element_name );
}


bool NC::NCMATData::hasCell() const
{
  return cell.lengths[0]!=0. || cell.lengths[1]!=0. || cell.lengths[2]!=0.
    || cell.angles[0]!=0. || cell.angles[1]!=0. || cell.angles[2]!=0.;
}

bool NC::NCMATData::hasUnitCell() const {
  nc_assert( hasCell()==hasAtomPos() && hasDebyeTemperature()==hasCell() );
  return hasCell();
}

void NC::NCMATData::validateCell() const
{
  if (!hasCell())
    return;
  if (cell.lengths[0]==0. && cell.lengths[1]==0. && cell.lengths[2]!=0.)
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" cell section is missing \"lengths\" data");
  if (cell.angles[0]==0. && cell.angles[1]==0. && cell.angles[2]!=0.)
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" cell section is missing \"angles\" data");
  for (int i = 0; i<3; ++i) {
    if ( !(cell.lengths[i]>0.0) || cell.lengths[i]> 1.0e4 )
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid lattice length specified");
    if ( !(cell.angles[i]>0.0) || cell.angles[i]>=180. )
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid lattice angle specified");
    if (ncmax(cell.angles[0],ncmax(cell.angles[1],cell.angles[2]))<=k2Pi)
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid lattice angles specified (perhaps they are in radians instead of the expected degrees?)");
  }
}

void NC::NCMATData::validateAtomDB() const
{
  //Not a complete validation, just a partial one (when used with an
  //AtomDBExtender, the entries will be validated more carefully).
  if (!hasAtomDB())
    return;
  for (auto&& e : enumerate(atomDBLines)) {
    const AtomDBLine& line = e.val;
    nc_assert(!line.empty());
    try {
      validateAtomDBLine( line );
    } catch (Error::BadInput&err) {
      NCRYSTAL_THROW2(BadInput,"Invalid entry in @ATOMDB section in the line: \""<<joinstr(line)<<"\". Error is: "<<err.what());
    }
    if (line.at(0)=="nodefaults") {
      if (e.idx!=0||line.size()!=1)
        NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" \"nodefaults\" keyword in @ATOMDB section can only appear in"
                        " the first line (where it must be alone)");
      continue;
    }
    //Not needed: validateElementName(line.at(0));
  }
}

void NC::NCMATData::validateElementNameByVersion(const std::string& s, unsigned version)
{
  nc_assert_always(version>0&&version<=3);
  AtomSymbol atomsymbol(s);
  if ( atomsymbol.isInvalid() )
    NCRYSTAL_THROW2(BadInput,"Invalid element name \""<<s<<"\"");//invalid in any version
  if (version==3)
    return;//All is supported in the v3.

  //Version-specific tests for older versions:
  if (atomsymbol.isCustomMarker())
    NCRYSTAL_THROW2(BadInput,"Invalid element name \""<<s
                    <<"\" (custom markers X, X1, X2, ..., X99 are only supported from NCMAT v3).");
  //"D" was introduced in v2 as the only new element over v1:
  if (s=="D") {
    if (version==1)
      NCRYSTAL_THROW(BadInput,"Element \"D\" is not supported in NCMAT v1 files (requires NCMAT v2 or later)");
    return;
  }
  if (atomsymbol.isIsotope())
    NCRYSTAL_THROW2(BadInput,"Invalid element name \""<<s<<"\" (general isotope markers are only supported from NCMAT v3).");
  //All ok for the older format:
  return;
}

void NC::NCMATData::validateElementName(const std::string& s) const
{
  try{
    validateElementNameByVersion(s,version);
  } catch (Error::BadInput&e) {
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" "<<e.what());
  }
}

void NC::NCMATData::validateAtomPos() const
{
  if (!hasAtomPos())
    return;
  decltype(atompos)::const_iterator it(atompos.begin()), itE(atompos.end());
  for (;it!=itE;++it) {
    validateElementName(it->first);
    if ( !(it->second[0]>=-1.0) || !(it->second[0]<=1.0)
         || !(it->second[1]>=-1.0) || !(it->second[1]<=1.0)
         || !(it->second[2]>=-1.0) || !(it->second[2]<=1.0) )
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid atomic position detected for element \""<<it->first<<"\" (all position coordinates must be in [-1.0,1.0]");
  }
}

void NC::NCMATData::validateSpaceGroup() const
{
  if (!hasSpaceGroup())
    return;
  if (spacegroup<1||spacegroup>230)
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid spacegroup number (expects a number from 1 to 230)");
}

void NC::NCMATData::validateDebyeTemperature() const
{
  if (!hasDebyeTemperature())
    return;
  if ( debyetemp_global && !debyetemp_perelement.empty() )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" specifies both global and per-element Debye temperatures");
  if ( !(debyetemp_global>=0.0) )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" specifies invalid value of global Debye temperature");
  if ( !debyetemp_perelement.empty() ) {
    std::set<std::string> seen;
    std::vector<std::pair<std::string,double> >::const_iterator it, itE(debyetemp_perelement.end());
    for (it=debyetemp_perelement.begin();it!=itE;++it) {
      validateElementName(it->first);
      if (seen.count(it->first))
        NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" specifies multiple per-element Debye temperatures for element "<<it->first);
      seen.insert(it->first);
      if ( !(it->second>=0.0) )
        NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" specifies invalid value of per-element Debye temperature for element "<<it->first);
    }
  }
}

void NC::NCMATData::validateDensity() const
{
  if ( density == 0 )
    return;
  if ( !(density>0.0) || ncisinf(density) )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" specifies invalid material density in the density section (negative, nan or inf)");
}

void NC::NCMATData::validate() const
{
  if ( ! ( version==1 || version==2 || version==3 ) )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" unsupported NCMAT format version "<<version);

  std::set<std::string> allElementNames;

  ////////////////////////////////////////////////////////////////////////////////
  //validate individual sections:
  validateCell();
  validateAtomPos();
  validateSpaceGroup();
  validateDebyeTemperature();
  validateDensity();
  validateAtomDB();
  const bool hasDebyeTemp = hasDebyeTemperature();
  for (std::size_t i = 0; i<dyninfos.size(); ++i) {
    if ( !hasDebyeTemp && dyninfos.at(i).dyninfo_type==DynInfo::VDOSDebye )
      NCRYSTAL_THROW2(BadInput,"dyninfo sections of type vdosdebye are only allowed when Debye temperatures are also specified.");
    try {
      dyninfos.at(i).validate();
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" problem in dyninfos["<<i<<"]: "<<e.what());
    }
    validateElementName(dyninfos.at(i).element_name);//more careful version-specific validation
    allElementNames.insert(dyninfos.at(i).element_name);
  }

  for (const auto& e : customSections) {
    if (e.first.empty() || ! contains_only(e.first,"ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid custom section name: \""<<e.first
                      <<"\" (must be non-empty and contain only capitalised letters A-Z)");
  }

  ////////////////////////////////////////////////////////////////////////////////
  //validate presence/absence of sections:

  if ( version==1 ) {
    std::string missing=(hasCell()?(hasAtomPos()?(hasDebyeTemperature()?"":"Debye temperatures"):"atom positions"):"cell");
    if (!missing.empty())
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid section configuration for data claiming to be in legacy NCMAT version 1 (missing section with "<<missing<<" info)");
    std::string toomuch=(hasDynInfo()?"dynamic":(hasDensity()?"density":""));
    if (!toomuch.empty())
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" invalid section configuration for data claiming to be in legacy NCMAT version 1 (can not have section with "<<toomuch<<" info)");
  }

  if (hasSpaceGroup()) {
    if (!hasCell())
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" missing Cell information (always required when specifying space group)");
  }

  if ( hasCell() != hasAtomPos() )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" either specify both Cell and Atom Position information, or specify none of them");

  bool hasunitcellinfo = hasCell();
  nc_assert(hasunitcellinfo==hasAtomPos());

  if (!hasunitcellinfo&&!hasDynInfo())
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" missing information about both unit cell and material dynamics so can not provide any interesting physics and can not even deduce the basic material composition");

  if ( hasDensity() && hasunitcellinfo )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" density should never be explicitly specified when unit cell information is present");

  if ( hasDynInfo() || hasDensity() ) {
    //dynamic information is (should be) present. The material is non-crystalline if Cell info is also absent.
    if ( hasDensity() && !hasDynInfo() )
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" direct density specification should never be needed except for files with dynamic information and no unit cell");
    if ( !hasunitcellinfo && !( hasDynInfo() && hasDensity() ) )
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" when not specifying unit cells, dynamic info sections must be present and density must be explicitly specified");
  } else {
    // Crystalline material with no special dynamic info => must provide unit cell as well as Debye temps
    nc_assert(hasunitcellinfo);
    if (!hasDebyeTemperature())
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" missing Debye temperature information for crystalline material");
  }

  if ( !hasunitcellinfo && hasDebyeTemperature() )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" Debye temperature information is only relevant for crystalline materials with a unit cell defined");
  if ( hasunitcellinfo && !hasDebyeTemperature() )
    NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" Debye temperature information is currently required for all crystalline materials with a unit cell defined");

  ////////////////////////////////////////////////////////////////////////////////
  //More detailed checks involving multiple sections.

  // First collect atomposition into elem->fraction map:
  std::map<std::string,unsigned> atompos_elems2count;
  if (!atompos.empty()) {
    decltype(atompos)::const_iterator it(atompos.begin()), itE(atompos.end());
    for (;it!=itE;++it)
      atompos_elems2count[it->first] += 1;
  }

  for (auto& e : atompos_elems2count)
    allElementNames.insert(e.first);

  // Fractions in DYNINFO sections should add up to unity and the same
  // element should not be in two different DYNINFO sections:

  std::map<std::string,double> dyninfo_elems2frac;
  if (!dyninfos.empty()) {
    std::vector<DynInfo>::const_iterator it(dyninfos.begin()), itE(dyninfos.end());
    double dyninfo_totfrac = 0.0;
    for (;it!=itE;++it) {
      if (dyninfo_elems2frac.count(it->element_name))
        NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" same element ("<<it->element_name<<") specified in two different dynamic info sections");
      dyninfo_elems2frac[it->element_name] = it->fraction;
      dyninfo_totfrac += it->fraction;
      //Check that all elements in dyninfo sections appear in atompos section:
      if ( !atompos.empty() && !atompos_elems2count.count(it->element_name) )
        NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" element ("<<it->element_name<<") specified in dynamic info sections does not appear in the list of atom positions");
    }
    double dist = ncabs(1.0-dyninfo_totfrac);
    if (dist>1e-9)
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" the fractions specified in the dyninfo sections do not add up to unity"
                      <<(dist<1e-3?" (they almost do, but more precision is required - note that it is possible"
                         " to specify fractions like 2/3 in NCMAT files instead of numbers like 0.666667)":""));

    //Check that all elements in atompos section appear in a dyninfo section:
    std::map<std::string,unsigned>::const_iterator itA(atompos_elems2count.begin()), itAE(atompos_elems2count.end());
    for (;itA!=itAE;++itA) {
      if (!dyninfo_elems2frac.count(itA->first))
        NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" element ("<<itA->first<<") appearing in list of atom"
                        " positions is not specified in dynamic info sections (either remove all dynamic info"
                        " sections or make sure there is a complete set of them)");
      // Check that fractions in DYNINFO sections are be consistent with
      // fractions inferred from ATOMPOSITIONS sections:
      double calcfraction = itA->second*1.0/atompos.size();
      double dyninfofraction = dyninfo_elems2frac[itA->first];
      double distf = ncabs(calcfraction-dyninfofraction);
      if (distf>1e-9) {
        std::ostringstream ss;
        ss << sourceFullDescr<<" the fraction for "<<itA->first<<" specified in the dyninfo section differs from the"
           << " same fraction calculated from the number of each atom appearing in the list of atompositions";
        if (distf<1e-3)
          ss << " (they almost do, but more precision is required - note that it is possible to specify"
            " fractions like "<<itA->second<<"/"<<atompos.size()<<" in NCMAT files instead of numbers like "<<dyninfofraction<<")";
        NCRYSTAL_THROW(BadInput,ss.str());
      }
    }
  }

  // Consistency between elements in ATOMPOSITIONS/DYNINFO sections and those with local Debye
  // temperatures:

  if ( !debyetemp_perelement.empty() ) {
    for (auto& e : debyetemp_perelement)
      allElementNames.insert(e.first);

    std::vector<std::pair<std::string,double> >::const_iterator itD(debyetemp_perelement.begin()), itDE(debyetemp_perelement.end());
    std::string refsec;
    if ( !atompos_elems2count.empty() ) {
      for (;itD!=itDE;++itD)
        if ( !atompos_elems2count.count(itD->first) )
          break;
      for (const auto& e : atompos_elems2count) {
        if (!std::any_of(debyetemp_perelement.cbegin(), debyetemp_perelement.cend(),
                         [&e](const std::pair<std::string,double>& e2){ return e.first==e2.first; })) {
          NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" Per-element Debye temperature specified for some elements"
                          " but missing for element "<<e.first<<" which occurs in @ATOMPOSITIONS");
        }
      }
    } else if ( !dyninfo_elems2frac.empty() ) {
      for (;itD!=itDE;++itD)
        if ( !dyninfo_elems2frac.count(itD->first) )
          break;
    }
    if (itD!=itDE)
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" element ("<<itD->first
                      <<") appearing in list of per-element Debye temperatures is not present in "
                      <<((!atompos_elems2count.empty())?"the list of atom positions":"any of the dynamic info sections"));
  }

  //Check consistency between temperature values specified in different parts of
  //the file (so far only dyninfo-sab has temperatures):
  double di_temp = 0.0;
  for (auto& di: dyninfos) {
    if ( di.dyninfo_type == DynInfo::ScatKnl ) {
      double tt = di.fields.at("temperature").at(0);
      nc_assert(tt>0.0);//already checked
      if ( di_temp>0.0 && di_temp != tt )
        NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" temperature values specified in different dynamic info sections are different")
    }
  }

  //Check that any custom markers used in atompos/debyetemp/dyninfo sections are defined in atomdb section:
  std::set<std::string> atomdb_custommarkers;
  for(auto& e : atomDBLines) {
    if (AtomSymbol(e.at(0)).isCustomMarker())
      atomdb_custommarkers.insert(e.at(0));
  }
  for (auto& e : allElementNames) {
    if ( AtomSymbol(e).isCustomMarker() && !atomdb_custommarkers.count(e))
      NCRYSTAL_THROW2(BadInput,sourceFullDescr<<" custom marker \""<<e<<"\" is used but has no definition in the @ATOMDB section")
  }
}
