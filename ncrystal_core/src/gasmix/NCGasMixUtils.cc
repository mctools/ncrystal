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

#include "NCrystal/internal/gasmix/NCGasMixUtils.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/atomdb/NCAtomDBExtender.hh"
#include <sstream>
namespace NC = NCrystal;
namespace NCGM = NCrystal::GasMix;

namespace NCRYSTAL_NAMESPACE {
  namespace GasMix {
    namespace {
      using PairSD = std::pair<std::string,double>;
      static constexpr Pressure default_pressure = Pressure::one_atm();
      static constexpr Temperature default_temperature = Temperature{ 293.15 };
      struct Unit {
        const char* name;
        double scale;
        double offset = 0.0;
        constexpr Unit( const char * name_,
                        double scale_,
                        double offset_ = 0.0 ) : name(name_), scale(scale_), offset(offset_) {}
      };
      static constexpr Unit pressure_units[6] = {
        Unit{ "Pa",   1.0 },//Be aware: Element 91 has symbol "Pa", so we should check for chem. formulas before pressure!
        Unit{ "mbar", 1.0e2 },
        Unit{ "hPa",  1.0e2 },
        Unit{ "kPa",  1.0e3 },
        Unit{ "bar",  1.0e5 },
        Unit{ "atm",  Pressure::one_atm_raw() }
      };
      static constexpr Unit density_units[2] = {
        Unit{ "gcm3", 1.0 },
        Unit{ "kgm3",  1.0e-3 }
      };
      static constexpr Unit temperature_units[3] = {
        Unit{ "K", 1.0, 0.0 },
        Unit{ "C", 1.0, 273.15 },
        Unit{ "F", 5.0/9.0, (273.15-32.*5.0/9.0) }
      };

      std::string toStringWithBestUnit( double value, Span<const Unit> units )
      {
        std::string best;
        for ( auto& unit : units ) {
          std::ostringstream ss;
          ss << fmt( (value-unit.offset) / unit.scale ) << unit.name;
          std::string res = ss.str();
          if ( best.empty() || best.size() > res.size() )
            best = res;
        }
        return best;
      }
      class fmtFlexUnit {
      private:
        std::string m_res;
      public:
        fmtFlexUnit( double value, Span<const Unit> units ) : m_res( toStringWithBestUnit(value,units) ) {}
        const std::string& result() const noexcept { return m_res; }
      };
      inline std::ostream& operator<<(std::ostream& os, const fmtFlexUnit& val) { return os << val.result(); }
    }
  }
}

std::string NCGM::requestToString( const GasMixRequest& gm )
{
  std::ostringstream ss;
  nc_assert_always( isOneOf(gm.fracType,
                            GasMixRequest::FracType::MolarFraction,
                            GasMixRequest::FracType::MassFraction) );

  nc_assert_always( !gm.components.empty() );
  const bool singleCompFormula( gm.components.size() == 1 && gm.components.front().first == 1.0 );
  bool first_comp(true);
  for ( auto& c : gm.components ) {
    if ( !singleCompFormula ) {
      if (!first_comp)
        ss << '+';
      first_comp=false;
      ss << fmt( c.first ) << 'x';
    }
    streamSimpleChemicalFormula(ss,c.second);
  }
  if ( gm.fracType == GasMixRequest::FracType::MassFraction )
    ss << "/massfractions";
  if ( gm.relativeHumidity > 0.0 )
    ss << '/' << fmt(gm.relativeHumidity) << "relhumidity";
  if ( !floateq(gm.temperature.dbl(),default_temperature.dbl(),1e-15,0.0) )
    ss << '/' << fmtFlexUnit(gm.temperature.dbl(),temperature_units);
  nc_assert_always( !gm.targetPresDens.empty() );
  if ( gm.targetPresDens.has_value<Pressure>() ) {
    auto pressure = gm.targetPresDens.get<Pressure>();
    if ( !floateq(pressure.dbl(),default_pressure.dbl(),1e-15,0.0) )
      ss << '/' << fmtFlexUnit(pressure.dbl(),pressure_units);
  } else {
    nc_assert( gm.targetPresDens.has_value<Density>() );
    auto density = gm.targetPresDens.get<Density>();
    ss << '/' << fmtFlexUnit(density.dbl(),density_units);
  }
  for ( auto& e : gm.atomDBLines )
    ss << '/' << joinstr( e, "_" );
  return ss.str();
}

NCGM::GasMixRequest NCGM::requestFromString( const std::string& request )
{
  auto request_sv = StrView(request);

  const char * errmsgprefix = "Syntax error in gasmix specification: ";
  if (!request_sv.isSimpleASCII(StrView::AllowTabs::No,StrView::AllowNewLine::No))
    NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Does not contain only simple ASCII characters :\""<<request<<"\"");
  if (request_sv.contains_any("\"'$@\\*<>"))
    NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Contains forbidden characters :\""<<request<<"\"");

  Optional<Density> density_req;
  Optional<Pressure> pressure_req;
  Optional<double> rel_humidity;
  Optional<Temperature> temperature;
  Optional<std::string> formula;
  Optional<GasMixRequest::FracType> fracType;
  SmallVector<VectS,5> atomdb_lines;

  auto extractValue = [&errmsgprefix]( Unit unit, StrView str ) -> Optional<double>
  {
    if (!str.endswith(unit.name))
      return NullOpt;
    auto val = str.substr(0,str.size()-StrView(unit.name).size()).rtrimmed().toDbl();
    if (!val.has_value())
      NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Could not parse value from: \""<<str<<"\"");
    val.value() *= unit.scale;
    val.value() += unit.offset;
    return val;
  };
  auto handlePressure = [&errmsgprefix,&pressure_req,&extractValue]( StrView entry ) -> bool
  {
    for ( auto& pu : pressure_units ) {
      auto v = extractValue(pu,entry);
      if ( v.has_value() ) {
        if ( !(v.value()>0.0) || !(v.value()<1e10) )
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Pressure value out of bounds: "<<entry);
        if ( pressure_req.has_value() )
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple pressure entries are not allowed.");
        pressure_req = Pressure{ v.value() };
        return true;
      }
    }
    return false;
  };
  auto handleDensity = [&errmsgprefix,&density_req,&extractValue]( StrView entry ) -> bool
  {
    for ( auto& pu : density_units ) {
      auto v = extractValue(pu,entry);
      if ( v.has_value() ) {
        if ( !(v.value()>0.0) || !(v.value()<1e10) )
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Density value out of bounds: "<<entry);
        if ( density_req.has_value() )
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple density entries are not allowed.");
        density_req = Density{ v.value() };
        return true;
      }
    }
    return false;
  };
  auto handleRelHumidity = [&errmsgprefix,&rel_humidity,&extractValue]( StrView entry ) -> bool
  {
    auto v = extractValue(Unit{"relhumidity",1.0},entry);
    if ( v.has_value() ) {
      if ( !(v.value()>=0.0) || !(v.value()<=1.0) )
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Relative humidity value out of bounds: "<<entry);
      if ( rel_humidity.has_value() )
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple relative humidity entries are not allowed.");
      rel_humidity = v.value();
      return true;
    } else {
      return false;
    }
  };
  auto handleFormula = [&errmsgprefix,&formula]( StrView entry ) -> bool
  {
    //Due to cfg-string rules, we can't use '*' sign for multiplication (as in the
    //nicely readable 0.3*CO2+0.7*Ar). We instead use a lower-case 'x'. Note that no
    //element name contains a lower case x, and only Xe contains an X at all (but
    //always upper case). We can detect + parse a formula by it never containing any underscores, and then:
    //
    // "gasmix::0.5xCO2+0.5xAr"  <- contains '+' and 'x', split on '+' then 'x'
    // "gasmix::1.0xCO2"  <-  no '+', but a lower case 'x' followed by a chemical formula.
    // "gasmix::CO2"  <- no '+', but starts with capital letter and consists of chemical formula
    if ( entry.contains('_') )
      return false;
    auto ok = [&entry,&formula,&errmsgprefix]()
    {
      if ( formula.has_value() )
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple gas mixture formulas specified.");
      formula = entry.to_string();
      return true;
    };
    if ( entry == "air" )
      return ok();
    const bool has_x = entry.contains('x');
    const bool has_plus = entry.contains('+');

    if ( has_x && has_plus )
      return ok();//Assume something like "0.5xCO2+0.5xAr"

    if ( has_x && !has_plus ) {
      //Accept formulas like "1.0xCO2" (single component gas)
      auto parts = entry.splitTrimmed('x');
      if ( parts.size()==2 && parts.front().toDbl().has_value() && tryDecodeSimpleChemicalFormula( parts.back().to_string() ).has_value() )
        return ok();
    }
    if ( !has_x && !has_plus && tryDecodeSimpleChemicalFormula( entry.to_string() ).has_value() ) {
      //Accept formulas like "CO2" (single component gas). Prefix with
      //"1.0x" to make handling of code easier.
      ok();
      formula = std::string("1.0x")+formula.value();
      return true;
    }
    return false;
  };
  auto handleFracType = [&errmsgprefix,&fracType]( StrView entry ) -> bool
  {
    Optional<GasMixRequest::FracType> ft;
    if ( entry == "massfractions" ) {
      ft.emplace(GasMixRequest::FracType::MassFraction);
    } else if ( entry == "molarfractions" ) {
      ft.emplace(GasMixRequest::FracType::MolarFraction);
    }
    if ( !ft.has_value() )
      return false;
    if (fracType.has_value() )
      NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple fraction type (massfractions/molarfractions) keywords are not allowed.");
    fracType = ft.value();
    return true;
  };
  auto handleTemp = [&errmsgprefix,&temperature,&extractValue]( StrView entry ) -> bool
  {
    for ( auto& tu : temperature_units ) {
      auto v = extractValue(tu,entry);
      if ( v.has_value() ) {
        if ( !(v.value()>0.0) || !(v.value()<1e10) )
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Temperature value out of bounds: "<<entry);
        if ( temperature.has_value() )
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple temperature entries are not allowed.");
        temperature = Temperature{ v.value() };
        return true;
      }
    }
    return false;
  };
  auto handleAtomDB = [&errmsgprefix,&atomdb_lines]( StrView e ) -> bool
  {
    constexpr static int ncmat_version = AtomDBExtender::latest_version;//latest NCMAT version
    //NB: Code duplicated between here and NCQuickFact.cc!!
    std::string tmp = e.to_string();
    //spaces, semicolons, or double underscores are disallowed as
    //separators, but we first remove them to provide a better error msg
    //for users who use them by mistake:
    if ( contains(tmp,':') )
      strreplace(tmp, ":","_");
    if ( contains(tmp,' ') )
      strreplace(tmp, " ","_");
    while ( contains(tmp,"__") )
      strreplace(tmp, "__","_");
    auto parts = split2(tmp,0,'_');
    auto np = parts.size();
    bool should_be_atomdb_line = ( ( np >= 3 && parts.at(1) == "is" && ( np == 3 || np % 2 == 0 ) )
                                   || ( np == 5 && endswith(parts.at(1),'u') ) );
    if (should_be_atomdb_line) {
      if ( e.contains(':') || e.contains(' ') )
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Do not use spaces or semicolons as separators in ATOMDB sections (use single underscores instead): "<<e<<"\".");
      if ( e.contains("__") )
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Do not use double underscores as separators in ATOMDB sections: "<<e<<"\".");
      try {
        validateAtomDBLine(parts, ncmat_version );
      } catch ( Error::BadInput& err ) {
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Invalid ATOMDB syntax: "<<err.what()<<"\".");
      }
      atomdb_lines.push_back(std::move(parts));
      return true;
    } else {
      return false;
    }
  };

  for ( auto e : request_sv.splitTrimmedNoEmpty('/') ) {
    if ( !handleFormula(e)  //handleFormula+atomdb first since element names
         && !handleAtomDB(e)//clashes with units for other fields (i.e Pa, K, C, ...)
         && !handleFracType(e)
         && !handlePressure(e)
         && !handleDensity(e)
         && !handleRelHumidity(e)
         && !handleTemp(e) )
      NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Unrecognised syntax: \""<<e<<"\".");
  }

  //Sanity check settings and set defaults:
  if ( density_req.has_value() && pressure_req.has_value() )
    NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Can not specify both pressure and density: \""<<request<<"\".");
  if ( ! formula.has_value() )
    NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Missing gas mixture formula specification: \""<<request<<"\".");
  if ( formula.value() == "air" ) {
    if ( fracType.has_value() && fracType.value() != GasMixRequest::FracType::MolarFraction )
      NCRYSTAL_THROW2(BadInput,errmsgprefix<<"\"air\" keyword only works with molar fractions");
    fracType.emplace( GasMixRequest::FracType::MolarFraction );
    //Air, from table at
    //https://www.engineeringtoolbox.com/air-composition-d_212.html (excluding
    //minor contributions mentioned below the table). Adjusted N2 and O2 content
    //slightly to get fractions to sum to unity (original values were 0.78084*N2+0.20946*O2+...):
    formula = "0.7807712xN2+0.20945xO2+0.00934xAr+0.000412xCO2+0.00001818xNe+0.00000524xHe+0.00000179xCH4+0.0000010xKr+0.0000005xH2+0.00000009xXe";
  }

  GasMixRequest result;
  result.fracType = fracType.value_or( GasMixRequest::FracType::MolarFraction );
  result.relativeHumidity = rel_humidity.value_or( 0.0 );
  result.temperature = temperature.value_or( default_temperature );
  if ( density_req.has_value() )
    result.targetPresDens = density_req.value();
  else
    result.targetPresDens = pressure_req.value_or( default_pressure );

  result.atomDBLines.reserve( atomdb_lines.size() );
  for ( auto& e : atomdb_lines )
    result.atomDBLines.push_back( std::move(e) );

  nc_assert(formula.has_value());
  StableSum fracsum;
  for ( auto e : StrView(formula.value()).splitTrimmed('+') ) {
    auto parts = e.splitTrimmed('x');
    double frac(-1.0);
    Optional<DecodedChemForm> cf;
    if ( parts.size() == 2 ) {
      frac = parts.at(0).toDbl().value_or(-1.0);
      cf = tryDecodeSimpleChemicalFormula( parts.at(1).to_string() );
    }
    if ( !cf.has_value() || !(frac>=0.0) || !(frac<=1.0) )
      NCRYSTAL_THROW2(BadInput,errmsgprefix<<"invalid component entry: \""<<e<<"\"");
    if ( frac > 0.0 ) {
      fracsum.add(frac);
      result.components.emplace_back( frac, std::move(cf.value()));
    }
  }

  const double totfrac = fracsum.sum();
  if ( totfrac != 1.0 ) {
    if ( !floateq(1.0,totfrac,1e-6,0.0) )
      NCRYSTAL_THROW2(BadInput,errmsgprefix<<"fractions do not sum up to unity");
    //snap to 1.0:
    for ( auto& e : result.components )
      e.first /= totfrac;
  }
  return result;
}

NCGM::GasMixResult NCGM::analyseGasMixRequest( const GasMixRequest& req )
{
  //1) Set up atom db:
  AtomDBExtender atomdb;//<-- the database
  for ( auto& a : req.atomDBLines ) {
    if ( a.size() == 1 && a.front() == "nodefaults" )
      NCRYSTAL_THROW(BadInput,"nodefaults keyword not allowed in gas mixture atomdb");
    atomdb.addData( a );
  }

  //2) Handle relativeHumidity:
  if ( req.relativeHumidity > 0 ) {
    //First find partial pressure from the water, and the mass of one water molecule:
    auto Pwater = Pressure{ req.relativeHumidity * saturatedVapourPressureOfWater( req.temperature ).dbl() };

    //Guard against heavy enrichments (if we want to improve this we could turn
    //to https://srd.nist.gov/jpcrdreprint/1.555627.pdf):
    auto atomH = atomdb.lookupAtomData("H");
    auto atomO = atomdb.lookupAtomData("O");
    const double threshold = 0.15;//somewhat arbitrarily chosen...
    if ( ncabs( atomH->averageMassAMU().dbl()/1.00784 - 1.0 ) > threshold
         || ncabs( atomO->averageMassAMU().dbl()/15.999 - 1.0 ) > threshold )
      NCRYSTAL_THROW(BadInput,"The saturated vapour pressure formula needed to decode relative humidity has not been validated for the case of H or O atoms being heavily enriched.");

    //Find partial density from the water:
    double water_mol_per_m3 = Pwater.dbl() / ( constant_gas_R * req.temperature.dbl() );
    double water_molecule_grampermole = constant_dalton2gpermol * ( 2.0 * atomH->averageMassAMU().dbl() + atomO->averageMassAMU().dbl() );
    Density RhoWater{ 1e-6 * water_molecule_grampermole * water_mol_per_m3 };

    //Figure out partial P/rho of the rest of the gas:

    //Handle request via recursion, removing the partial pressure or density
    //associated with the water:
    GasMixRequest req2 = req;
    req2.relativeHumidity = 0.0;
    if ( req.targetPresDens.has_value<Pressure>() ) {
      auto Ptot = req.targetPresDens.get<Pressure>();
      Pressure PNonWater{ Ptot.dbl() - Pwater.dbl() };
      if ( PNonWater.dbl() <= 0.0 )
        NCRYSTAL_THROW2(BadInput,"Requested total pressure of " << Ptot
                        <<" must be higher than the partial pressure ("<<Pwater
                        <<") induced by the water molecules at the requested T="<<req.temperature
                        <<" and rel. humidity of "<<req.relativeHumidity);
      req2.targetPresDens = PNonWater;
    } else {
      nc_assert( req.targetPresDens.has_value<Density>() );
      auto Rhotot = req.targetPresDens.get<Density>();
      Density RhoNonWater{ Rhotot.dbl() - RhoWater.dbl() };
      if ( RhoNonWater.dbl() <= 0.0 )
        NCRYSTAL_THROW2(BadInput,"Requested total density of " << Rhotot
                        <<" must be higher than the partial density ("<<RhoWater
                        <<") of the water molecules at the requested T="<<req.temperature
                        <<" and rel. humidity of "<<req.relativeHumidity);
      req2.targetPresDens = RhoNonWater;
    }

    GasMixResult res = analyseGasMixRequest( req2 );
    //Make sure there was not water added explicitly as well in the component list:
    auto cf_h2o = decodeSimpleChemicalFormula("H2O");
    for ( auto& e : res.components )
      if ( e.second == cf_h2o )
        NCRYSTAL_THROW2(BadInput,"Water (H2O) molecules should not be added explicitly in the component list when setting a non-zero relhumidity value");

    //Now add back the water:
    const double molarfrac_water = Pwater.dbl() / ( Pwater.dbl() + res.pressure.dbl() );
    for ( auto& e : res.components )
      e.first *= ( 1.0 - molarfrac_water );

    res.components.emplace_back( molarfrac_water, cf_h2o );
    std::stable_sort(res.components.begin(),res.components.end());

    auto addElemIfAbsent = [&res]( const AtomSymbol& atom, const AtomDataSP& atomdatasp )
    {
      for ( auto& e : res.atomDB )
        if ( e.first == atom )
          return;
      res.atomDB.emplace_back(atom,atomdatasp);
    };
    addElemIfAbsent( AtomSymbol("H"), atomH );
    addElemIfAbsent( AtomSymbol("O"), atomO );
    std::stable_sort( res.atomDB.begin(), res.atomDB.end() );
    res.pressure.dbl() += Pwater.dbl();
    res.density.dbl() += RhoWater.dbl();
    //And avoid numerical issues:
    if ( req.targetPresDens.has_value<Pressure>() ) {
      nc_assert( floateq( res.pressure.dbl(), req.targetPresDens.get<Pressure>().dbl(), 1e-6, 1e-10 ) );
      res.pressure = req.targetPresDens.get<Pressure>();
    } else {
      nc_assert( floateq( res.density.dbl(), req.targetPresDens.get<Density>().dbl(), 1e-6, 1e-10 ) );
      res.density = req.targetPresDens.get<Density>();
    }
    return res;
  }

  nc_assert(req.relativeHumidity==0.0);

  //Temperature field of result is simply as requested:
  NCGM::GasMixResult res;
  res.temperature = req.temperature;

  //Normalise components:
  auto normaliseComponents = []( const ComponentList& comps )
  {
    //merge existing entries:
    ComponentList newcomps;
    for ( auto& e : comps ) {
      if ( !(e.first>0.0) )
        continue;//skip frac=0 entries
      auto norm_cf = normaliseSimpleChemicalFormula( e.second );
      bool found(false);
      for ( auto& o : newcomps ) {
        if ( o.second == norm_cf ) {
          o.first += e.first;
          found = true;
          break;
        }
      }
      if (!found)
        newcomps.emplace_back( e.first, norm_cf );
    }
    //sort:
    std::stable_sort(newcomps.begin(),newcomps.end());
    return newcomps;
  };
  auto components = normaliseComponents( req.components );

  //Find molecular masses and fill the atomDB field at the same time:
  SmallVector<AtomMass,6> molecular_masses;
  auto lookupAtom = [&atomdb,&res]( const AtomSymbol& atom )
  {
    //simple linear search for simplicity:
    for ( auto& e : res.atomDB )
      if ( e.first == atom )
        return e.second;
    nc_assert(atom.isElement());
    auto atomdata = atomdb.lookupAtomData(elementZToName(atom.Z()));
    res.atomDB.emplace_back( atom, atomdata );
    std::stable_sort(res.atomDB.begin(),res.atomDB.end());
    return atomdata;
  };

  for ( auto& ecomp : components ) {
    StableSum mass;
    for ( auto& e : ecomp.second)
      mass.add( e.first * lookupAtom(e.second)->averageMassAMU().dbl() );
    molecular_masses.push_back( AtomMass{ mass.sum() } );
  }

  //Prepare res.components with molar-fractions:

  nc_assert(res.components.empty());
  if ( req.fracType == GasMixRequest::FracType::MassFraction ) {
    //Convert mass fractions to molar fractions!
    StableSum mwsum;
    for ( auto i : ncrange( molecular_masses.size() ) ) {
      const auto& mass_frac = components.at(i).first;
      nc_assert( molecular_masses.at(i).dbl() > 0.0 );
      const double molar_weight = mass_frac / molecular_masses.at(i).dbl();
      res.components.emplace_back( molar_weight, components.at(i).second );
      mwsum.add( molar_weight );
    }
    nc_assert( mwsum.sum() > 0.0 );
    const double mwnorm( 1.0 / mwsum.sum() );
    for ( auto& e : res.components )
      e.first *= mwnorm;
    components.clear();
  } else {
    //already molar fractions:
    nc_assert( req.fracType == GasMixRequest::FracType::MolarFraction );
    res.components = std::move(components);
  }

  //Average molecular mass:
  nc_assert( molecular_masses.size() == res.components.size() );
  StableSum massSum;
  for ( auto i : ncrange( molecular_masses.size() ) ) {
    massSum.add( res.components.at(i).first * molecular_masses.at(i).dbl() );
  }
  const AtomMass averageMass{ massSum.sum() };

  //total rho and P are related by the ideal gas law, using the total molar
  //density (n/V) and average molecular mass. Thus: rho = Mavr*n/V = Mavr*P/RT =
  //(Mavr/RT)*P = c_rhoP * P, where we can calculate c_rhoP in convenient units:
  const double c_rhoP = 1e-6 * averageMass.dbl() * constant_dalton2gpermol / ( constant_gas_R * req.temperature.dbl() );//in [g/(cm3*Pa)]
  nc_assert_always( c_rhoP > 0.0 );

  //Now, calculate final rho and P:
  if ( req.targetPresDens.has_value<Pressure>() ) {
    res.pressure = req.targetPresDens.get<Pressure>();
    res.density = Density{ c_rhoP * res.pressure.dbl() };
  } else {
    nc_assert( req.targetPresDens.has_value<Density>() );
    res.density = req.targetPresDens.get<Density>();
    res.pressure = Pressure{ res.density.dbl() / c_rhoP };
  }

  //We might need to resort, but we do it here at the end (since the
  //res.components and molecular_masses vectors should stay in same order
  //above):
  std::stable_sort( res.components.begin(),res.components.end() );

  return res;
}

NC::Pressure NCGM::saturatedVapourPressureOfWater( Temperature temp )
{
  temp.validate();
  //Using a formula from https://doi.org/10.1175/JAMC-D-17-0334.1 (seems to be
  //at least precise for -100C to 200C, but possibly also giving somewhat
  //meaningful results at temperatures outside this range.)
  const double T_C = temp.dbl() - 273.15;//celcius
  if ( T_C >= 0.0 ) {
    return Pressure{ std::exp(34.494-4924.99/(T_C+237.1)) / std::pow(T_C+105.0,1.57) };
  } else {
    return Pressure{ std::exp(43.494-6545.8/(T_C+278.0)) / ncsquare(T_C+868.0) };
  }
}

std::ostream& NCGM::operator<<( std::ostream& os, const GasMixResult& gm )
{
  os << "GasMixResult{T="
     << gm.temperature
     << ", P=" << fmtFlexUnit( gm.pressure.dbl(), pressure_units )
     << ", Rho=" << fmtFlexUnit( gm.density.dbl(), density_units ) << ';';
  bool first(true);
  for ( auto& e : gm.components ) {
    if (!first)
      os<<'+';
    first = false;
    os << fmt(e.first)<<'x';
    streamSimpleChemicalFormula( os, e.second );
  }
  os << "}";
  return os;
}

NCGM::AtomicComponentList NCGM::flattenComponentList( const GasMixResult& gm )
{
  NCGM::AtomicComponentList res;
  auto addEntry = [&res]( double weight, const AtomSymbol& atom )
  {
    for ( auto& e : res ) {
      if ( e.second == atom ) {
        e.first += weight;
        return;
      }
    }
    res.emplace_back(weight,atom);
  };
  for ( auto& e : gm.components )
    for ( auto& N_atom : e.second )
      addEntry( N_atom.first*e.first, N_atom.second );

  StableSum sum;
  for ( auto& e : res )
    sum.add( e.first );
  nc_assert( sum.sum() > 0.0 );
  const double scale( 1.0 / sum.sum() );
  for ( auto& e : res )
    e.first *= scale;

  std::stable_sort(res.begin(),res.end());
  return res;
}

std::ostream& NC::operator<<( std::ostream& os, const NCGM::AtomicComponentList& cl )
{
  bool first(true);
  for ( auto& e : cl ) {
    if (!first)
      os<<'+';
    first = false;
    nc_assert_always(e.second.isElement());
    os << fmt(e.first)<<'*'<<elementZToName(e.second.Z());
  }
  return os;

}
