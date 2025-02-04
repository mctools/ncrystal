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

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/atomdb/NCAtomDBExtender.hh"

// "Quick" TextData factories which rather than working from pre-existing file
// data, simply generates NCMAT data on the basis of information encoded directly
// in the request key itself. This information consists of material densities and
// simple chemical formulas for material composition. Two instances of the
// factory exist, "freegas" and "solid". The former will model materials with the
// free gas model, and the latter as amorphous solids assuming a Debye model
// phonon DOS and a Debye temperature (defaults to 300K).

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    enum class QuickType { FreeGas, Solid };
    std::string generateQuickNCMAT( std::ostream& os,
                                    QuickType qtype,
                                    const std::string& request )
    {
      constexpr static int ncmat_version = AtomDBExtender::latest_version;//latest NCMAT version

      //accepted units and their NCMAT version (NB: Put kgcm3 before gcm3 since kgcm3 ends with gcm3):
      std::array<PairSS,3> density_units
        = { PairSS{ "perAa3", "atoms_per_aa3"},
            PairSS{ "kgm3",  "kg_per_m3"},
            PairSS{ "gcm3",   "g_per_cm3"} };

      const char * errmsgprefix = ( qtype == QuickType::FreeGas
                                    ? "Syntax error in freegas specification: "
                                    : "Syntax error in solid specification: " );
      if (!isSimpleASCII(request,AllowTabs::No,AllowNewLine::No))
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Does not contain only simple ASCII characters :\""<<request<<"\"");
      if (contains_any(request,"\"'$@\\*<>"))
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Contains forbidden characters :\""<<request<<"\"");

      Optional<std::pair<std::string,std::string>> density_line;
      Optional<std::pair<std::string,DecodedChemForm>> chemform;
      std::map<AtomSymbol,std::string> debyeTemps;
      Optional<std::string> globalDebyeTemp;
      SmallVector<VectS,5> atomdb_lines;

      std::string request_normalised;
      request_normalised.reserve(request.size());
      for ( auto e : trimEntries(split2(request,0,'/')) ) {
        if ( e.empty() )
          continue;
        if ( !request_normalised.empty() )
          request_normalised += "/";
        request_normalised += e;
        //Debye temperatures (solids only):
        if ( startswith(e,"TDebye") ) {
          if ( qtype == QuickType::FreeGas ) {
            NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Debye temperature entries like \""
                            <<e<<"\" are only allowed for a solid.");
          }
          //Allow per-element or global debye temperatures like:
          //   "solid:Al2O3/3gcm/TDebye200K_Al"
          //   "solid:Al2O3/3gcm/TDebye200K"
          auto parts = trimEntries(split2(e.substr(6),0,'_'));
          if ( parts.empty() || parts.size() > 2 || !(endswith(parts.front(),'K'))
               || [&parts]()
               {
                 for ( auto& ep : parts )
                   if ( ep.empty() )
                     return true;
                 return false;
               }() ) {
            NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Invalid Debye temperature syntax: \""<<e<<"\".");
          }
          auto tmp_valstr = parts.front().substr(0,parts.front().size()-1);
          double tmp_val;
          if ( !safe_str2dbl( tmp_valstr, tmp_val ) || std::isinf(tmp_val) || std::isnan(tmp_val) || !(tmp_val>0.0) || !(tmp_val<1e9) )
            NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Invalid value in Debye temperature specification: \""<<tmp_valstr<<"\".");
          if ( parts.size() == 1 ) {
            if ( globalDebyeTemp.has_value() )
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"More than one global Debye temperature specified.");
            globalDebyeTemp = tmp_valstr;
          } else {
            nc_assert_always(parts.size()==2);
            AtomSymbol as(parts.back());
            if (!as.isElement())
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Not a valid natural element \""<<parts.back()<<"\"");
            auto test = debyeTemps.insert( decltype(debyeTemps)::value_type{ as, tmp_valstr } );
            if ( !test.second )
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"More than one Debye temperature for"
                              " element \""<<parts.back()<<"\" specified.");
          }
          continue;
        }
        //Density:
        bool handled = false;
        for ( auto& du : density_units ) {
          if ( endswith(e,du.first) ) {
            if ( density_line.has_value() )
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple density entries are not allowed.");
            double tmp_val;
            std::string tmps = e.substr( 0, e.size()-du.first.size() );
            if ( !safe_str2dbl(tmps,tmp_val) || std::isinf(tmp_val) || std::isnan(tmp_val) || !(tmp_val>0.0) || !(tmp_val<1e10) )
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Density value is invalid or out-of-range: \""<<e<<"\"");
            density_line.emplace( e, tmps + " " + du.second);
            handled = true;
            break;
          }
        }
        if ( handled )
          continue;
        //Chemical formula:
        auto tmp_chemform = tryDecodeSimpleChemicalFormula( e );
        if ( tmp_chemform.has_value() ) {
          if ( chemform.has_value() )
            NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Multiple chemical formula entries.");
          chemform.emplace( e, std::move(tmp_chemform.value()) );
          continue;
        }
        //AtomDB:
        {
          //NB: Code duplicated between here and NCGasMixUtils.cc!!
          std::string tmp(e);
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
            if ( contains(e,':') || contains(e,' ') )
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Do not use spaces or semicolons as separators in ATOMDB sections (use single underscores instead): "<<e<<"\".");
            if ( contains(e,"__") )
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Do not use double underscores as separators in ATOMDB sections: "<<e<<"\".");
            try {
              validateAtomDBLine(parts, ncmat_version );
            } catch ( Error::BadInput& err ) {
              NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Invalid ATOMDB syntax: "<<err.what()<<"\".");
            }
            atomdb_lines.push_back(std::move(parts));
            continue;
          }
        }

        //Nothing:
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Unrecognised syntax: \""<<e<<"\".");
      }

      if ( !chemform.has_value() )
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Missing chemical formula.");
      if ( !density_line.has_value() )
        NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Missing density specification.");

      if ( qtype == QuickType::Solid ) {
        //using DecodedChemForm = SmallVector<std::pair<std::uint_least32_t,AtomSymbol>,4>;
        for ( const auto& e : debyeTemps ) {
          bool found = false;
          for ( auto ec : chemform.value().second ) {
            if ( ec.second == e.first ) {
              found = true;
              break;
            }
          }
          if (!found)
            NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Debye temperature specified for invalid (or unused) atom.");
        }
        //User should not specify both global and per-element Debye temperatures:
        if ( !debyeTemps.empty() && globalDebyeTemp.has_value() )
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Do not specify both per-element and global Debye temperatures.");
        //Nonetheless, we (ab)use the globalDebyeTemp variable to hold the fall-back choice:
        if ( !globalDebyeTemp.has_value() )
          globalDebyeTemp = "300";
      }

      auto atomsymb2name = [](const AtomSymbol& as) { nc_assert(as.isElement()); return elementZToName(as.Z()); };
      os << "NCMAT v"<<ncmat_version<<"\n#\n# Automatically generated NCMAT data for "<<chemform.value().first
         <<" ("<<(qtype==QuickType::FreeGas?"free gas model":"amorphous solid, Debye model")<<")\n";
      os << "#\n# Request: \""<<(qtype == QuickType::FreeGas?"freegas::":"solid::")<<request_normalised<<"\"\n#\n";

      os << "@STATEOFMATTER\n  "
         <<(qtype==QuickType::FreeGas?"gas":"solid")<<"\n";
      os << "@DENSITY\n  "
         <<density_line.value().second<<"\n";

      if ( !atomdb_lines.empty() ) {
        os << "@ATOMDB\n";
        for ( auto& e : atomdb_lines)
          os << "  "<<joinstr(e)<<"\n";
      }

      std::uint64_t ntot = 0;
      for ( const auto& e : chemform.value().second )
        ntot += e.first;

      for ( const auto& e : chemform.value().second ) {
        if (!e.second.isElement())
          NCRYSTAL_THROW2(BadInput,errmsgprefix<<"Chemical formula \""<<chemform.value().first
                          <<"\" does not exclusively contain natural elements");
        os << "@DYNINFO\n  element  "<<atomsymb2name(e.second)
           <<"\n  fraction ";
        if ( e.first==ntot )
          os<<"1";
        else
          os<<e.first<<"/"<<ntot;
        if ( qtype == QuickType::FreeGas ) {
          os <<"\n  type     freegas\n";
        } else {
          auto it = debyeTemps.find(e.second);
          nc_assert(globalDebyeTemp.has_value());
          std::string dt_str = ( it == debyeTemps.end() ) ? globalDebyeTemp.value() : it->second;
          os <<"\n  type     vdosdebye\n  debye_temp "<<dt_str<<"\n";
        }
      }
      return request_normalised;
    }

    class QuickFact final : public NC::FactImpl::TextDataFactory
    {
      QuickType m_qtype;
    public:
      QuickFact( QuickType qtype ) : m_qtype(qtype) {}
      const char * name() const noexcept override
      {
        return m_qtype == QuickType::FreeGas ? "freegas" : "solid";
      }
      Priority query( const TextDataPath& ) const override
      {
        return Priority::OnlyOnExplicitRequest;
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        std::string request = p.path();
        std::ostringstream src_descr;
        src_descr << "<automatically generated content from \""<<p.toString()<<"\">";
        std::ostringstream ncmat_data;
        std::string request_normalised = generateQuickNCMAT( ncmat_data, m_qtype, p.path() );
        return TextDataSource::createFromInMemData( RawStrData(ncmat_data.str(),src_descr.str().c_str()), "ncmat",
                                                    std::move(request_normalised) );
      }
      std::vector<BrowseEntry> browse() const override {
        //Can't really browse ALL possibilities of course, but we can instead
        //provide a few examples.
        std::vector<BrowseEntry> b;
        Priority p{Priority::OnlyOnExplicitRequest};
        const char * descr = "examples of on-demand unstructured materials";
        if ( m_qtype == QuickType::FreeGas ) {
          b.reserve(4);
          b.push_back( {"CF4/3.72kgm3", descr,p} );
          b.push_back( {"CO2/1.98kgm3", descr,p} );
          b.push_back( {"He/0.17kgm3/He_is_He3", descr,p} );
          b.push_back( {"Ar/2.5e-5perAa3", descr,p} );
        } else {
          b.reserve(6);
          b.push_back( {"CH2/1gcm3", descr,p} );
          b.push_back( {"Gd2O3/7.07gcm3", descr,p} );
          b.push_back( {"B4C/2.52gcm3/B_is_0.95_B10_0.05_B11", descr,p} );
          b.push_back( {"Al2O3/4gcm3", descr,p} );
          b.push_back( {"Al2O3/4gcm3/TDebye900K", descr,p} );
          b.push_back( {"Al2O3/4gcm3/TDebye750K_Al/TDebye1000K_O", descr,p} );
        }
        return b;
      }
    };

  }
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_quickgasmix_factory)();

extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_quick_factory)()
{
  NC::FactImpl::registerFactory( std::make_unique<NC::QuickFact>(NC::QuickType::FreeGas) );
  NC::FactImpl::registerFactory( std::make_unique<NC::QuickFact>(NC::QuickType::Solid) );
  //the gasmix factory should be part of the same "stdquick" plugin.
  NCRYSTAL_APPLY_C_NAMESPACE(register_quickgasmix_factory)();
}
