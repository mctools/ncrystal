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
#include "NCrystal/internal/gasmix/NCGasMixUtils.hh"
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

    static constexpr auto gasmixquickfactname = "gasmix";
    std::string generateGasMixNCMAT( std::ostream& os,
                                     const std::string& str_request )
    {
      constexpr static int ncmat_version = AtomDBExtender::latest_version;//latest NCMAT version

      auto req = GasMix::requestFromString( str_request );
      auto str_request_normalised = requestToString(req);
      auto gm = GasMix::analyseGasMixRequest(req);

      os << "NCMAT v"<<ncmat_version<<"\n#\n# Automatically generated NCMAT data for gas mixture\n";
      os << "#\n# Request: \""<<gasmixquickfactname<<"::"<<str_request_normalised<<"\"\n#\n";
      os << "# Resulting gas parameters:\n#\n";
      os << "#  T   = "<<fmt(gm.temperature.dbl())<<"K\n";
      os << "#  P   = "<<fmt(gm.pressure.dbl())<<"Pa\n";
      os << "#  Rho = "<<fmt(gm.density.dbl()*1000.0)<<"kg/m^3\n#\n";
      os << "@STATEOFMATTER\n  gas\n";
      os << "@TEMPERATURE\n  "<<fmt(gm.temperature.dbl())<<"\n";
      os << "@DENSITY\n  "<< fmt(1000.0*gm.density.dbl()) << " kg_per_m3\n";
      if ( !req.atomDBLines.empty() ) {
        os << "@ATOMDB\n";
        for ( auto& e : req.atomDBLines )
          os << "  "<<joinstr(e)<<"\n";
      }
      auto flat_comp = GasMix::flattenComponentList( gm );
      for ( auto& comp : flat_comp ) {
        nc_assert_always(comp.second.isElement());
        os << "@DYNINFO\n  element  "
           <<elementZToName(comp.second.Z())
           <<"\n  fraction " << fmt_frac(comp.first );
        os <<"\n  type     freegas\n";
      }
      return str_request_normalised;
    }

    class QuickFactGasMix final : public NC::FactImpl::TextDataFactory
    {
    public:
      const char * name() const noexcept override
      {
        return gasmixquickfactname;
      }
      Priority query( const TextDataPath& ) const override
      {
        return Priority::OnlyOnExplicitRequest;
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        std::ostringstream ncmat_data;
        auto str_request_normalised = generateGasMixNCMAT( ncmat_data, p.path() );
        std::ostringstream src_descr;
        src_descr << "<automatically generated content from \""<<str_request_normalised<<"\">";
        return TextDataSource::createFromInMemData( RawStrData(ncmat_data.str(),src_descr.str().c_str()), "ncmat",
                                                    std::move(str_request_normalised) );
      }
      std::vector<BrowseEntry> browse() const override {
        //Can't really browse ALL possibilities of course, but we can instead
        //provide a few examples.
        std::vector<BrowseEntry> b;
        Priority p{Priority::OnlyOnExplicitRequest};
        const char * descr = "examples of on-demand gas mixtures";
        const char* examples[] = {
          "CO2",
          "He/10bar",
          "He/1.64kgm3",
          "0.7xCO2+0.3xAr/1.5atm/250K",
          "0.72xCO2+0.28xAr/massfractions/1.5atm/250K",
          "BF3/2atm/25C/B_is_0.95_B10_0.05_B11",
          "air",
          "air/-10C/0.8atm/0.30relhumidity",
          "0.7xCO2+0.3xAr/0.001relhumidity",
        };

        auto nExamples = sizeof(examples)/sizeof(*examples);
        b.reserve(nExamples);
        for ( auto i : ncrange(nExamples) )
          b.push_back( {examples[i], descr,p} );
        return b;
      }
    };
  }
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_quickgasmix_factory)()
{
  //Called from ncrystal_register_quick_factory(), to be part of the same
  //"stdquick" plugin.
  NC::FactImpl::registerFactory( std::make_unique<NC::QuickFactGasMix>() );
}
