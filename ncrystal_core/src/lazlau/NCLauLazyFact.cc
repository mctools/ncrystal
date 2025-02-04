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
#include "NCrystal/factories/NCDataSources.hh"
#include "NCLazy.hh"

namespace NC = NCrystal;

///////////////////////////////////////////////////////
//                                                   //
// A small test factory providing a simple .lau/.laz //
// reader (probably not fully functional).           //
//                                                   //
///////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace {

    class AltLauFact final : public FactImpl::InfoFactory {
    public:
      const char * name() const noexcept override { return "stdlaz"; }

      Priority query( const FactImpl::InfoRequest& request ) const override
      {
        std::string ext = request.getDataType();
        return (ext=="laz"||ext=="lau") ? Priority{100} : Priority::Unable;
      }

      InfoPtr produce( const FactImpl::InfoRequest& request ) const override
      {
        return Lazy::buildInfoFromLazyData(request);
      }
    };
  }
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdlaz_factory)()
{
  NC::FactImpl::registerFactory(std::make_unique<NC::AltLauFact>());
  NC::DataSources::addRecognisedFileExtensions("laz");
  NC::DataSources::addRecognisedFileExtensions("lau");
}
