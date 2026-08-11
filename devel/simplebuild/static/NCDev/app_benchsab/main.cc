////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include <iostream>


namespace NC = NCrystal;

void doBenchmark( int knllux, NC::shared_obj<const NC::SABData> sab )
{
  auto cfg = NC::SABCfg::createConfig(knllux);
  using SABProcessor = NC::SABUtils::SABProcessor;
  const auto do_sample = SABProcessor::SampleSupport::YES;
  const auto do_diagnostics = SABProcessor::StoreExtraDiagnostics::NO;
  std::shared_ptr<const NC::VectD> egrid = nullptr;//fixme?

  SABProcessor( cfg, std::move(sab), egrid, do_sample, do_diagnostics );
}

int main( int argc, char** argv )
{
  //Parse args:
  if ( argc!=4 ) {
    std::cout << "Please provide args: [cfgstr] [displaylabel] [nload]"
              << std::endl;
    std::cout << "NB: displaylabel can be empty (\"\") for monoatomic materials"
              << std::endl;
    return 1;
  }
  const std::string cfgstr(argv[1]);
  NC::MatCfg matcfg(cfgstr);
  const auto vdoslux = matcfg.get_vdoslux();
  const auto knllux = matcfg.get_knllux();
  const std::string displaylbl(argv[2]);
  auto opt_nload = NC::StrView(argv[3]).toUInt64();
  nc_assert_always(opt_nload.has_value());
  const auto nload = opt_nload.value();

  //Find SABData:
  auto sab =[&cfgstr,&displaylbl,vdoslux]()
  {
    NC::Optional<std::string> dlbl;
    if ( !displaylbl.empty() )
      dlbl = displaylbl;
    auto info = NC::FactImpl::createInfo(cfgstr);
    auto di = NC::InfoUtils::findDynInfo( info, dlbl );
    auto di_knl = dynamic_cast<const NC::DI_ScatKnl*>(di);
    if (!di_knl)
      NCRYSTAL_THROW(BadInput,"Selected component does not provide a SAB knl");
    return NC::extractSABDataFromDynInfo( di_knl, vdoslux );
  }();

  //Initiate SABProcessor:
  for ( std::uint64_t i = 0; i < nload; ++i )
    doBenchmark( knllux, sab );

  return 0;
}
