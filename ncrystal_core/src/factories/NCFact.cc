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

#include "NCrystal/factories/NCFact.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCDataSources.hh"
#include "NCrystal/interfaces/NCRNG.hh"

namespace NC = NCrystal;

NC::Scatter NC::createScatter( const MatCfg& cfg )
{
  auto rngproducer = getDefaultRNGProducer();
  auto rng = rngproducer->produce();
  return Scatter( std::move(rngproducer),
                  std::move(rng),
                  FactImpl::createScatter( cfg ) );
}

NC::Scatter NC::createScatter_RNGByIdx( const MatCfg& cfg, RNGStreamIndex rngidx )
{
  auto rngproducer = getDefaultRNGProducer();
  auto rng = rngproducer->produceByIdx(rngidx);
  return Scatter( std::move(rngproducer),
                  std::move(rng),
                  FactImpl::createScatter( cfg ) );
}

NC::Scatter NC::createScatter_RNGForCurrentThread( const MatCfg& cfg )
{
  auto rngproducer = getDefaultRNGProducer();
  auto rng = rngproducer->produceForCurrentThread();
  return Scatter( std::move(rngproducer),
                  std::move(rng),
                  FactImpl::createScatter( cfg ) );
}

NC::Absorption NC::createAbsorption( const MatCfg& cfg )
{
  return Absorption( FactImpl::createAbsorption( cfg ) );
}

NC::shared_obj<const NC::Info> NC::createInfo( const MatCfg& cfg )
{
  return FactImpl::createInfo(cfg);
}

void NC::registerInMemoryFileData( std::string virtualFileName,
                                    std::string&& data )
{
  DataSources::registerInMemoryFileData( std::move(virtualFileName),
                                         std::move(data) );
}

void NC::registerInMemoryStaticFileData( std::string virtualFileName,
                                          const char* static_data )
{
  DataSources::registerInMemoryStaticFileData( std::move(virtualFileName),
                                               static_data );
}
