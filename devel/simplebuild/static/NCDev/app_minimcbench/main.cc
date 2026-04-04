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

#include "NCMiniMC/NCMMC_Query.hh"
#include "NCMiniMC/NCMMC_Utils.hh"
#include "NCMiniMC/NCMMC_RunSim.hh"
#include "NCMiniMC/NCMMC_Geom.hh"
#include "NCMiniMC/NCMMC_Source.hh"
#include "NCMiniMC/NCMMC_StdTallies.hh"
#include "NCMiniMC/NCMMC_SimEngine.hh"
#include "NCMiniMC/NCMMC_EngineOpts.hh"
#include "NCThreads/NCFactThreads.hh"
#include "NCUtils/NCString.hh"
#include <iostream>
#include <chrono>

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

int the_main_fct( int argc, char ** argv ) {

  const NC::ThreadCount nthreads{ argc > 1
                                  ? (unsigned)NC::str2int( argv[1] )
                                  : 1 };
  const char * src_cfgstr = ( argc > 2
                              ? argv[2]
                              : "constant;wl=1.8;n=1e6" );

  const char * geom_cfgstr = ( argc > 3
                               ? argv[3]
                               : "sphere;r=0.0128" );//unit r is meter

  const char * mat_cfgstr = ( argc > 4
                               ? argv[4]
                               : "Be3N2_sg206_BerylliumNitride.ncmat" );

  //Note: Do NOT enable the threadpool if we need the mmcscaling script to be
  //able to figure out the overhead and the Amdahl's p value for just the MiniMC
  //part.

  if ( false ) {
    NC::FactoryThreadPool::enable( nthreads );
  }

  std::string enginecfg;
  {
    std::ostringstream ss;
    ss << "tally=mu;nthreads="<<nthreads.get();
    enginecfg = ss.str();
  }

  auto matdef = NCMMC::MatDef{ mat_cfgstr };
  auto geom = NCMMC::createGeometry(geom_cfgstr);
  auto src = NCMMC::createSource(src_cfgstr);
  auto eopts = NCMMC::parseEngineOpts( enginecfg );
  auto srcmd = src->metaData();
  auto tally = NC::makeSO<NCMMC::TallyStdHists>( eopts,
                                                 src->metaData(),
                                                 matdef.matTemp );
  nc_assert_always( nthreads.get() > 0 );

  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

  auto resmd = NCMMC::runSim_StdEngine( geom, src, tally, matdef, eopts );

  nc_assert_always( resmd.provided.weight == resmd.provided.count );
  const auto nsrcparticles = resmd.provided.count;
  nc_assert_always( resmd.miss.count == 0 );
  nc_assert_always( resmd.miss.weight == 0.0 );
  nc_assert_always( resmd.tallied.count > resmd.provided.count );//splitting
  nc_assert_always( resmd.tallied.weight < resmd.provided.weight );//attenuation

  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  const double dt = (std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()) * 1e-6;
  std::cout << "    Number of threads used:  "<<nthreads<<std::endl;
  std::cout << "Number of source particles:  "<<resmd.provided.count<<std::endl;
  const double count_tallied = resmd.tallied.count;
  std::cout << "Number of tallied particles: "<<count_tallied<<std::endl;
  std::cout << "ntally/nsrc factor:          "<<count_tallied/double(nsrcparticles)<<std::endl;
  std::cout << "Time spent in sim-loop (s):  "<<dt<<std::endl;
  auto nthreads_min1 = (nthreads.get()?nthreads.get():1);
  std::cout << "            Time * nthread (s): "<<dt * nthreads_min1<<std::endl;
  std::cout << "       src particle rate (MHz): "<<(nsrcparticles*1e-6/dt)<<std::endl;
  std::cout << "            ..per thread (MHz): "<<(nsrcparticles*1e-6/dt)/nthreads_min1<<std::endl;
  std::cout << "     tally particle rate (MHz): "<<(count_tallied*1e-6/dt)<<std::endl;
  std::cout << "            ..per thread (MHz): "<<(count_tallied*1e-6/dt)/nthreads_min1<<std::endl;

  return 0;
}

int main( int argc, char ** argv ) {
  const unsigned nrepeat_for_profiling = 1;
  for (unsigned i = 0; i < nrepeat_for_profiling; ++i) {
    int ec = the_main_fct(argc,argv);
    if (ec!=0)
      return ec;
  };
  return 0;
}
