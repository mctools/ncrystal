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

#include "NCFactory_NXS.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCAtomUtils.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCNXSLib.hh"
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace NC = NCrystal;

namespace NCrystal {
  struct FakeFileForNXSLoader {
    //Trying to mimic the behaviour of the fgets C-function, so we can load in-memory .nxs files
    FakeFileForNXSLoader( const TextData& td ) : m_it(td.begin()), m_itE(td.end()), m_dataDescr(td.description()) {}
    TextData::Iterator m_it;
    TextData::Iterator m_itE;
    std::string m_dataDescr;
  };
  //Need to pass function pointer to C doe. Thus we use a global function with
  //static data, and protect usage with a mutex:
  static FakeFileForNXSLoader * s_fakeFileForNXSLoader = nullptr;
  static std::mutex s_fakeFileForNXSLoader_mtx;
  char* fake_fgets_for_nxs( char* str,int count )
  {
    nc_assert(s_fakeFileForNXSLoader);
    auto& FF = *s_fakeFileForNXSLoader;
    nc_assert(count>=3);
    nc_assert(count<=1000000000);
    if ( FF.m_it == FF.m_itE ) {
      str[0] = '\0';
      return nullptr;
    }
    const std::string& nextLine = *FF.m_it;
    //std::fgets expected behaviour:
    //   Reads at most count - 1 characters from the given file stream and
    //   stores them in the character array pointed to by str. Parsing stops if
    //   a newline character is found, in which case str will contain that
    //   newline character, or if end-of-file occurs. If bytes are read and no
    //   errors occur, writes a null character at the position immediately after
    //   the last character written to str.
    //
    //Note that nextLine does NOT contain a newline character.
    const std::size_t n = FF.m_it->size();
    if  ( n+2 >= static_cast<std::size_t>(count) )
      NCRYSTAL_THROW2(DataLoadError,"Too long line encountered in .nxs input data: "<<FF.m_dataDescr);
    std::memcpy( str, &nextLine[0], n );
    str[n]='\n';
    str[n+1]='\0';
    ++FF.m_it;//advance to next line
    return str;
  }

  void initNXS( nxs::NXS_UnitCell* uc,
                const TextData& textData,
                //const std::string& nxs_file,
                double temperature_kelvin,
                unsigned maxhkl,
                bool fixpolyatom )
  {
    const char * old_SgError = nxs::SgError;
    nxs::SgError = 0;

    nxs::NXS_AtomInfo *atomInfoList;
    std::string dataDescr = textData.description();
    int numAtomInfos = 0;
    {
      FakeFileForNXSLoader fake( textData );
      NCRYSTAL_LOCK_GUARD(s_fakeFileForNXSLoader_mtx);
      s_fakeFileForNXSLoader = &fake;
      numAtomInfos = nxs::nxs_readParameterFile( fake_fgets_for_nxs, uc, &atomInfoList);
      s_fakeFileForNXSLoader = nullptr;
    }

    // Alternatively, we could work with on-disk files like this:
    // textData.verifyOnDiskFileUnchanged();
    // int numAtomInfos = nxs::nxs_readParameterFile( nxs_file.c_str(), uc, &atomInfoList);
    // textData.verifyOnDiskFileUnchanged();

    if (numAtomInfos==NXS_ERROR_READINGFILE) {
      NCRYSTAL_THROW2(DataLoadError,"Problems parsing input data: "<<dataDescr);
      //Was (but now we provide the content in-memory, so actual errors are pure syntax errors now):
      //NCRYSTAL_THROW2(FileNotFound,"Could not find and open input data: "<<dataDescr);
    }
    if (numAtomInfos<=0)
      NCRYSTAL_THROW2(DataLoadError,
                     "Could not read crystal information from data: "<<dataDescr);
    if( NXS_ERROR_OK != nxs::nxs_initUnitCell(uc) )
      NCRYSTAL_THROW2(DataLoadError,
                     "Could not initialise unit cell based on parameters in data: "<<dataDescr);

    uc->temperature = temperature_kelvin;
    for( int i=0; i< numAtomInfos; i++ )
      nxs::nxs_addAtomInfo( uc, atomInfoList[i] );
    free(atomInfoList);
    atomInfoList = 0;
    uc->maxHKL_index = maxhkl;
    int fix_incoh_xs = ( fixpolyatom ? 1 : 0 );
    nxs::nxs_initHKL( uc, fix_incoh_xs );
    if (nxs::SgError) {
      nxs::SgError = old_SgError;
      NCRYSTAL_THROW2(DataLoadError,
                      "Could not initialise unit cell from data \""<<dataDescr
                      <<"\" due to NXS errors: \""<<nxs::SgError<<"\"");
    }
    nxs::SgError = old_SgError;

    //Validate input, and replace b=0 and c=0 entries with a, for crystals where
    //it is allowed (this is needed since .nxs files allow b=0 or c=0 as
    //shorthand in cases where those are equal to a).
    checkAndCompleteLattice( uc->sgInfo.TabSgName->SgNumber, uc->a, uc->b, uc->c );

  }
  void deinitNXS_partly(nxs::NXS_UnitCell*uc)
  {
    if (!uc->hklList)
      return;
    nxs::NXS_HKL *it = &(uc->hklList[0]);
    nxs::NXS_HKL *itE = it + uc->nHKL;
    for (;it!=itE;++it)
      free(it->equivHKL);
    free(uc->hklList);
    uc->hklList = 0;
    free(uc->sgInfo.ListSeitzMx);
    uc->sgInfo.ListSeitzMx = 0;
  }
  void deinitNXS(nxs::NXS_UnitCell*uc)
  {
    deinitNXS_partly(uc);
    free(uc->atomInfoList);
    uc->atomInfoList = 0;
  }

  struct XSectProvider_NXS final : private MoveOnly {
    XSectProvider_NXS(bool bkgdlikemcstas)
      : m_bkgdlikemcstas(bkgdlikemcstas)
    {
      std::memset(&nxs_uc,0,sizeof(nxs_uc));
    }
    double xsectScatNonBragg(const double& lambda) const;
    ~XSectProvider_NXS()
    {
      deinitNXS(&nxs_uc);
    }
    nxs::NXS_UnitCell nxs_uc;
  private:
    bool m_bkgdlikemcstas;
  };
}

double NC::XSectProvider_NXS::xsectScatNonBragg(const double& lambda) const
{

  nxs::NXS_UnitCell* ucpar = const_cast<nxs::NXS_UnitCell*>(&nxs_uc);
  double xsect_cell;
  if ( m_bkgdlikemcstas )
    xsect_cell = ( nxs::nxs_IncoherentElastic(lambda, ucpar )
                   + nxs::nxs_IncoherentInelastic(lambda, ucpar )
                   + nxs::nxs_CoherentInelastic(lambda, ucpar ) );
  else
    xsect_cell = ( nxs::nxs_SinglePhonon(lambda, ucpar )
                   + nxs::nxs_MultiPhonon_COMBINED(lambda, ucpar )
                   + nxs::nxs_IncoherentElastic(lambda, ucpar ) );
  return xsect_cell > 0.0 ? xsect_cell / nxs_uc.nAtoms : 0.0;//protect against negative numbers and NaNs propagating from nxslib code.
}

NC::Info NC::loadNXSCrystal( const TextData& textData,
                             Temperature temperature,
                             double dcutoff_lower_aa,
                             double dcutoff_upper_aa,
                             bool bkgdlikemcstas,
                             bool fixpolyatom )
{
  std::string dataDescr = textData.description();

  const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
  if (verbose)
    std::cout<<"NCrystal::NCNXSFactory::invoked loadNXSCrystal("<< dataDescr
             <<", temp="<<temperature
             <<", dcutoff="<<dcutoff_lower_aa
             <<", dcutoffup="<<dcutoff_upper_aa
             <<", bkgdlikemcstas="<<bkgdlikemcstas
             <<", fixpolyatom="<<fixpolyatom
             <<")"<<std::endl;

  ///////////////////////////////
  // Create Info object //
  ///////////////////////////////

  Info info;

  //////////////////////////////////////////////////////////////////
  // Load and init NXS info (twice to figure out adequate maxhkl) //
  //////////////////////////////////////////////////////////////////

  struct NXSXSectProviderWrapper {
    //Dummy struct needed since std::function can only accept copy-able function
    //objects.
    std::shared_ptr<XSectProvider_NXS> shptr_xsprov_nxs;
    CrossSect operator()(NeutronEnergy ekin) const { nc_assert(!!shptr_xsprov_nxs); return CrossSect{ shptr_xsprov_nxs->xsectScatNonBragg(ekin.wavelength().dbl()) }; }
  };

  NXSXSectProviderWrapper xsect_provider{std::make_shared<XSectProvider_NXS>(bkgdlikemcstas)};
  nxs::NXS_UnitCell& nxs_uc = xsect_provider.shptr_xsprov_nxs->nxs_uc;

  const double fsquare_cut = 100.0 * 1e-5 ;//remove reflections with vanishing contribution
                                           //(NB: Hardcoded to same value as in .ncmat factory).
                                           //factor 100.0 is to convert to nxs units.

  int maxhkl = 1;
  if (verbose)
    std::cout<<"NCrystal::NCNXSFactory::calling nxslib initHKL with maxhkl="<<maxhkl<<" (to init non-hkl info)"<<std::endl;
  initNXS(&nxs_uc, textData, temperature.get(), maxhkl, fixpolyatom);

  const bool enable_hkl(dcutoff_lower_aa!=-1);
  if (enable_hkl) {
    RotMatrix rec_lat = getReciprocalLatticeRot( nxs_uc.a, nxs_uc.b, nxs_uc.c,
                                                 nxs_uc.alpha * kDeg, nxs_uc.beta * kDeg, nxs_uc.gamma * kDeg );

    if (dcutoff_lower_aa==0) {
      //have to determine appropriate dcutoff for this crystal. Aim for
      //maxhkl=20, but not outside [0.1,0.5].
      dcutoff_lower_aa = ncmax( 0.1, ncmin( 0.5, estimateDCutoff( 20, rec_lat ) ));
      std::string cmt;
      if (dcutoff_lower_aa>=dcutoff_upper_aa) {
        cmt = " (lower than usual due to value of dcutoffup)";
        dcutoff_lower_aa = 0.8*dcutoff_upper_aa;
      }
      if (verbose)
        std::cout<<"NCrystal::NCNXSFactory::automatically selected dcutoff level "<< dcutoff_lower_aa << " Aa"<<cmt<<std::endl;
    }

    int max_h, max_k, max_l;
    estimateHKLRange( dcutoff_lower_aa, rec_lat, max_h, max_k, max_l );

    maxhkl = ncmax(max_h,ncmax(max_k,max_l)) + 1;//+1 for safety
    deinitNXS(&nxs_uc);
    if (maxhkl>50)
      NCRYSTAL_THROW2(CalcError,"Combinatorics too great to reach requested dcutoff = "<<dcutoff_lower_aa<<" Aa");

    if (verbose)
      std::cout<<"NCrystal::NCNXSFactory::calling nxslib initHKL with maxhkl="<<maxhkl
               <<" (for all info including hkl)"<<std::endl;
    initNXS(&nxs_uc, textData, temperature.get(), maxhkl, fixpolyatom );
  }

  info.setXSectProvider(xsect_provider);

  //////////////////////
  // ... add HKL info //
  //////////////////////

  if (enable_hkl) {
    info.enableHKLInfo(dcutoff_lower_aa,dcutoff_upper_aa);
    nxs::NXS_HKL *it = &(nxs_uc.hklList[0]);
    nxs::NXS_HKL *itE = it + nxs_uc.nHKL;
    for (;it!=itE;++it) {
      if(it->dhkl < dcutoff_lower_aa || it->dhkl>dcutoff_upper_aa) //cut off d-spacing
        continue;
      if(it->FSquare < fsquare_cut) //remove reflections with vanishing contribution (left due to rounding errors?)
        continue;
      HKLInfo hi;
      hi.h = it->h;
      hi.k = it->k;
      hi.l = it->l;
      hi.multiplicity = it->multiplicity;
      hi.dspacing = it->dhkl;
      hi.fsquared = 0.01 * it->FSquare;
      info.addHKL(std::move(hi));
    }
    //We used to emit a warning here, but decided not to (user should be allowed
    //to deliberately exclude all bragg edges via the dcutoff parameter without
    //getting warnings):
    // if (!info.nHKL())
    //   printf("NCrystal::loadNXSCrystal WARNING: No HKL planes selected from file \"%s\"\n",nxs_file);
  }

  ////////////////////////////
  // ... add structure info //
  ////////////////////////////

  StructureInfo si;
  si.spacegroup = nxs_uc.sgInfo.TabSgName->SgNumber;
  si.lattice_a = nxs_uc.a;
  si.lattice_b = nxs_uc.b;
  si.lattice_c = nxs_uc.c;
  si.alpha = nxs_uc.alpha;
  si.beta  = nxs_uc.beta;
  si.gamma = nxs_uc.gamma;
  si.volume = nxs_uc.volume;
  si.n_atoms = nxs_uc.nAtoms;

  info.setStructInfo(si);

  /////////////////////////////////////////////
  // ... add cross section info and atom info//
  /////////////////////////////////////////////

  double sigma_free = 0;
  double sigma_abs = 0;
  unsigned ntot(0);

  std::vector<AtomInfo> aivec;
  aivec.reserve(nxs_uc.nAtomInfo);//correct unless multiple entries for same element
  double average_atomic_mass_amu(0.0);

  struct AtomInfoCollector {
    //Helper class collecting info needed for AtomInfo constructor
    OptionalAtomDataSP atomDataSP;
    Optional<AtomIndex> index;
    AtomInfo::AtomPositions positions;
    AtomInfo moveToAtomInfo( double raw_debyetemp_value )
    {
      nc_assert_always( index.has_value() );
      nc_assert_always( atomDataSP != nullptr );
      nc_assert_always( !positions.empty() );
      Optional<DebyeTemperature> dt;
      if ( raw_debyetemp_value > 0.0 )
        dt = DebyeTemperature{ raw_debyetemp_value };
      return AtomInfo( IndexedAtomData{ std::move(atomDataSP),index.value() },
                       std::move(positions), dt, NullOpt );
    }

  };

  std::map<unsigned,AtomInfoCollector> zval_2_atominfo;

  for( unsigned i=0; i<nxs_uc.nAtomInfo; i++ ) {
    nxs::NXS_AtomInfo & nxs_ai = nxs_uc.atomInfoList[i];
    if (!nxs_ai.nAtoms)
      NCRYSTAL_THROW(BadInput,"Encountered NXS atom data with unit cell count of 0");

    //We can not just get the atomic number from
    //nxs_ai.elementNumber, since it is never filled by the NXS
    //code, so we use NCrystal's elementNameToZ:
    const std::string label(nxs_ai.label);
    unsigned Zval = elementNameToZ(label);
    if (Zval==0)
      NCRYSTAL_THROW2(BadInput,"unknown or unsupported element symbol specified in .nxs file: \""<<label
                      <<"\" (should be written like H, He, Li, ...)");
    const double nxsatom_capturexs = nxs_ai.sigmaAbsorption;
    const double nxsatom_cohscatlen = nxs_ai.b_coherent*0.1;//unit fm-> 10fm=sqrt(barn) [c.f. note about coherentScatLen unit in NCAtomData.hh]
    const double nxsatom_incxs = nxs_ai.sigmaIncoherent;
    const double nxsatom_massamu = nxs_ai.M_m * const_neutron_atomic_mass;
    AtomData newatomdata( SigmaBound{nxsatom_incxs}, nxsatom_cohscatlen,
                          nxsatom_capturexs, AtomMass{nxsatom_massamu}, Zval );

    auto& ai = zval_2_atominfo[Zval];
    if ( ai.atomDataSP == nullptr ) {
      //first time we saw this Z-value:
      ai.atomDataSP = std::make_shared<const AtomData>(std::move(newatomdata));
      nc_assert_always( zval_2_atominfo.size() > 0 );
      nc_assert_always( zval_2_atominfo.size() < std::numeric_limits<AtomIndex::value_type>::max() );
      ai.index = AtomIndex{ static_cast<AtomIndex::value_type>( zval_2_atominfo.size()-1 ) };
    } else {
      //We saw this Z-value before, demand that it was specified with identical parameters;
      if ( ! ai.atomDataSP->sameValuesAs(newatomdata,1e-15,1e-15) ) {
        //NB: NCrystal can in principle support this, but not really sure how consistently NXS handles that internally...
        NCRYSTAL_THROW2(DataLoadError,"Inconsistent parameters specified in data\""<<dataDescr
                        <<" (multiple entries for same element have different values of physical constants)");
      }
    }

    nc_assert( ai.atomDataSP != nullptr );
    const unsigned multiplicity = nxs_ai.nAtoms;

    sigma_abs += ai.atomDataSP->captureXS() * multiplicity;
    sigma_free += ai.atomDataSP->freeScatteringXS().get()  * multiplicity;
    average_atomic_mass_amu += ai.atomDataSP->averageMassAMU().dbl() * multiplicity;
    ntot += multiplicity;

    ai.positions.reserve(ai.positions.size()+multiplicity);
    for ( auto ipos : ncrange(multiplicity) )
      ai.positions.emplace_back(nxs_ai.x[ipos],nxs_ai.y[ipos],nxs_ai.z[ipos]);

  }

  if (ntot) {
    average_atomic_mass_amu /= ntot;
    sigma_abs /= ntot;
    sigma_free /= ntot;
    info.setXSectAbsorption( SigmaAbsorption{ sigma_abs } );
    info.setXSectFree( SigmaFree{ sigma_free } );
  }

  for (auto& e : zval_2_atominfo)
    info.addAtom( e.second.moveToAtomInfo( nxs_uc.debyeTemp ) );//same Debye temp for all atoms

  //////////////////////////////
  // ... add temperature info //
  //////////////////////////////

  if (nxs_uc.temperature>0.0)
    info.setTemperature(Temperature{nxs_uc.temperature});

  /////////////////////
  // ... add density //
  /////////////////////

  info.setDensity( Density{ nxs_uc.density } );

  //1e27 in next line converts kg/Aa^3 to g/cm^3:
  const double numberdensity_2_density = 1e27 * average_atomic_mass_amu * constant_dalton2kg;
  info.setNumberDensity( NumberDensity{ nxs_uc.density / numberdensity_2_density } );


  ///////////
  // Done! //
  ///////////

  //A bit dirty, but we don't really need to keep the hkllist around from this
  //point on, so might as well save some memory:
  deinitNXS_partly(&nxs_uc);

  info.objectDone();
  return info;
}
