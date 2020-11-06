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

#include "NCFactory_NXS.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCAtomUtils.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCNXSLib.hh"
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace NCrystal {
  void initNXS( nxs::NXS_UnitCell* uc,
                const char * nxs_file,
                double temperature_kelvin,
                unsigned maxhkl,
                bool fixpolyatom )
  {
    const char * old_SgError = nxs::SgError;
    nxs::SgError = 0;

    nxs::NXS_AtomInfo *atomInfoList;
    int numAtomInfos = nxs::nxs_readParameterFile( nxs_file, uc, &atomInfoList);

    if (numAtomInfos==NXS_ERROR_READINGFILE)
      NCRYSTAL_THROW2(FileNotFound,"Could not find and open input file \""<<nxs_file<<"\"");

    if (numAtomInfos<=0)
      NCRYSTAL_THROW2(DataLoadError,
                     "Could not read crystal information from file \""<<nxs_file<<"\"");
    if( NXS_ERROR_OK != nxs::nxs_initUnitCell(uc) )
      NCRYSTAL_THROW2(DataLoadError,
                     "Could not initialise unit cell based on parameters in file \""<<nxs_file<<"\"");

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
                      "Could not initialise unit cell from file \""<<nxs_file
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

  struct XSectProvider_NXS final : public MoveOnly {
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

double NCrystal::XSectProvider_NXS::xsectScatNonBragg(const double& lambda) const
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

const NCrystal::Info * NCrystal::loadNXSCrystal( const char * nxs_file,
                                                 double temperature_kelvin,
                                                 double dcutoff_lower_aa,
                                                 double dcutoff_upper_aa,
                                                 bool bkgdlikemcstas,
                                                 bool fixpolyatom )
{
  const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
  if (verbose)
    std::cout<<"NCrystal::NCNXSFactory::invoked loadNXSCrystal("<< nxs_file
             <<", temp="<<temperature_kelvin
             <<", dcutoff="<<dcutoff_lower_aa
             <<", dcutoffup="<<dcutoff_upper_aa
             <<", bkgdlikemcstas="<<bkgdlikemcstas
             <<", fixpolyatom="<<fixpolyatom
             <<")"<<std::endl;

  ///////////////////////////////
  // Create CrystalInfo object //
  ///////////////////////////////

  Info * crystal = new Info();
  RCGuard guard(crystal);//prevent leaks in case of exceptions thrown

  //////////////////////////////////////////////////////////////////
  // Load and init NXS info (twice to figure out adequate maxhkl) //
  //////////////////////////////////////////////////////////////////

  struct NXSXSectProviderWrapper {
    //Dummy struct needed since std::function can only accept copy-able function
    //objects.
    std::shared_ptr<XSectProvider_NXS> shptr_xsprov_nxs;
    double operator()(double wl) const { nc_assert(!!shptr_xsprov_nxs); return shptr_xsprov_nxs->xsectScatNonBragg(wl); }
  };

  NXSXSectProviderWrapper xsect_provider{std::make_shared<XSectProvider_NXS>(bkgdlikemcstas)};
  nxs::NXS_UnitCell& nxs_uc = xsect_provider.shptr_xsprov_nxs->nxs_uc;

  const double fsquare_cut = 100.0 * 1e-5 ;//remove reflections with vanishing contribution
                                           //(NB: Hardcoded to same value as in .ncmat factory).
                                           //factor 100.0 is to convert to nxs units.

  int maxhkl = 1;
  if (verbose)
    std::cout<<"NCrystal::NCNXSFactory::calling nxslib initHKL with maxhkl="<<maxhkl<<" (to init non-hkl info)"<<std::endl;
  initNXS(&nxs_uc, nxs_file, temperature_kelvin, 1, fixpolyatom);

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
    initNXS(&nxs_uc, nxs_file, temperature_kelvin, maxhkl, fixpolyatom );
  }

  crystal->setXSectProvider(xsect_provider);

  //////////////////////
  // ... add HKL info //
  //////////////////////

  if (enable_hkl) {
    crystal->enableHKLInfo(dcutoff_lower_aa,dcutoff_upper_aa);
    nxs::NXS_HKL *it = &(nxs_uc.hklList[0]);
    nxs::NXS_HKL *itE = it + nxs_uc.nHKL;
    for (;it!=itE;++it) {
      if(it->dhkl < dcutoff_lower_aa || it->dhkl>dcutoff_upper_aa) //cut off d-spacing
        continue;
      if(it->FSquare < fsquare_cut) //remove reflections with vanishing contribution (left due to rounding errors?)
        continue;
      NCrystal::HKLInfo hi;
      hi.h = it->h;
      hi.k = it->k;
      hi.l = it->l;
      hi.multiplicity = it->multiplicity;
      hi.dspacing = it->dhkl;
      hi.fsquared = 0.01 * it->FSquare;
      crystal->addHKL(std::move(hi));
    }
    //We used to emit a warning here, but decided not to (user should be allowed
    //to deliberately exclude all bragg edges via the dcutoff parameter without
    //getting warnings):
    // if (!crystal->nHKL())
    //   printf("NCrystal::loadNXSCrystal WARNING: No HKL planes selected from file \"%s\"\n",nxs_file);
  }

  ////////////////////////////
  // ... add structure info //
  ////////////////////////////

  NCrystal::StructureInfo si;
  si.spacegroup = nxs_uc.sgInfo.TabSgName->SgNumber;
  si.lattice_a = nxs_uc.a;
  si.lattice_b = nxs_uc.b;
  si.lattice_c = nxs_uc.c;
  si.alpha = nxs_uc.alpha;
  si.beta  = nxs_uc.beta;
  si.gamma = nxs_uc.gamma;
  si.volume = nxs_uc.volume;
  si.n_atoms = nxs_uc.nAtoms;

  crystal->setStructInfo(si);

  /////////////////////////////////////////////
  // ... add cross section info and atom info//
  /////////////////////////////////////////////

  double sigma_free = 0;
  double sigma_abs = 0;
  unsigned ntot(0);

  std::vector<AtomInfo> aivec;
  aivec.reserve(nxs_uc.nAtomInfo);//correct unless multiple entries for same element
  double average_atomic_mass_amu(0.0);

  std::map<unsigned,AtomInfo> zval_2_atominfo;

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
                          nxsatom_capturexs, nxsatom_massamu, Zval );

    AtomInfo& ai = zval_2_atominfo[Zval];
    if (!ai.atom.atomDataSP) {
      //first time we saw this Z-value:
      ai.atom.atomDataSP = std::make_shared<const AtomData>(std::move(newatomdata));
      ai.atom.index.value = zval_2_atominfo.size()-1;
    } else {
      //We saw this Z-value before, demand that it was specified with identical parameters;
      if ( ! ai.atom.atomDataSP->sameValuesAs(newatomdata,1e-15,1e-15) ) {
        //NB: NCrystal can in principle support this, but not really sure how consistently NXS handles that internally...
        NCRYSTAL_THROW2(DataLoadError,"Inconsistent parameters specified in "<<nxs_file
                        <<" (multiple entries for same element have different values of physical constants)");
      }
    }

    nc_assert(ai.atom.atomDataSP!=nullptr);
    const unsigned multiplicity = nxs_ai.nAtoms;

    ai.number_per_unit_cell += multiplicity;
    sigma_abs += ai.atom.atomDataSP->captureXS() * multiplicity;
    sigma_free += ai.atom.atomDataSP->freeScatteringXS().val  * multiplicity;
    average_atomic_mass_amu += ai.atom.atomDataSP->averageMassAMU() * multiplicity;
    ntot += multiplicity;

    for (unsigned ipos = 0; ipos <multiplicity;++ipos)
      ai.positions.emplace_back(nxs_ai.x[ipos],nxs_ai.y[ipos],nxs_ai.z[ipos]);

  }

  if (ntot) {
    average_atomic_mass_amu /= ntot;
    sigma_abs /= ntot;
    sigma_free /= ntot;
    crystal->setXSectAbsorption(sigma_abs);
    crystal->setXSectFree(sigma_free);
  }

  for (auto& e : zval_2_atominfo)
    crystal->addAtom(std::move(e.second));

  //////////////////////////////
  // ... add temperature info //
  //////////////////////////////

  if (nxs_uc.temperature>0.0)
    crystal->setTemperature(nxs_uc.temperature);
  if (nxs_uc.debyeTemp>0.0)
    crystal->setGlobalDebyeTemperature(nxs_uc.debyeTemp);

  /////////////////////
  // ... add density //
  /////////////////////

  crystal->setDensity(nxs_uc.density);

  //1e27 in next line converts kg/Aa^3 to g/cm^3:
  const double numberdensity_2_density = 1e27 * average_atomic_mass_amu * constant_dalton2kg;
  crystal->setNumberDensity( nxs_uc.density / numberdensity_2_density );


  ///////////
  // Done! //
  ///////////

  //A bit dirty, but we don't really need to keep the hkllist around from this
  //point on, so might as well save some memory:
  deinitNXS_partly(&nxs_uc);

  crystal->objectDone();

  crystal->ref();
  guard.clear();
  crystal->unrefNoDelete();
  return crystal;
}
