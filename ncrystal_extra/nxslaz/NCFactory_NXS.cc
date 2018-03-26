////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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
#include "NCrystal/NCException.hh"
#include "NCMath.hh"
#include "NCLatticeUtils.hh"
#include "NCNXSLib.hh"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

namespace NCrystal {
  void initNXS(nxs::NXS_UnitCell* uc, const char * nxs_file, double temperature_kelvin, unsigned maxhkl,
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

  struct XSectProvider_NXS : public XSectProvider {
    XSectProvider_NXS(bool bkgdlikemcstas)
      : m_bkgdlikemcstas(bkgdlikemcstas)
    {
      std::memset(&nxs_uc,0,sizeof(nxs_uc));
    }
    virtual double xsectScatNonBragg(const double& lambda) const;
    virtual ~XSectProvider_NXS()
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

  XSectProvider_NXS * xsect_provider = new XSectProvider_NXS(bkgdlikemcstas);
  crystal->setXSectProvider(xsect_provider);//register immediately for proper memory cleanup in case of errors
  nxs::NXS_UnitCell& nxs_uc = xsect_provider->nxs_uc;

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
                                                 nxs_uc.alpha * M_PI/180, nxs_uc.beta * M_PI/180, nxs_uc.gamma * M_PI/180 );

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

  //////////////////////
  // ... add HKL info //
  //////////////////////

  if (enable_hkl) {
    crystal->enableHKLInfo(dcutoff_lower_aa,dcutoff_upper_aa);
    NCrystal::HKLInfo hi;
    nxs::NXS_HKL *it = &(nxs_uc.hklList[0]);
    nxs::NXS_HKL *itE = it + nxs_uc.nHKL;
    for (;it!=itE;++it) {
      if(it->dhkl < dcutoff_lower_aa || it->dhkl>dcutoff_upper_aa) //cut off d-spacing
        continue;
      if(it->FSquare < fsquare_cut) //remove reflections with vanishing contribution (left due to rounding errors?)
        continue;
      hi.h = it->h;
      hi.k = it->k;
      hi.l = it->l;
      hi.multiplicity = it->multiplicity;
      hi.dspacing = it->dhkl;
      hi.fsquared = 0.01 * it->FSquare;
      crystal->addHKL(hi);
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

  //TODO for NC2: Some duplication with NCNeutronSCL in the following table,
  //which should likely exist somewhere else in a light weight form like this
  //(can then also be used in G4NCMatHelper in place of GetNistElementNames()):
  static const char* s_elementsSymbols[] = {"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
                                            "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
                                            "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
                                            "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
                                            "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
                                            "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
                                            "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
                                            "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
                                            "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
                                            "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
                                            "Ds", "Rg"};

  double sigma_free = 0;
  double sigma_abs = 0;
  unsigned ntot(0);
  std::vector<AtomInfo> aivec;
  for( unsigned i=0; i<nxs_uc.nAtomInfo; i++ ) {
    nxs::NXS_AtomInfo & nxs_ai = nxs_uc.atomInfoList[i];
    if (!nxs_ai.nAtoms)
      NCRYSTAL_THROW(BadInput,"Encountered NXS atom data with unit cell count of 0");

    sigma_abs += nxs_ai.sigmaAbsorption * nxs_ai.nAtoms;

    double massceof = (1+nxs_ai.M_m)/nxs_ai.M_m;
    sigma_free += (nxs_ai.b_coherent*nxs_ai.b_coherent*0.01*4*M_PI + nxs_ai.sigmaIncoherent)/(massceof*massceof) * nxs_ai.nAtoms;
    ntot += nxs_ai.nAtoms;

    AtomInfo ai;
    ai.number_per_unit_cell = nxs_ai.nAtoms;
    ai.debye_temp = 0;//not available

    //We can not just get the atomic number from
    //nxs_ai.elementNumber, since it is never filled by the NXS
    //code. Instead we look in the table above by brute-force:
    ai.atomic_number = 0;
    std::string label(nxs_ai.label);
    for (unsigned ie = 0; ie < sizeof(s_elementsSymbols)/sizeof(char*); ++ie) {
      if ( label == s_elementsSymbols[ie] ) {
        ai.atomic_number = ie+1;
        break;
      }
    }
    if (!ai.atomic_number)
      NCRYSTAL_THROW2(BadInput,"unknown element symbol specified in .nxs file: \""<<label
                      <<"\" (should be written like H, He, Li, ...)");


    std::vector<AtomInfo>::iterator itai;
    for(itai=aivec.begin();itai!=aivec.end();++itai)
    {
      if(itai->atomic_number==ai.atomic_number)
      {
        itai->number_per_unit_cell+=ai.number_per_unit_cell;
        break;
      }
    }
    if(itai==aivec.end())
      aivec.push_back(ai);
  }

  for(std::vector<AtomInfo>::iterator itai=aivec.begin();itai!=aivec.end();++itai)
  {
    crystal->addAtom(*itai);
  }

  if (ntot) {
    sigma_abs /= ntot;
    sigma_free /= ntot;
    crystal->setXSectAbsorption(sigma_abs);
    crystal->setXSectFree(sigma_free);
  }

  //////////////////////////////
  // ... add temperature info //
  //////////////////////////////

  if (nxs_uc.temperature>0.0)
    crystal->setTemperature(nxs_uc.temperature);
  if (nxs_uc.debyeTemp>0.0)
    crystal->setDebyeTemperature(nxs_uc.debyeTemp);

  //////////////////////
  // ... and density //
  //////////////////////
  crystal->setDensity(nxs_uc.density);


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
