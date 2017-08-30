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

#include "NCrystal/NCInfo.hh"
#include "NCNCMatLoader.hh"
#include "NCPhononDebye.hh"
#include "NCrystal/NCCalcBase.hh"
#include "NCrystal/NCException.hh"
#include "NCFile.hh"
#include "NCNeutronSCL.hh"
#include "NCReflections.hh"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace NCrystal {

  struct XSectProvider_NCMAT : public XSectProvider {
    //Energy range from 1e-5eV to ncmin(1,log(DBL_MAX)*kt*2)).  Outside this
    //interval, appropriate extrapolations are used (see below)
    XSectProvider_NCMAT(double kt)
      : XSectProvider(),m_en_inel(logspace(-5,log10(ncmin(1,log(std::numeric_limits<double>::max())*kt*2)), 200)),
        m_xs_inel(m_en_inel.size(),0.), m_saturated_xs(-1.0),
        m_k(0.), m_k2(0.0)
    {
    }
    virtual double xsectScatNonBragg(const double& lambda) const;
    virtual ~XSectProvider_NCMAT()
    {
    }

    //internal data and methods used just by factory code below for initialisation:
    void accumInelastic(const std::vector<double> & xs, double frac);
    void setSaturatedXSAndInit(double xs)
    {
      nc_assert_always(m_saturated_xs==-1.0&&xs>=0.0);
      m_saturated_xs = xs;
      double lambda_edge = ekin2wl(m_en_inel.back());
      double xs_edge = m_xs_inel.back();
      nc_assert(lambda_edge>0);
      //see the implementation of xsectScatNonBragg for an explanation and usage of these constants:
      m_k = (xs_edge-m_saturated_xs)/(lambda_edge*lambda_edge);
      m_k2 = m_xs_inel.front() / ekin2wl( m_en_inel.front() );
      size_t n = m_xs_inel.size();
      m_k3.reserve(n);
      m_k3.push_back(std::make_pair(0.,0.));
      for (size_t u = 1; u < n; ++u) {
        size_t l = u-1;
        double Eu = m_en_inel.at(u);
        double Xu = m_xs_inel.at(u);
        double El = m_en_inel.at(l);
        double Xl = m_xs_inel.at(l);
        double dd = (Xu-Xl)/(Eu-El);
        m_k3.push_back(std::make_pair(Xl-El*dd,dd));
      }
      std::vector<double>().swap(m_xs_inel);//release m_xs_inel memory since no longer needed.
    }
    std::vector<double>& getEnergyVector() { return m_en_inel; }
    std::vector<double>& getXSVector() { return m_xs_inel; }
  private:
    std::vector<double> m_en_inel, m_xs_inel;
    std::vector<std::pair<double,double> > m_k3;
    double m_saturated_xs;
    double m_k;
    double m_k2;
  };

}

double NCrystal::XSectProvider_NCMAT::xsectScatNonBragg(const double& lambda) const
{
  double kiEn = wl2ekin(lambda); // wl2ekin costs 1 branch + 1 division + 1 multiplication
  std::vector<double>::const_iterator upper = std::upper_bound (m_en_inel.begin(), m_en_inel.end(), kiEn);

  if(upper == m_en_inel.end())
  {
    //interpolate xs between m_en_inel.back() and 0 Angstrom. In this region,
    //the contribution from Bragg diffraction will tend to zero as lambda^2 -
    //this is assuming we are at wavelengths low enough that all significant
    //Bragg edges already contribute (which we should hopefully be if our phonon
    //expansion had enough terms). Since we should approach a flat total
    //scattering cross-section, the phonon background cross-section should
    //consequently increase as m_saturated_xs+k*lambda^2 for some (usually
    //negative) constant k, which we fix by requiring continuity at lambda_edge
    //(see the implementation of setSaturatedXSAndInit for how we calculate k):
    return m_saturated_xs + m_k * lambda * lambda;
  }
  else if (upper == m_en_inel.begin() )
  {
    //Extrapolate to low energies using single-phonon 1/v law (m_k2 fixed by
    //continuity condition):
    return m_k2 * lambda;
  }
  else
  {
    //inside calculated range = interpolate in bins
    const std::pair<double,double>& interp_consts = m_k3[upper-m_en_inel.begin()];
    return interp_consts.first + interp_consts.second * kiEn;
  }

}

void NCrystal::XSectProvider_NCMAT::accumInelastic(const std::vector<double>& xs, double frac)
{
  nc_assert_always(m_saturated_xs==-1.0);

  if(xs.size()!=m_en_inel.size())
    NCRYSTAL_THROW(CalcError,"XSectProvider_NCMAT::accumInelastic input vector has different size than the xs vector");

  for(unsigned i=0;i<xs.size();i++)
    m_xs_inel[i] += xs[i]*frac;
}



namespace NCrystal {

  const Info * loadNCMATCrystal( const char * ncmat_file, double temperature_kelvin,
                                 double dcutoff_lower_aa, double dcutoff_upper_aa,
                                 int nphonon, bool expandhkl)
{

  const bool verbose = (std::getenv("NCRYSTAL_DEBUGINFO") ? true : false);
  if (verbose)
    std::cout<<"NCrystal::NCMATFactory::invoked loadNCMATCrystal("<< ncmat_file
             <<", temp="<<temperature_kelvin
             <<", dcutoff="<<dcutoff_lower_aa
             <<", dcutoffupper="<<dcutoff_upper_aa
             <<", nphonon="<<nphonon
             <<", expandhkl="<<expandhkl<<")"<<std::endl;

  NCMatLoader loader(ncmat_file);
  nc_assert_always(loader.m_atomnum>0);

  ////////////////////////
  // Create Info object //
  ////////////////////////

  Info * crystal = new Info();
  NeutronSCL *nscl = NeutronSCL::instance();

  //1. structure info
  StructureInfo si;
  si.spacegroup = loader.getSpacegroupNum();
  si.alpha = loader.m_alpha;
  si.beta = loader.m_beta;
  si.gamma = loader.m_gamma;
  si.lattice_a = loader.m_a;
  si.lattice_b = loader.m_b;
  si.lattice_c = loader.m_c;
  double volume = loader.m_lattice.m_cell.colX().cross(loader.m_lattice.m_cell.colY()).dot(loader.m_lattice.m_cell.colZ());
  si.volume = volume;
  si.n_atoms = loader.m_atomnum;
  crystal->setStructInfo(si);

  //2. temperature info
  if (temperature_kelvin>0.0)
    crystal->setTemperature(temperature_kelvin);
  if (loader.m_debye_temp>0.0)
    crystal->setDebyeTemperature(loader.m_debye_temp);

  //3. cross sections and atom info
  double sigma_abs=0;
  XSectProvider_NCMAT * xsect_provider = 0;

  std::map<std::string, double> msd_map;
  std::vector<double> xs;

  double cellsPerGram = 1e24/volume;//1e24 is the conversion factor from Aa**3 to cm**3
  double density =0.;
  double totalFreeScaXS = 0;

  std::map<std::string, double>& debye_map = loader.m_debye_map;

  for(std::map<std::string, std::vector<Vector> >::const_iterator it = loader.m_lattice.m_atomic_pos.begin();
      it != loader.m_lattice.m_atomic_pos.end();
      ++it )
  {
    //if per-element debye temps exists, they are preferable and should be
    //available for all elements. Otherwise a global debye temp should be
    //available (asserting here as a LogicError, since the loader should already
    //have checked this):
    nc_assert_always( debye_map.empty() ? (loader.m_debye_temp>0.0) : (debye_map.find(it->first)!=debye_map.end()));

    //phonon inelastic (create with nphonon=1 even when input nphonon=-1, since we always need the MSD map).
    PhononDebye inel( ( debye_map.empty() ? loader.m_debye_temp*constant_boltzmann : debye_map.at(it->first)*constant_boltzmann ),
                      temperature_kelvin*constant_boltzmann, it->first, ( nphonon==-1 ? 1 : nphonon) );
    msd_map[it->first] = inel.getMSD();

    if(nphonon!=-1)
    {
      if (verbose&&nphonon==0)
        std::cout<<"NCrystal::NCMATFactory::automatically selected nphonon level "<< inel.getMaxPhononNum()<<std::endl;

      if (!xsect_provider)
        xsect_provider = new XSectProvider_NCMAT(temperature_kelvin*constant_boltzmann);

      inel.doit(xsect_provider->getEnergyVector(), xs);
      xsect_provider->accumInelastic(xs,double(it->second.size())/loader.m_atomnum);
    }


    msd_map[it->first] = inel.getMSD();

    //capture, bound incoherent and coherent
    unsigned atnum = it->second.size();
    sigma_abs += atnum* nscl->getCaptureXS(it->first);

    //free sacttering xs
    totalFreeScaXS += atnum* nscl->getFreeXS(it->first);
    //density
    density += it->second.size()*nscl->getNeutronWeightedMass(it->first)*const_neutron_mass;

    //atominfo
    AtomInfo ai;
    ai.atomic_number = nscl->getAtomicNumber(it->first);
    ai.number_per_unit_cell = atnum;
    if (!loader.m_debye_map.empty()) {
      nc_assert_always(loader.m_debye_map.find(it->first)!=loader.m_debye_map.end());//should have been checked by loader
      ai.debye_temp = loader.m_debye_map[it->first];
    } else {
      ai.debye_temp = 0;
    }
    crystal->addAtom(ai);
  }
  density*=cellsPerGram;

  sigma_abs /= loader.m_atomnum;
  double sigma_free = totalFreeScaXS/ loader.m_atomnum;
  crystal->setXSectAbsorption(sigma_abs);
  crystal->setXSectFree(sigma_free);

  if (xsect_provider) {
    //Shorten xsect_provider's energy/cross-section vectors to the point where
    //the phonon spectrum drops off due to insufficient terms in the expansion
    //(querying beyond the vector triggers extrapolation towards the
    //saturated/free cross-section at 0Aa).:
    std::vector<double> & v_e  = xsect_provider->getEnergyVector();
    std::vector<double> & v_xs = xsect_provider->getXSVector();
    nc_assert(v_e.size()==v_xs.size()&&v_xs.size()>1);
    size_t n = v_e.size()-1;
    while (n) {
      double delta_xs = v_xs.at(n)-v_xs.at(n-1);
      double delta_wl = ekin2wl(v_e.at(n)) - ekin2wl(v_e.at(n-1));
      double dxsdwl = delta_xs/delta_wl;
      if ( dxsdwl <= 0 && n+1 == v_e.size() )
        break;//XS decreasing with wl at edge - don't perform truncation at all in this case.
      if (dxsdwl < 0.075*sigma_free)//XS no longer increasing rapidly, stop here ( 0.075 barn/Aa is a hand-tuned value).
        break;
      --n;
    }
    if (!n)
      NCRYSTAL_THROW(CalcError,"Cross sections from phonon expansion keeps increasing rapidly over all wavelengths");
    //perform truncation:
    v_e.resize(n);
    v_xs.resize(n);
#if __cplusplus >= 201103L
    v_e.shrink_to_fit();
    v_xs.shrink_to_fit();
#endif
    if (verbose)
      std::cout<<"NCrystal::NCMATFactory::multi-phonon cross-section extrapolated below "
               <<ekin2wl(v_e.back()) <<" Aa"<<std::endl;

  }

  if (xsect_provider) {
    //finish up initialisation and register:
    xsect_provider->setSaturatedXSAndInit(sigma_free);
    crystal->setXSectProvider(xsect_provider);
  }

  //4. HKL info
  if ( dcutoff_lower_aa != -1 )
  {
    if(dcutoff_lower_aa==0) {
    //Very simple heuristics here for now, can likely be improved (specifically
    //we needed 0.4Aa for expensive Y2O3 with ~80 atoms/cell):
      if (loader.m_atomnum>40)
        dcutoff_lower_aa = 0.4;
      else
        dcutoff_lower_aa = 0.15;
      std::string cmt;
      if (dcutoff_lower_aa>=dcutoff_upper_aa) {
        cmt = " (lower than usual due to value of dcutoffupper)";
        dcutoff_lower_aa = 0.5*dcutoff_upper_aa;
      }
      if (verbose)
        std::cout<<"NCrystal::NCMATFactory::automatically selected dcutoff level "<< dcutoff_lower_aa << " Aa"<<cmt<<std::endl;
    }
    crystal->enableHKLInfo(dcutoff_lower_aa,dcutoff_upper_aa);
    loader.fillHKL(*crystal,  msd_map , dcutoff_lower_aa , dcutoff_upper_aa, expandhkl);
  }

  //5. density
  crystal->setDensity(density);


  ///////////
  // Done! //
  ///////////

  crystal->objectDone();
  return crystal;

}
}
///////////////////////////////////////////////////////////////////
//Finally, the code which can be used to enable our NCMAT loader //
///////////////////////////////////////////////////////////////////

#include "NCrystal/NCFactoryRegistry.hh"
#include "NCrystal/NCMatCfg.hh"

namespace NCrystal {

  class NCMATFactory : public FactoryBase {
  public:
    const char * getName() const { return "NCrystalNCMATFactory"; }

    virtual int canCreateInfo( const MatCfg& cfg ) const {
      return cfg.getDataFileExtension()=="ncmat" ? 100 : 0;
    }
    virtual const Info * createInfo( const MatCfg& cfg ) const
    {
      nc_assert_always(canCreateInfo(cfg));
      return loadNCMATCrystal( cfg.getDataFile().c_str(), cfg.get_temp(),
                               cfg.get_dcutoff(), cfg.get_dcutoffupper(),
                               cfg.get_nphonon(), cfg.get_expandhkl() );
    }
  };

}


//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_ncmat_factory()
{
  if (!NCrystal::hasFactory("NCrystalNCMATFactory"))
    registerFactory(new NCrystal::NCMATFactory);
}
