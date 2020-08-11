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

#include "NCrystal/ncrystal.h"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCAbsorption.hh"
#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/NCFactory.hh"
#include "NCrystal/NCFactoryRegistry.hh"
#include "NCDynInfoUtils.hh"
#include "NCrystal/NCDump.hh"
#include "NCMath.hh"
#include "NCNeutronSCL.hh"
#include <cstring>
#include <cstdio>
#include <cstdlib>

namespace NCrystal {

  namespace NCCInterface {

    static RCHolder<RandomBase> saved_rng;

    static int quietonerror = 0;
    static int haltonerror = 1;
    static int waserror = 0;
    static char errmsg[512];
    static char errtype[64];
    static void (*custom_error_handler)(char *,char*) = 0;

    void setError(const char *msg, const char * etype = 0) throw() {
      if (!etype)
        etype="ncrystal_c-interface";
      strncpy(errmsg,msg,sizeof(errmsg)-1);
      strncpy(errtype,etype,sizeof(errtype)-1);
      //Ensure final null-char in case of very long input strings:
      errmsg[sizeof(errmsg)-1]='\0';
      errtype[sizeof(errtype)-1]='\0';
      if (custom_error_handler) {
        (*custom_error_handler)(errtype,errmsg);
      }
      waserror = 1;
      if (!quietonerror)
        printf("NCrystal ERROR [%s]: %s\n",errtype,errmsg);
      if (haltonerror) {
        printf("NCrystal terminating due to ERROR\n");
        exit(1);
      }
    }

    void handleError(const std::exception &e) throw() {
      const Error::Exception* nce = dynamic_cast<const Error::Exception*>(&e);
      if (nce) {
        setError(nce->what(),nce->getTypeName());
        return;
      }
      const std::runtime_error* stdrte = dynamic_cast<const std::runtime_error*>(&e);
      if (stdrte)
        setError(stdrte->what(),"std::runtime_error");
      else
        setError("<unknown>","std::exception");
    }

    class RandFctWrapper : public RandomBase {
    public:
      RandFctWrapper(double (*rg)()) : m_rg(rg) {}
      virtual double generate() { return m_rg(); }
    protected:
      virtual ~RandFctWrapper();
      double (*m_rg)();
    };
    RandFctWrapper::~RandFctWrapper(){}

    void * & internal(void*o) {
      //object is here a pointer to a struct like ncrystal_xxx_t (which one is
      //not important since they all have the same layout):
      return ((ncrystal_info_t*)o)->internal;
    }
    Scatter * extract_scatter(ncrystal_scatter_t o) {
      nc_assert(o.internal);
      return reinterpret_cast<Scatter*>(o.internal);
    }
    Process * extract_process(ncrystal_process_t o) {
      nc_assert(o.internal);
      return reinterpret_cast<Process*>(o.internal);
    }
    Info * extract_info(ncrystal_info_t o) {
      nc_assert(o.internal);
      return reinterpret_cast<Info*>(o.internal);
    }
    RCBase * extract_rcbase(void* o) {
      nc_assert(internal(o));
      return reinterpret_cast<RCBase*>(internal(o));
    }
  }
}

namespace ncc = NCrystal::NCCInterface;
namespace NC = NCrystal;

void ncrystal_seterrhandler(void (*handler)(char*,char*))
{
  ncc::custom_error_handler = handler;
}

int ncrystal_error()
{
  return ncc::waserror;
}

const char * ncrystal_lasterror()
{
  return ncc::waserror ? ncc::errmsg : 0;
}

const char * ncrystal_lasterrortype()
{
  return ncc::waserror ? ncc::errtype : 0;
}

void ncrystal_clearerror()
{
  ncc::waserror = 0;
}

int ncrystal_setquietonerror(int q)
{
  int old = ncc::quietonerror;
  ncc::quietonerror = q;
  return old;
}

int ncrystal_sethaltonerror(int h)
{
  int old = ncc::haltonerror;
  ncc::haltonerror = h;
  return old;
}

int ncrystal_valid(void* object)
{
  if (!object)
    return 0;
  void *& i = ncc::internal(object);
  return i ? 1 : 0;
}

ncrystal_process_t ncrystal_cast_scat2proc(ncrystal_scatter_t s)
{
  //static cast, should always be possible:
  ncrystal_process_t p;
  p.internal = s.internal;
  return p;
}

ncrystal_process_t ncrystal_cast_abs2proc(ncrystal_absorption_t a)
{
  //static cast, should always be possible:
  ncrystal_process_t p;
  p.internal = a.internal;
  return p;
}

#define NCCATCH catch (std::exception& e) { ncc::handleError(e); }

int ncrystal_refcount(void* o)
{
  int rc(-999);
  if (!ncrystal_valid(o)) {
    ncc::setError("ncrystal_refcount called with invalid object");
    return rc;
  }
  try {
    int rc2 = ncc::extract_rcbase(o)->refCount();
    rc = rc2;
  } NCCATCH;
  return rc;
}

void ncrystal_ref(void* o)
{
  if (!ncrystal_valid(o)) {
    ncc::setError("ncrystal_ref called with invalid object");
    return;
  }
  try {
    ncc::extract_rcbase(o)->ref();
  } NCCATCH;
}

void ncrystal_unref(void* o)
{
  if (!ncrystal_valid(o)) {
    ncc::setError("ncrystal_unref called with invalid object");
    return;
  }
  try {
    NC::RCBase* rcb = ncc::extract_rcbase(o);
    unsigned rc = rcb->refCount();
    if (rc==1)
      ncc::internal(o) = 0;//make sure passed ncrystal_xxx_t is now invalid
    rcb->unref();
  } NCCATCH;
}

void ncrystal_unrefnodelete(void* o)
{
  if (!ncrystal_valid(o)) {
    ncc::setError("ncrystal_unrefnodelete called with invalid object");
    return;
  }
  try {
    NC::RCBase* rcb = ncc::extract_rcbase(o);
    rcb->unrefNoDelete();
  } NCCATCH;
}

void ncrystal_invalidate(void* o)
{
  if (!ncrystal_valid(o))
    return;
  ncc::internal(o) = 0;
}

void ncrystal_dump(ncrystal_info_t ci)
{
  if (!ncrystal_valid(&ci)) {
    ncc::setError("ncrystal_dump called with invalid info object");
    return;
  }
  try {
    NC::dump(ncc::extract_info(ci));
  } NCCATCH;
}


int ncrystal_info_getstructure( ncrystal_info_t ci_t,
                                unsigned* spacegroup,
                                double* lattice_a, double* lattice_b, double* lattice_c,
                                double* alpha, double* beta, double* gamma,
                                double* volume, unsigned* n_atoms )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getstructure called with invalid info object");
    return 0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    if (!ci->hasStructureInfo())
      return 0;
    const NC::StructureInfo& si = ci->getStructureInfo();
    *spacegroup = si.spacegroup;
    *lattice_a = si.lattice_a;
    *lattice_b = si.lattice_b;
    *lattice_c = si.lattice_c;
    *alpha = si.alpha;
    *beta = si.beta;
    *gamma = si.gamma;
    *volume = si.volume;
    *n_atoms = si.n_atoms;
    return 1;
  } NCCATCH;
  return 0;
}

double ncrystal_info_gettemperature( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_gettemperature called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasTemperature() ? ci->getTemperature() : -1;
  } NCCATCH;
  return -1;
}

double ncrystal_info_getxsectabsorption( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getxsectabsorption called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasXSectAbsorption() ? ci->getXSectAbsorption() : -1;
  } NCCATCH;
  return -1;
}

double ncrystal_info_getxsectfree( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getxsectfree called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasXSectFree() ? ci->getXSectFree() : -1;
  } NCCATCH;
  return -1;
}

double ncrystal_info_getglobaldebyetemp( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getglobaldebyetemp called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasGlobalDebyeTemperature() ? ci->getGlobalDebyeTemperature() : -1;
  } NCCATCH;
  return -1;
}

double ncrystal_info_getdensity( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getdensity called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasDensity() ? ci->getDensity() : -1;
  } NCCATCH;
  return -1;
}

double ncrystal_info_getnumberdensity( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getnumberdensity called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasNumberDensity() ? ci->getNumberDensity() : -1;
  } NCCATCH;
  return -1;
}

int ncrystal_info_nhkl( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_nhkl called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasHKLInfo() ? ci->nHKL() : -1;
  } NCCATCH;
  return -1;
}

double ncrystal_info_hkl_dlower( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_hkl_dlower called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hklDLower();
  } NCCATCH;
  return -1;
}

double ncrystal_info_hkl_dupper( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_hkl_dupper called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hklDUpper();
  } NCCATCH;
  return -1;
}

void ncrystal_info_gethkl( ncrystal_info_t ci_t, int idx,
                           int* h, int* k, int* l, int* multiplicity,
                           double * dspacing, double* fsquared )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_gethkl called with invalid info object");
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    NC::HKLList::const_iterator it = ci->hklBegin() + idx;
    nc_assert(it<ci->hklEnd());
    *h = it->h;
    *k = it->k;
    *l = it->l;
    *multiplicity = it->multiplicity;
    *dspacing = it->dspacing;
    *fsquared = it->fsquared;
  } NCCATCH;
}


unsigned ncrystal_info_ndyninfo( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_ndyninfo called with invalid info object");
    return 0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return static_cast<unsigned>(ci->getDynamicInfoList().size());
  } NCCATCH;
  return 0;
}


void ncrystal_dyninfo_base( ncrystal_info_t ci_t,
                            unsigned idyninfo,
                            double* fraction,
                            const char** elementname,
                            double* temperature,
                            unsigned* ditypeid )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_dyninfo_base called with invalid info object");
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    auto& di = ci->getDynamicInfoList().at(idyninfo);
    *fraction = di->fraction();
    *temperature = di->temperature();
    *elementname = di->elementName().c_str();
    if (dynamic_cast<const NC::DI_Sterile*>(di.get()))
      *ditypeid = 0;
    else if (dynamic_cast<const NC::DI_FreeGas*>(di.get()))
      *ditypeid = 1;
    else if (dynamic_cast<const NC::DI_ScatKnlDirect*>(di.get()))
      *ditypeid = 2;
    else if (dynamic_cast<const NC::DI_VDOS*>(di.get()))
      *ditypeid = 3;
    else if (dynamic_cast<const NC::DI_VDOSDebye*>(di.get()))
      *ditypeid = 4;
    else
      *ditypeid = 99;
  } NCCATCH;
}

void ncrystal_dyninfo_extract_scatknl( ncrystal_info_t ci_t,
                                       unsigned idyninfo,
                                       unsigned vdoslux,
                                       double * suggestedEmax,
                                       unsigned* negrid,
                                       unsigned* nalpha,
                                       unsigned* nbeta,
                                       const double** egrid,
                                       const double** alphagrid,
                                       const double** betagrid,
                                       const double** sab )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_dyninfo_extract_scatknl called with invalid info object");
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    auto& di = ci->getDynamicInfoList().at(idyninfo);
    nc_assert_always(!!di);
    std::shared_ptr<const NC::SABData> shptr_sabdata;
    NC::DI_ScatKnl::EGridShPtr shptr_egrid;
    auto di_sk = dynamic_cast<const NC::DI_ScatKnl*>(di.get());
    if (di_sk) {
      shptr_sabdata = NC::extractSABDataFromDynInfo(di_sk,vdoslux);
      shptr_egrid = di_sk->energyGrid();
      //In case the sabdata factory does not keep strong references, we must
      //keep the newly created object in shptr_sabdata alive when returning to
      //C/Python code without shared pointers. For now we do this by adding to a
      //global static array here:
      static std::vector<std::shared_ptr<const NC::SABData>> s_keepAlive;
      static std::mutex s_keepAlive_mutex;
      std::lock_guard<std::mutex> guard(s_keepAlive_mutex);
      s_keepAlive.push_back(shptr_sabdata);
      static bool first = true;
      if (first) {
        //Register for clearance by global clearCaches function:
        first = false;
        NC::registerCacheCleanupFunction([](){ s_keepAlive.clear(); });
      }

    }
    if (!!shptr_sabdata) {
      const auto& sabdata = *shptr_sabdata;
      unsigned na = static_cast<unsigned>(sabdata.alphaGrid().size());
      unsigned nb = static_cast<unsigned>(sabdata.betaGrid().size());
      unsigned nsab = static_cast<unsigned>(sabdata.sab().size());
      nc_assert_always(na>1&&nb>1&&na*nb==nsab);
      *nalpha = na;
      *nbeta = nb;
      *alphagrid = &sabdata.alphaGrid()[0];
      *betagrid = &sabdata.betaGrid()[0];
      *sab = &sabdata.sab()[0];
      *suggestedEmax = sabdata.suggestedEmax();
    } else {
      *nalpha = 0;
      *nbeta = 0;
      *alphagrid = nullptr;
      *betagrid = nullptr;
      *sab = nullptr;
      *suggestedEmax = 0.0;
    }
    if ( !shptr_egrid || shptr_egrid->empty() ) {
      *negrid = 0;
      //make sure egrid does at least point to valid memory (and avoid actually
      //using 0-length arrays);
      static const double dummy[1] = { 0.0 };
      *egrid = &dummy[0];
    } else {
      *negrid = shptr_egrid->size();
      *egrid = &(*shptr_egrid)[0];
    }
  } NCCATCH;
}

void ncrystal_dyninfo_extract_vdos( ncrystal_info_t ci_t,
                                    unsigned idyninfo,
                                    double * egridMin,
                                    double * egridMax,
                                    unsigned * vdos_ndensity,
                                    const double ** vdos_density )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_dyninfo_extract_vdos called with invalid info object");
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    auto& di = ci->getDynamicInfoList().at(idyninfo);
    auto di_vdos = dynamic_cast<const NC::DI_VDOS*>(di.get());
    if (di_vdos) {
      const auto& vdosData = di_vdos->vdosData();
      const auto& v_egrid = vdosData.vdos_egrid();
      const auto& v_density = vdosData.vdos_density();
      nc_assert_always(v_density.size()<=std::numeric_limits<unsigned>::max());
      *egridMin = v_egrid.first;
      *egridMax = v_egrid.second;
      *vdos_ndensity = static_cast<unsigned>(v_density.size());
      *vdos_density = &v_density[0];
    } else {
      static const double dummy[1] = { 0.0 };//avoid pointing to invalid mem
      *egridMin = 0;
      *egridMax = 0;
      *vdos_ndensity = 0;
      *vdos_density = &dummy[0];
    }
  } NCCATCH;
}

NCRYSTAL_API void ncrystal_dyninfo_extract_vdos_input( ncrystal_info_t ci_t,
                                                       unsigned idyninfo,
                                                       unsigned* vdos_negrid,
                                                       const double ** vdos_egrid,
                                                       unsigned* vdos_ndensity,
                                                       const double ** vdos_density )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_dyninfo_extract_vdos_input called with invalid info object");
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    auto& di = ci->getDynamicInfoList().at(idyninfo);
    auto di_vdos = dynamic_cast<const NC::DI_VDOS*>(di.get());
    static const double dummy[1] = { 0.0 };//avoid pointing to invalid mem
    *vdos_negrid = 0;
    *vdos_ndensity = 0;
    *vdos_egrid = &dummy[0];
    *vdos_density = &dummy[0];
    if (di_vdos) {
      const NC::VectD& egrid = di_vdos->vdosOrigEgrid();
      const NC::VectD& density = di_vdos->vdosOrigDensity();
      nc_assert_always(density.size()<=std::numeric_limits<unsigned>::max());
      if ( !egrid.empty() && !density.empty() ) {
        *vdos_egrid = &egrid[0];
        *vdos_density = &density[0];
        *vdos_negrid = static_cast<unsigned>(egrid.size());
        *vdos_ndensity = static_cast<unsigned>(density.size());
      }
    }
  } NCCATCH;
}

void ncrystal_dyninfo_extract_vdosdebye( ncrystal_info_t ci_t,
                                         unsigned idyninfo,
                                         double * debye_temp )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_dyninfo_extract_vdosdebye called with invalid info object");
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    auto& di = ci->getDynamicInfoList().at(idyninfo);
    auto di_vdosdebye = dynamic_cast<const NC::DI_VDOSDebye*>(di.get());
    *debye_temp = di_vdosdebye ? di_vdosdebye->debyeTemperature() : 0.0;
  } NCCATCH;
}


double ncrystal_info_dspacing_from_hkl( ncrystal_info_t ci_t, int h, int k, int l )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_dspacing_from_hkl called with invalid info object");
    return 0.0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->dspacingFromHKL(h,k,l);
  } NCCATCH;
  return 0.0;
}

void ncrystal_setrandgen( double (*rg)() )
{
  try {
    NC::setDefaultRandomGenerator( rg ? new ncc::RandFctWrapper(rg) : 0);
  } NCCATCH;
}

void ncrystal_save_randgen()
{
  try {
    if (ncc::saved_rng.obj()) {
      ncc::setError("ncrystal_save_randgen called when a state is already saved");
      return;
    }
    ncc::saved_rng = NC::defaultRandomGenerator(false);
  } NCCATCH;
}

void ncrystal_restore_randgen()
{
  try {
    NC::RandomBase * rng = ncc::saved_rng.obj();
    NC::RCGuard guard(rng);
    ncc::saved_rng = 0;
    NC::setDefaultRandomGenerator( rng );
  } NCCATCH;
}

void ncrystal_setbuiltinrandgen()
{
  try {
    NC::RandomBase * rng = new NC::RandXRSR();
    NC::RCGuard guard(rng);
    NC::setDefaultRandomGenerator( rng );
  } NCCATCH;
}


double ncrystal_wl2ekin( double wl )
{
  return NC::wl2ekin(wl);//doesnt throw
}

double ncrystal_ekin2wl( double ekin )
{
  return NC::ekin2wl(ekin);//doesnt throw
}

void ncrystal_natelemdata( unsigned z, const char ** name,
                           double* mass_amu, double* sigma_inc,
                           double* sigma_coh, double* sigma_abs )
{
  try {
    auto db = NC::NeutronSCL::instance();
    const std::string& name_str = db->getAtomName(z);
    *name = name_str.c_str();
    if ( name_str.empty() ) {
      *mass_amu = 0.0;
      *sigma_inc = 0.0;
      *sigma_coh = 0.0;
      *sigma_abs = 0.0;
    } else {
      *mass_amu  = db->getAtomicMass(name_str);
      *sigma_inc = db->getIncoherentXS(name_str);
      *sigma_coh = db->getCoherentXS(name_str);
      *sigma_abs = db->getCaptureXS(name_str);
    }
  } NCCATCH;
}



int ncrystal_isnonoriented(ncrystal_process_t o)
{
  NC::Process * process = ncc::extract_process(o);
  if (!process) {
    ncc::setError("ncrystal_isnonoriented called with invalid object");
    return 0;
  }
  try {
    return process->isOriented() ? 0 : 1;
  } NCCATCH;
  return 0;
}

const char * ncrystal_name(ncrystal_process_t o)
{
  NC::Process * process = ncc::extract_process(o);
  if (!process) {
    ncc::setError("ncrystal_name called with invalid object");
    return 0;
  }
  try {
    return process->getCalcName();
  } NCCATCH;
  return 0;
}

void ncrystal_domain( ncrystal_process_t o,
                      double* ekin_low, double* ekin_high)
{
  NC::Process * process = ncc::extract_process(o);
  if (!process) {
    ncc::setError("ncrystal_domain called with invalid object");
    return;
  }
  try {
    process->domain(*ekin_low,*ekin_high);
  } NCCATCH;
}


void ncrystal_crosssection_nonoriented( ncrystal_process_t o, double ekin, double* result)
{
  NC::Process * process = ncc::extract_process(o);
  if (!process) {
    ncc::setError("ncrystal_crosssection_nonoriented called with invalid object");
    return;
  }
  try {
    *result = process->crossSectionNonOriented(ekin);
  } NCCATCH;
}

void ncrystal_genscatter_nonoriented( ncrystal_scatter_t o, double ekin, double* result_angle, double* result_dekin )
{
  NC::Scatter * scatter = ncc::extract_scatter(o);
  if (!scatter) {
    ncc::setError("ncrystal_genscatter_nonoriented called with invalid object");
    return;
  }
  try {
    scatter->generateScatteringNonOriented( ekin, *result_angle, *result_dekin );
  } NCCATCH;
}

void ncrystal_genscatter_nonoriented_many( ncrystal_scatter_t o,
                                           const double * ekin,
                                           unsigned long n_ekin,
                                           unsigned long repeat,
                                           double* results_angle,
                                           double* results_dekin )
{
  NC::Scatter * scatter = ncc::extract_scatter(o);
  if (!scatter) {
    ncc::setError("ncrystal_genscatter_nonoriented_many called with invalid object");
    return;
  }
  try {
    while (repeat--) {
      for (unsigned long i = 0; i < n_ekin; ++i)
        scatter->generateScatteringNonOriented(ekin[i],*results_angle++,*results_dekin++);
    }
  } NCCATCH;
}

void ncrystal_crosssection_nonoriented_many( ncrystal_process_t o,
                                             const double * ekin,
                                             unsigned long n_ekin,
                                             unsigned long repeat,
                                             double* results )
{
  NC::Process * process = ncc::extract_process(o);
  if (!process) {
    ncc::setError("ncrystal_crosssection_nonoriented_many called with invalid object");
    return;
  }
  try {
    while (repeat--) {
      for (unsigned long i = 0; i < n_ekin; ++i)
        *results++ = process->crossSectionNonOriented(ekin[i]);
    }
  } NCCATCH;
}

void ncrystal_crosssection( ncrystal_process_t o, double ekin, const double (*direction)[3], double* result)
{
  *result = -1.0;
  NC::Process * process = ncc::extract_process(o);
  if (!process) {
    ncc::setError("ncrystal_crosssection called with invalid object");
    return;
  }
  try {
    double r = process->crossSection( ekin, *direction );
    *result = r;
  } NCCATCH;
}

void ncrystal_genscatter( ncrystal_scatter_t o, double ekin, const double (*direction)[3],
                          double (*result_direction)[3], double* result_deltaekin )
{
  NC::Scatter * scatter = ncc::extract_scatter(o);
  if (!scatter) {
    ncc::setError("ncrystal_genscatter called with invalid object");
    *result_direction[0] = *result_direction[1] = *result_direction[2] = 0.0;
    *result_deltaekin = 0.0;
    return;
  }
  try {
    scatter->generateScattering( ekin, *direction, *result_direction, *result_deltaekin );
  } catch (std::exception& e) {
    *result_direction[0] = *result_direction[1] = *result_direction[2] = 0.0;
    *result_deltaekin = 0.0;
    ncc::handleError(e);
  }
}

ncrystal_info_t ncrystal_create_info( const char * cfgstr )
{
  ncrystal_info_t o;
  o.internal = 0;
  try {
    const NC::Info * info = NC::createInfo(cfgstr);
    nc_assert(info);
    info->ref();
    o.internal = (void*)info;
  } NCCATCH;
  return o;
}

double ncrystal_decodecfg_packfact( const char * cfgstr )
{
  try {
    NC::MatCfg cfg(cfgstr);
    return cfg.get_packfact();
  } NCCATCH;
  return -1.0;
}

unsigned ncrystal_decodecfg_vdoslux( const char * cfgstr )
{
  try {
    NC::MatCfg cfg(cfgstr);
    return cfg.get_vdoslux();
  } NCCATCH;
  return 999;
}

ncrystal_scatter_t ncrystal_create_scatter( const char * cfgstr )
{
  ncrystal_scatter_t o;
  o.internal = 0;
  try {
    const NC::Scatter * scatter = NC::createScatter(cfgstr);
    nc_assert(scatter);
    scatter->ref();
    o.internal = (void*)scatter;
  } NCCATCH;
  return o;
}

ncrystal_absorption_t ncrystal_create_absorption( const char * cfgstr )
{
  ncrystal_absorption_t o;
  o.internal = 0;
  try {
    const NC::Absorption * absorption = NC::createAbsorption(cfgstr);
    nc_assert(absorption);
    absorption->ref();
    o.internal = (void*)absorption;
  } NCCATCH;
  return o;
}

void ncrystal_clear_info_caches()
{
  try {
    NC::clearInfoCaches();
  } NCCATCH;
}

void ncrystal_clear_caches()
{
  try {
    NC::clearCaches();
  } NCCATCH;
}

void ncrystal_disable_caching()
{
  try {
    NC::disableCaching();
  } NCCATCH;
}

void ncrystal_enable_caching()
{
  try {
    NC::enableCaching();
  } NCCATCH;
}

void ncrystal_clear_factory_registry()
{
  try {
    NC::clearFactoryRegistry();
  } NCCATCH;
}

int ncrystal_has_factory( const char* name )
{
  int res = 0;
  try {
    res = NC::hasFactory(name) ? 1 : 0;
  } NCCATCH;
  return res;
}

int ncrystal_version()
{
  return NCRYSTAL_VERSION;
}

const char * ncrystal_version_str()
{
  return NCRYSTAL_VERSION_STR;
}
