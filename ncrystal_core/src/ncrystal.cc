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
#include "NCrystal/internal/NCDynInfoUtils.hh"
#include "NCrystal/NCDump.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCAtomUtils.hh"
#include "NCrystal/internal/NCAtomDB.hh"
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

    struct AtomWrapper : public RCBase {

      //Ref-counted instance of AtomDataSP. This allows C/Python side to
      //par-take in life-time management. The wrapper is associated with a
      //particular Info instance, and therefore also knows the displayLabel
      //(unless a sub-component).
      AtomDataSP atomDataSP;

      //Cache display labels and descriptions (without values) here (note that
      //the displayLabels are only valid in connection with a particular Info
      //object, and only "top-level" atoms, the ones with an AtomIndex on the
      //Info object, have a displayLabel.

      const std::string& displayLabel() const
      {
        static std::string s_empty;
        return !displayLabel_ptr ? s_empty : *displayLabel_ptr;
      }
      const std::string& description() const
      {
        nc_assert( description_ptr != nullptr || displayLabel_ptr != nullptr );
        return description_ptr!=nullptr ? *description_ptr : *displayLabel_ptr;
      }

      std::unique_ptr<std::string> displayLabel_ptr;
      std::unique_ptr<std::string> description_ptr;

      //NB: Consider MT safety carefully, createInfo(..) is protected by mutex,
      //but the C-interface code is not (but could be...)  We could have
      //InfoWrapper with a mutex and caches?
    };

    class RandFctWrapper final : public RandomBase {
    public:
      RandFctWrapper(double (*rg)()) : m_rg(rg) {}
      virtual double generate() { return m_rg(); }
    protected:
      virtual ~RandFctWrapper() = default;
      double (*m_rg)();
    };

    void * & internal(void*o) {
      //object is here a pointer to a struct like ncrystal_xxx_t (which one is
      //not important since they all have the same layout):
      return ((ncrystal_info_t*)o)->internal;
    }
    template <class T,class THandle>
    T* doExtract(THandle o)
    {
      //Extract any of the handles (relying on the fact that they all wrap
      //something derived from RCBase). In debug builds we can detect errors in
      //user-code where the wrong kind of handle is provided to the function.
#ifndef NDEBUG
      if ( !o.internal || !dynamic_cast<T*>(reinterpret_cast<RCBase*>(o.internal)) )
        return nullptr;
#endif
      return reinterpret_cast<T*>(o.internal);
    }
    Scatter * extract_scatter(ncrystal_scatter_t o) {
      return doExtract<Scatter,ncrystal_scatter_t>(o);
    }
    Process * extract_process(ncrystal_process_t o) {
      return doExtract<Process,ncrystal_process_t>(o);
    }
    Info * extract_info(ncrystal_info_t o) {
      return doExtract<Info,ncrystal_info_t>(o);
    }
    RCBase * extract_rcbase(void* o) {
      nc_assert(o!=nullptr);
      return reinterpret_cast<RCBase*>(internal(o));
    }
    AtomWrapper * extract_atomwrapper(ncrystal_atomdata_t o) {
      return doExtract<AtomWrapper,ncrystal_atomdata_t>(o);
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
                            unsigned* atomdataindex,
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
    *atomdataindex = di->atom().index.value;
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
  return NC::wl2ekin(wl);//doesn't throw
}

double ncrystal_ekin2wl( double ekin )
{
  return NC::ekin2wl(ekin);//doesn't throw
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

void ncrystal_genscatter_many( ncrystal_scatter_t o,
                               double ekin,
                               const double (*indir)[3],
                               unsigned long repeat,
                               double * results_dirx,
                               double * results_diry,
                               double * results_dirz,
                               double * results_dekin )
{
  NC::Scatter * scatter = ncc::extract_scatter(o);
  if (!scatter) {
    ncc::setError("ncrystal_genscatter_many called with invalid object");
    return;
  }
  try {
    double outdir[3];
    while (repeat--) {
      scatter->generateScattering(ekin,*indir,outdir,*results_dekin++);
      *results_dirx++ = outdir[0];
      *results_diry++ = outdir[1];
      *results_dirz++ = outdir[2];
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

NCRYSTAL_API unsigned ncrystal_info_natominfo( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_natominfo called with invalid info object");
    return 0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    if (!ci->hasAtomInfo())
      return 0;
    return static_cast<unsigned>(std::distance(ci->atomInfoBegin(),ci->atomInfoEnd()));
  } NCCATCH;
  return 0;
}

NCRYSTAL_API int ncrystal_info_hasatompos( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_hasatompos called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasAtomPositions() ? 1 : 0;
  } NCCATCH;
  return 0;
}

NCRYSTAL_API int ncrystal_info_hasatommsd( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_hasatommsd called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasAtomMSD() ? 1 : 0;
  } NCCATCH;
  return 0;
}

NCRYSTAL_API void ncrystal_info_getatominfo( ncrystal_info_t ci_t, unsigned iatom,
                                             unsigned* atomdataindex,
                                             unsigned* number_per_unit_cell,
                                             double* debye_temp, double* msd )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getatominfo called with invalid info object");
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    if (iatom >= std::distance(ci->atomInfoBegin(),ci->atomInfoEnd()))
      NCRYSTAL_THROW(BadInput,"ncrystal_info_getatominfo iatom is out of bounds");
    const NC::AtomInfo& ai = *std::next(ci->atomInfoBegin(),iatom);
    *atomdataindex = ai.atom.index.value;
    *number_per_unit_cell = ai.number_per_unit_cell;
    *debye_temp = ai.debye_temp;
    if ( *debye_temp == 0.0 && ci->hasGlobalDebyeTemperature() )
      *debye_temp = ci->getGlobalDebyeTemperature();//convenience - not yet done in C++ [todo]
    *msd = ai.mean_square_displacement;
  } NCCATCH;
}

NCRYSTAL_API void ncrystal_info_getatompos( ncrystal_info_t ci_t,
                                            unsigned iatom, unsigned ipos,
                                            double* x, double* y, double* z )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getatompos called with invalid info object");
    *x = *y = *z = -999.0;
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    if (iatom >= std::distance(ci->atomInfoBegin(),ci->atomInfoEnd()))
      NCRYSTAL_THROW(BadInput,"ncrystal_info_getatominfo iatom is out of bounds");
    const NC::AtomInfo& ai = *std::next(ci->atomInfoBegin(),iatom);
    if (ai.positions.empty())
      NCRYSTAL_THROW(BadInput,"ncrystal_info_getatompos called but positions not available");
    if ( ! (ipos<ai.positions.size()) )
      NCRYSTAL_THROW(BadInput,"ncrystal_info_getatominfo ipos is out of bounds");
    const auto& pos = ai.positions[ipos];
    *x = pos.x;
    *y = pos.y;
    *z = pos.z;
    return;
  } NCCATCH;
  *x = *y = *z = -999.0;
}

NCRYSTAL_API int ncrystal_info_hasanydebyetemp( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_hasanydebyetemp called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->hasAnyDebyeTemperature() ? 1 : 0;
  } NCCATCH;
  return -1;
}

NCRYSTAL_API double ncrystal_info_getdebyetempbyelement( ncrystal_info_t ci_t,
                                                         unsigned atomdataindex )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getdebyetempbyelemname called with invalid info object");
    return -1;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return ci->getDebyeTemperatureByElement(NC::AtomIndex{atomdataindex});
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

unsigned ncrystal_info_ncustomsections( ncrystal_info_t ci_t )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_ncustomsections called with invalid info object");
    return 0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    return static_cast<unsigned>( ci->getAllCustomSections().size() );
  } NCCATCH;
  return 0;
}

const char* ncrystal_info_customsec_name( ncrystal_info_t ci_t, unsigned isection )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_customsec_name called with invalid info object");
    return "";
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    const auto& cd = ci->getAllCustomSections();
    return cd.at(isection).first.c_str();
  } NCCATCH;
  return "";
}

unsigned ncrystal_info_customsec_nlines( ncrystal_info_t ci_t, unsigned isection )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_customsec_nlines called with invalid info object");
    return 0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    const auto& cd = ci->getAllCustomSections();
    return static_cast<unsigned>( cd.at(isection).second.size() );
  } NCCATCH;
  return 0;
}

unsigned ncrystal_info_customline_nparts( ncrystal_info_t ci_t, unsigned isection, unsigned iline )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_customline_nparts called with invalid info object");
    return 0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    const auto& cd = ci->getAllCustomSections();
    return static_cast<unsigned>( cd.at(isection).second.at(iline).size() );
  } NCCATCH;
  return 0;
}

const char* ncrystal_info_customline_getpart( ncrystal_info_t ci_t, unsigned isection, unsigned iline, unsigned ipart )
{
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_customline_getpart called with invalid info object");
    return "";
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    const auto& cd = ci->getAllCustomSections();
    return cd.at(isection).second.at(iline).at(ipart).c_str();
  } NCCATCH;
  return "";
}

void ncrystal_register_in_mem_file_data(const char* virtual_filename, const char* cdata)
{
  try {
    std::string name(virtual_filename);
    NC::registerInMemoryFileData(name,std::string(cdata));
  } NCCATCH;
}


ncrystal_atomdata_t ncrystal_create_atomdata( ncrystal_info_t ci_t,
                                              unsigned atomdataindex )
{
  ncrystal_atomdata_t o;
  o.internal = nullptr;
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_create_atomdata called with invalid info object");
    return o;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    NC::RCHolder<ncc::AtomWrapper> wrapper_holder(new ncc::AtomWrapper);
    auto& wrapper = *wrapper_holder.obj();
    wrapper.atomDataSP = ci->atomDataSP(NC::AtomIndex{atomdataindex});
    nc_assert(!!wrapper.atomDataSP);
    auto dlbl = ci->displayLabel(NC::AtomIndex{atomdataindex});
    auto descr = wrapper.atomDataSP->description(false);
    wrapper.displayLabel_ptr = std::make_unique<std::string>(dlbl);
    if ( dlbl!=descr )
      wrapper.description_ptr = std::make_unique<std::string>(descr);
    nc_assert(wrapper.displayLabel()==dlbl);
    nc_assert(wrapper.description()==descr);
    //return with ref count of 1:
    wrapper.ref();
    o.internal = (void*)wrapper_holder.obj();
  } NCCATCH;
  return o;
}

void ncrystal_atomdata_getfields( ncrystal_atomdata_t o,
                                  const char** displaylabel,
                                  const char** description,
                                  double* mass, double *incxs,
                                  double* cohsl_fm, double* absxs,
                                  unsigned* ncomponents,
                                  unsigned* zval, unsigned* aval )
{
  ncc::AtomWrapper * atomwrapper = ncc::extract_atomwrapper(o);
  if (!atomwrapper) {
    ncc::setError("ncrystal_atomdata_getfields called with invalid object");
    *displaylabel = *description = nullptr;
    *mass = *incxs = *cohsl_fm = *absxs = 0.0;
    *ncomponents = *zval = *aval = 0;
    return;
  }
  try {
    *displaylabel = atomwrapper->displayLabel().c_str();
    *description = atomwrapper->description().c_str();
    nc_assert(!!atomwrapper->atomDataSP);
    const NC::AtomData& data = *atomwrapper->atomDataSP;
    *mass = data.averageMassAMU();
    *cohsl_fm = data.coherentScatLenFM();
    *incxs = data.incoherentXS().val;
    *absxs = data.captureXS();
    *zval = data.isElement() ? data.Z() : 0;
    *aval = data.isSingleIsotope() ? data.A() : 0;
    *ncomponents = data.nComponents();
  } NCCATCH;
}

ncrystal_atomdata_t ncrystal_create_atomdata_subcomp( ncrystal_atomdata_t ad_t,
                                                      unsigned icomponent,
                                                      double* fraction )
{
  ncrystal_atomdata_t o;
  o.internal = nullptr;
  *fraction = -1.0;
  ncc::AtomWrapper * atomwrapper_parent = ncc::extract_atomwrapper(ad_t);
  if (!atomwrapper_parent) {
    ncc::setError("ncrystal_create_atomdata_subcomp called with invalid object");
    return o;
  }
  try {
    nc_assert(!!atomwrapper_parent->atomDataSP);
    const auto& comp = atomwrapper_parent->atomDataSP->getComponent(icomponent);
    NC::RCHolder<ncc::AtomWrapper> wrapper_holder(new ncc::AtomWrapper);
    auto& wrapper = *wrapper_holder.obj();
    wrapper.atomDataSP = comp.data;
    nc_assert(!!wrapper.atomDataSP);
    auto descr = wrapper.atomDataSP->description(false);
    wrapper.description_ptr = std::make_unique<std::string>(descr);
    nc_assert(wrapper.displayLabel().empty());
    nc_assert(wrapper.description()==descr);
    //return with ref count of 1:
    wrapper.ref();
    o.internal = (void*)wrapper_holder.obj();
    *fraction = comp.fraction;
  } NCCATCH;
  return o;
}

unsigned ncrystal_info_ncomponents( ncrystal_info_t ci_t )
{
  //0 means !hasComposition()
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_ncomponents called with invalid info object");
    return 0;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    if (!ci->hasComposition())
      return 0;
    auto n = ci->getComposition().size();
    nc_assert( n>0 && n < std::numeric_limits<unsigned>::max() );
    return static_cast<unsigned>(n);
  } NCCATCH;
  return 0;
}

void ncrystal_info_getcomponent( ncrystal_info_t ci_t, unsigned icomponent,
                                 unsigned* atomdataindex, double* fraction )
{
  *atomdataindex = 999999;
  *fraction = -1.0;
  if (!ncrystal_valid(&ci_t)) {
    ncc::setError("ncrystal_info_getcomponent called with invalid info object");
    return;
  }
  try {
    NC::Info * ci = ncc::extract_info(ci_t);
    auto n = ci->hasComposition() ? ci->getComposition().size() : 0;
    if ( ! ( icomponent<n) )
      NCRYSTAL_THROW(BadInput,"Requested component index is out of bounds");
    const auto& comp = ci->getComposition().at(icomponent);
    *atomdataindex = comp.atom.index.value;
    *fraction = comp.fraction;
  } NCCATCH;
}

ncrystal_atomdata_t ncrystal_create_atomdata_fromdb( unsigned z, unsigned a )
{
  ncrystal_atomdata_t o;
  o.internal = nullptr;
  try {
    NC::RCHolder<ncc::AtomWrapper> wrapper_holder(new ncc::AtomWrapper);
    auto& wrapper = *wrapper_holder.obj();
    wrapper.atomDataSP = NC::AtomDB::getIsotopeOrNatElem(z,a);
    if (!wrapper.atomDataSP)
      return o;//return invalid.
    auto descr = wrapper.atomDataSP->description(false);
    wrapper.description_ptr = std::make_unique<std::string>(descr);
    nc_assert(wrapper.displayLabel().empty());
    nc_assert(wrapper.description()==descr);
    //return with ref count of 1:
    wrapper.ref();
    o.internal = (void*)wrapper_holder.obj();
  } NCCATCH;
  return o;
}

ncrystal_atomdata_t ncrystal_create_atomdata_fromdbstr( const char* name )
{
  unsigned z(0),a(0);
  try {
    nc_assert(name);
    NC::AtomSymbol symb(name);
    if (symb.isElement()||symb.isIsotope()) {
      z = symb.Z();
      a = symb.A();
    }
  } NCCATCH;
  if (z)
    return ncrystal_create_atomdata_fromdb( z,a );
  return ncrystal_atomdata_t{nullptr};
}

unsigned ncrystal_atomdatadb_getnentries()
{
  return NC::AtomDB::getAllEntriesCount();
}

void ncrystal_atomdatadb_getallentries( unsigned* zvals,
                                        unsigned* avals )
{
  std::vector<std::pair<unsigned,unsigned>> all = NC::AtomDB::getAllEntries();
  nc_assert( static_cast<unsigned>(all.size()) == ncrystal_atomdatadb_getnentries() );
  for ( auto& e : all ) {
    *zvals++ = e.first;
    *avals++ = e.second;
  }
}
