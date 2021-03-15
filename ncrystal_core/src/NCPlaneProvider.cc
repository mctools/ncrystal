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

#include "NCrystal/internal/NCPlaneProvider.hh"
#include "NCrystal/internal/NCOrientUtils.hh"
#include "NCrystal/NCMatInfo.hh"
#include "NCrystal/internal/NCRotMatrix.hh"
#include "NCrystal/internal/NCEqRefl.hh"

namespace NCrystal {

  PlaneProvider::PlaneProvider() = default;
  PlaneProvider::~PlaneProvider() = default;

  class PlaneProviderStd final : public PlaneProvider {
  public:

    PlaneProviderStd(shared_obj<const MatInfo>);
    PlaneProviderStd( const MatInfo * );
    virtual ~PlaneProviderStd() = default;

    bool canProvide() const final;
    void prepareLoop() final;
    bool getNextPlane(double& dspacing, double& fsq, Vector& demi_normal) final;

  private:
    optional_shared_obj<const MatInfo> m_info_strongref;
    const MatInfo* m_info;
    enum{ STRAT_MISSING, STRAT_DEMINORMAL, STRAT_EXPHKL, STRAT_SPACEGROUP } m_strategy;
    //outer loop:
    HKLList::const_iterator m_it_hklE;
    HKLList::const_iterator m_it_hkl;
    //inner loop counter (common for STRAT_DEMINORMAL + STRAT_EXPHKL)
    size_t m_ii;
    RotMatrix m_reci_lattice;
    //needed just for STRAT_SPACEGROUP:
    struct StrSG;
    std::unique_ptr<StrSG> m_sg;
    bool gnp_de(double& dspacing, double& fsq, Vector& normal);
    bool gnp_eh(double& dspacing, double& fsq, Vector& normal);
    bool gnp_sg(double& dspacing, double& fsq, Vector& normal);
  };

  struct PlaneProviderStd::StrSG {
    StrSG(int spacegroup) : m_eqreflcalc(spacegroup) {}
    void prepareLoop(int h, int k, int l, unsigned expected_multiplicity ) {
      const std::set<EqRefl::HKL>& el = m_eqreflcalc.getEquivalentReflections(h,k,l);
      if ( el.size() * 2 != expected_multiplicity ) {
        NCRYSTAL_THROW2(MissingInfo,"Incomplete information for selected modeling: Neither"
                        " HKL normals nor expanded HKL info available, and the HKL grouping in the"
                        " input does not appear to have the multiplicities expected of symmetry"
                        " equivalent families ( h,k,l="<<h<<","<<k<<","<<l
                        <<" had multiplicity of "<<expected_multiplicity<<" where "
                        <<el.size() * 2<<" was expected).");
      }
      it = el.begin();
      itE = el.end();
    }
    std::set<EqRefl::HKL>::const_iterator it,itE;
  private:
    EqRefl m_eqreflcalc;
  };

  PlaneProviderStd::PlaneProviderStd(shared_obj<const MatInfo> cinfo)
    : PlaneProviderStd(cinfo.get())
  {
    m_info_strongref = std::move(cinfo);
  }

  PlaneProviderStd::PlaneProviderStd(const MatInfo* cinfo)
    : PlaneProvider(),
      m_info(cinfo),
      m_strategy(STRAT_MISSING),
      m_ii(0)
  {
    nc_assert(m_info);
    if (m_info->hasHKLInfo()) {
      m_it_hkl  = m_info->hklBegin();
      m_it_hklE = m_info->hklEnd();
      if ( m_info->hasHKLDemiNormals() ) {
        m_strategy = STRAT_DEMINORMAL;
      } else if ( m_info->hasExpandedHKLInfo() ) {
        m_strategy = STRAT_EXPHKL;
      } else if ( m_info->hasStructureInfo() && m_info->getStructureInfo().spacegroup ) {
        m_strategy = STRAT_SPACEGROUP;
        if (m_it_hkl!=m_it_hklE)
          m_sg = std::make_unique<StrSG>(m_info->getStructureInfo().spacegroup);
      }
    }
    if ( m_strategy == STRAT_EXPHKL || m_strategy == STRAT_SPACEGROUP )
      m_reci_lattice = getReciprocalLatticeRot( *m_info );
    if (canProvide())
      prepareLoop();
  }

  bool PlaneProviderStd::canProvide() const
  {
    return m_strategy!=STRAT_MISSING;
  }

  void PlaneProviderStd::prepareLoop()
  {
    if (!canProvide())
      NCRYSTAL_THROW(MissingInfo,"Insufficient information for plane normals: Neither"
                     " HKL normals, expanded HKL info, or spacegroup number is available.");
    nc_assert(m_info);
    m_ii = 0;
    m_it_hkl  = m_info->hklBegin();
    m_it_hklE = m_info->hklEnd();
    if ( m_sg ) {
      nc_assert(m_strategy == STRAT_SPACEGROUP);
      nc_assert(m_it_hkl!=m_it_hklE);
      m_sg->prepareLoop(m_it_hkl->h, m_it_hkl->k, m_it_hkl->l, m_it_hkl->multiplicity);
    }
  }

  bool PlaneProviderStd::getNextPlane(double& dspacing, double& fsq, Vector& demi_normal)
  {
    switch(m_strategy) {
    case STRAT_DEMINORMAL: return gnp_de(dspacing,fsq,demi_normal);
    case STRAT_EXPHKL: return gnp_eh(dspacing,fsq,demi_normal);
    case STRAT_SPACEGROUP: return gnp_sg(dspacing,fsq,demi_normal);
    case STRAT_MISSING:
      NCRYSTAL_THROW(MissingInfo,"Insufficient information for plane normals: Neither"
                     " HKL normals, expanded HKL info, or spacegroup number is available.");
    };
    return false;
  }

  bool PlaneProviderStd::gnp_de(double& dspacing, double& fsq, Vector& demi_normal)
  {
    if (m_it_hkl == m_it_hklE)
      return false;
    if (m_ii == m_it_hkl->demi_normals.size()) {
      ++m_it_hkl;
      m_ii = 0;
      return gnp_de(dspacing,fsq,demi_normal);
    }
    const HKLInfo::Normal & nn = m_it_hkl->demi_normals.at(m_ii++);
    dspacing = m_it_hkl->dspacing;
    fsq = m_it_hkl->fsquared;
    demi_normal = nn.as<Vector>();
    return true;
  }

  bool PlaneProviderStd::gnp_eh(double& dspacing, double& fsq, Vector& demi_normal)
  {
    if (m_it_hkl == m_it_hklE)
      return false;
    nc_assert_always( m_it_hkl->eqv_hkl );
    nc_assert_always( m_it_hkl->multiplicity%2==0 );
    if (m_ii * 2 == m_it_hkl->multiplicity) {
      ++m_it_hkl;
      m_ii = 0;
      return gnp_eh(dspacing,fsq,demi_normal);
    }
    size_t jj( (m_ii++) * 3 );
    fsq = m_it_hkl->fsquared;
    dspacing = m_it_hkl->dspacing;
    demi_normal = m_reci_lattice*Vector(m_it_hkl->eqv_hkl[jj],m_it_hkl->eqv_hkl[jj+1],m_it_hkl->eqv_hkl[jj+2]);
    demi_normal.normalise();
    return true;
  }

  bool PlaneProviderStd::gnp_sg(double& dspacing, double& fsq, Vector& demi_normal)
  {
    if (m_it_hkl == m_it_hklE)
      return false;
    nc_assert(!!m_sg);
    if (m_sg->it == m_sg->itE) {
      if (++m_it_hkl != m_it_hklE) {
        m_sg->prepareLoop(m_it_hkl->h, m_it_hkl->k, m_it_hkl->l,m_it_hkl->multiplicity);
      }
      return gnp_sg(dspacing,fsq,demi_normal);
    }
    fsq = m_it_hkl->fsquared;
    dspacing = m_it_hkl->dspacing;
    demi_normal = m_reci_lattice * Vector(m_sg->it->h,m_sg->it->k,m_sg->it->l);
    demi_normal.normalise();
    ++(m_sg->it);
    return true;
  }

  std::unique_ptr<PlaneProvider> createStdPlaneProvider(shared_obj<const MatInfo> info)
  {
    return std::make_unique<PlaneProviderStd>(std::move(info));
  }
  std::unique_ptr<PlaneProvider> createStdPlaneProvider(const MatInfo* info)
  {
    return std::make_unique<PlaneProviderStd>(info);
  }

}
