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

#include "NCrystal/internal/lcbragg/NCLCBragg.hh"
#include "NCrystal/internal/extd_utils/NCLCUtils.hh"
#include "NCrystal/internal/extd_utils/NCLCRefModels.hh"
#include "NCrystal/internal/scbragg/NCSCBragg.hh"
#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  struct LCBragg::pimpl {

    pimpl(LCBragg * lcbragg, LCAxis lcaxis, int mode,
          SCOrientation sco, const Info& cinfo, PlaneProvider * plane_provider,
          MosaicityFWHM mosaicity, double delta_d, double prec,double ntrunc)
      : m_ekin_low(-1)
    {
      nc_assert_always(lcbragg);
      if (!cinfo.hasStructureInfo())
        NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks structure information.");

      //Convert lcaxis to lab frame:
      const StructureInfo& si = cinfo.getStructureInfo();
      RotMatrix reci_lattice = getReciprocalLatticeRot( si );
      RotMatrix cry2lab = getCrystal2LabRot( sco, reci_lattice );
      LCAxis lcaxis_labframe = (cry2lab * lcaxis.as<Vector>()).unit().as<LCAxis>();

      if (mode==0) {
        nc_assert_always(delta_d==0);//mode=0 does not currently support delta_d!=0

        std::unique_ptr<PlaneProvider> stdpp;
        if (!plane_provider) {
          stdpp = createStdPlaneProvider(&cinfo);
          plane_provider = stdpp.get();
        }

        m_lchelper = std::make_unique<LCHelper>( lcaxis.as<Vector>().unit().as<LCAxis>(),
                                                 lcaxis_labframe,
                                                 mosaicity,
                                                 si.volume * si.n_atoms,
                                                 plane_provider,
                                                 prec, ntrunc);

        m_ekin_low = wl2ekin( m_lchelper->braggThreshold() );

      } else {
        auto scbragg = makeSO<SCBragg>(cinfo,sco,mosaicity,delta_d,plane_provider,prec, ntrunc);
        if (mode>0) {
          m_scmodel = std::make_shared<LCBraggRef>(scbragg, lcaxis_labframe, mode);
        } else {
          int nsample = -mode;
          nc_assert(nsample>0);
          m_scmodel = std::make_shared<LCBraggRndmRot>(scbragg, lcaxis_labframe, nsample);
        }
        m_ekin_low = m_scmodel->domain().elow.get();
        nc_assert(ncisinf(m_scmodel->domain().ehigh.get()));
      }

      nc_assert(m_ekin_low>0);
    }

    ~pimpl() = default;

    double m_ekin_low;
    std::unique_ptr<LCHelper> m_lchelper;
    ProcImpl::OptionalProcPtr m_scmodel;
  };

}

NC::LCBragg::LCBragg( const Info& ci, const SCOrientation& sco, MosaicityFWHM mosaicity,
                      const LCAxis& lcaxis, int mode, double delta_d, PlaneProvider * plane_provider,
                      double prec, double ntrunc)
  : m_pimpl(std::make_unique<pimpl>(this,lcaxis,mode,sco,ci,plane_provider,mosaicity,delta_d,prec,ntrunc))
{
  nc_assert_always(bool(m_pimpl->m_lchelper)!=bool(m_pimpl->m_scmodel!=nullptr));
}

NC::LCBragg::~LCBragg() = default;

NC::EnergyDomain NC::LCBragg::domain() const noexcept
{
  return { NeutronEnergy{m_pimpl->m_ekin_low}, NeutronEnergy{kInfinity} };
}

NC::CrossSect NC::LCBragg::crossSection(NC::CachePtr& cp, NC::NeutronEnergy ekin, const NC::NeutronDirection& indir ) const
{
  if ( ekin.get() < m_pimpl->m_ekin_low )
    return CrossSect{ 0.0 };

  if (! m_pimpl->m_scmodel ) {
    NeutronWavelength wl{ekin};
    if (!(wl.get()>0.0))
      return CrossSect{ 0.0 };
    const Vector& indirv = indir.as<Vector>().unit();
    return CrossSect{ m_pimpl->m_lchelper->crossSection( accessCache<LCHelper::Cache>(cp), wl.get(), indirv ) };
  } else {
    return CrossSect{ m_pimpl->m_scmodel->crossSection( cp, ekin, indir ) };
  }
}

NC::ScatterOutcome NC::LCBragg::sampleScatter(NC::CachePtr& cp, NC::RNG& rng, NC::NeutronEnergy ekin, const NC::NeutronDirection& indir ) const
{
  if ( ekin.get() < m_pimpl->m_ekin_low )
    return { ekin, indir };

  if (! m_pimpl->m_scmodel ) {
    NeutronWavelength wl{ekin};
    if (!(wl.get()>0.0))
      return { ekin, indir };

    const Vector& indirv = indir.as<Vector>().unit();
    Vector outdir;
    m_pimpl->m_lchelper->genScatter( accessCache<LCHelper::Cache>(cp), rng, wl.get(), indirv, outdir );
    return { ekin, outdir.as<NeutronDirection>() };
  } else {
    return m_pimpl->m_scmodel->sampleScatter( cp, rng, ekin, indir );
  }
}
