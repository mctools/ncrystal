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

#include "NCrystal/internal/NCLCBragg.hh"
#include "NCrystal/internal/NCLCUtils.hh"
#include "NCrystal/internal/NCLCRefModels.hh"
#include "NCrystal/internal/NCSCBragg.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCrystal/internal/NCOrientUtils.hh"
#include "NCrystal/internal/NCPlaneProvider.hh"

namespace NCrystal{

  struct LCBragg::pimpl {

    pimpl(LCBragg * lcbragg, Vector lcaxis, int mode,
          SCOrientation sco, const Info* cinfo, PlaneProvider * plane_provider,
          double mosaicity, double delta_d, double prec,double ntrunc)
      : m_ekin_low(-1)
    {
      nc_assert_always(lcbragg&&cinfo);

      //Convert lcaxis to lab frame:
      RotMatrix reci_lattice = getReciprocalLatticeRot( *cinfo );
      RotMatrix cry2lab = getCrystal2LabRot( sco, reci_lattice );
      Vector lcaxis_labframe = (cry2lab * lcaxis).unit();

      if (mode==0) {
        nc_assert_always(delta_d==0);//mode=0 does not currently support delta_d!=0
        if (!cinfo->hasStructureInfo())
          NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks structure information.");
        nc_assert_always(cinfo&&cinfo->hasStructureInfo());
        const StructureInfo& si = cinfo->getStructureInfo();

        std::unique_ptr<PlaneProvider> stdpp;
        if (!plane_provider) {
          stdpp = createStdPlaneProvider(cinfo);
          plane_provider = stdpp.get();
        }

        m_lchelper = std::make_unique<LCHelper>( lcaxis.unit(),
                                                 lcaxis_labframe,
                                                 mosaicity,
                                                 si.volume * si.n_atoms,
                                                 plane_provider,
                                                 prec, ntrunc);

        m_ekin_low = wl2ekin( m_lchelper->braggThreshold() );

      } else {
        auto scbragg = makeRC<SCBragg>(cinfo,sco,mosaicity,delta_d,plane_provider,prec, ntrunc);
        if (mode>0) {
          m_scmodel = makeRC<LCBraggRef>(scbragg.obj(), lcaxis_labframe, mode);
        } else {
          int nsample = -mode;
          nc_assert(nsample>0);
          m_scmodel = makeRC<LCBraggRndmRot>(scbragg.obj(), lcaxis_labframe, nsample);
        }
        lcbragg->registerSubCalc(m_scmodel.obj());
        nc_assert(m_scmodel->isSubCalc(scbragg.obj()));
        double ekin_high;
        m_scmodel->domain(m_ekin_low,ekin_high);
        nc_assert(ekin_high==kInfinity);
      }

      nc_assert(m_ekin_low>0);
    }

    ~pimpl() = default;

    double m_ekin_low;
    std::unique_ptr<LCHelper> m_lchelper;
    LCHelper::Cache m_lchcache;//TODO: do something more elaborate (or MT
                               //safe?) than this? Caching the last N calls,
                               //might speed up some scenarios.
    RCHolder<Scatter> m_scmodel;
  };

}

NCrystal::LCBragg::LCBragg( const Info* ci, const SCOrientation& sco, double mosaicity,
                            const double (&lcaxis)[3], int mode, double delta_d, PlaneProvider * plane_provider,
                            double prec, double ntrunc)
  : Scatter("LCBragg"),
    m_pimpl(std::make_unique<pimpl>(this,asVect(lcaxis),mode,sco,ci,plane_provider,mosaicity,delta_d,prec,ntrunc))
{
  nc_assert_always(ci);
  nc_assert_always(bool(m_pimpl->m_lchelper)!=bool(m_pimpl->m_scmodel.obj()));
  validate();
}

NCrystal::LCBragg::~LCBragg() = default;

void NCrystal::LCBragg::domain(double& ekin_low, double& ekin_high) const
{
  nc_assert(m_pimpl->m_ekin_low>0);
  ekin_low = m_pimpl->m_ekin_low;
  ekin_high = kInfinity;
}


void NCrystal::LCBragg::generateScattering( double ekin,
                                            const double (&indir)[3],
                                            double (&outdir)[3],
                                            double& delta_ekin ) const
{
  delta_ekin = 0;
  if ( ekin < m_pimpl->m_ekin_low ) {
    asVect(outdir) = asVect(indir);
    return;
  }

  if (! m_pimpl->m_scmodel ) {
    RandomBase* rng = this->getRNG();

    m_pimpl->m_lchelper->genScatter( m_pimpl->m_lchcache, rng, ekin2wl(ekin), asVect(indir), asVect(outdir) );
    return;
  } else {
    m_pimpl->m_scmodel->generateScattering(ekin,indir,outdir,delta_ekin);
  }
}

double NCrystal::LCBragg::crossSection( double ekin,
                                        const double (&indir)[3] ) const
{
  if ( ekin < m_pimpl->m_ekin_low )
    return 0.0;

  if (! m_pimpl->m_scmodel ) {
    return m_pimpl->m_lchelper->crossSection( m_pimpl->m_lchcache, ekin2wl(ekin), asVect(indir) );
  } else {
    return m_pimpl->m_scmodel->crossSection(ekin,indir);
  }
}
