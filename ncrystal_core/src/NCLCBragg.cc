////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCrystal/NCLCBragg.hh"
#include "NCLCUtils.hh"
#include "NCLCRefModels.hh"
#include "NCrystal/NCSCBragg.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCDefs.hh"
#include "NCVector.hh"
#include "NCLatticeUtils.hh"
#include "NCOrientUtils.hh"
#include "NCPlaneProvider.hh"

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

        UniquePtr<PlaneProvider> stdpp;
        if (!plane_provider) {
          stdpp = createStdPlaneProvider(cinfo);
          plane_provider = stdpp.obj();
        }

        m_lchelper = new LCHelper( lcaxis.unit(),
                                   lcaxis_labframe,
                                   mosaicity,
                                   si.volume * si.n_atoms,
                                   plane_provider,
                                   prec, ntrunc);


        m_ekin_low = wl2ekin( m_lchelper.obj()->braggThreshold() );

      } else {
        RCHolder<Scatter> scbragg(new SCBragg(cinfo,sco,mosaicity,delta_d,plane_provider,prec, ntrunc));
        if (mode>0) {
          m_scmodel = new LCBraggRef(scbragg.obj(), lcaxis_labframe, mode);
        } else {
          int nsample = -mode;
          nc_assert(nsample>0);
          m_scmodel = new LCBraggRndmRot(scbragg.obj(), lcaxis_labframe, nsample);
        }
        lcbragg->registerSubCalc(m_scmodel.obj());
        nc_assert(m_scmodel.obj()->isSubCalc(scbragg.obj()));
        double ekin_high;
        m_scmodel.obj()->domain(m_ekin_low,ekin_high);
        nc_assert(ekin_high==infinity);
      }

      nc_assert(m_ekin_low>0);
    }

    ~pimpl()
    {
    }

    double m_ekin_low;
    UniquePtr<LCHelper> m_lchelper;
    LCHelper::Cache m_lchcache;//TODO for NC2: do something more elaborate (or MT
                               //safe?) than this? Caching the last N calls,
                               //might speed up some scenarios.
    RCHolder<Scatter> m_scmodel;
  };

}

NCrystal::LCBragg::LCBragg( const Info* ci, const SCOrientation& sco, double mosaicity,
                            const double (&lcaxis)[3], int mode, double delta_d, PlaneProvider * plane_provider,
                            double prec, double ntrunc)
  : Scatter("LCBragg"),
    m_pimpl(new pimpl(this,asVect(lcaxis),mode,sco,ci,plane_provider,mosaicity,delta_d,prec,ntrunc))
{
  nc_assert_always(ci);
  nc_assert_always(bool(m_pimpl->m_lchelper())!=bool(m_pimpl->m_scmodel()));
  validate();
}

NCrystal::LCBragg::~LCBragg()
{
  delete m_pimpl;
}

void NCrystal::LCBragg::domain(double& ekin_low, double& ekin_high) const
{
  nc_assert(m_pimpl->m_ekin_low>0);
  ekin_low = m_pimpl->m_ekin_low;
  ekin_high = infinity;
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

    m_pimpl->m_lchelper.obj()->genScatter( m_pimpl->m_lchcache, rng, ekin2wl(ekin), asVect(indir), asVect(outdir) );
    return;
  } else {
    m_pimpl->m_scmodel.obj()->generateScattering(ekin,indir,outdir,delta_ekin);
  }
}

double NCrystal::LCBragg::crossSection( double ekin,
                                        const double (&indir)[3] ) const
{
  if ( ekin < m_pimpl->m_ekin_low )
    return 0.0;

  if (! m_pimpl->m_scmodel ) {
    return m_pimpl->m_lchelper.obj()->crossSection( m_pimpl->m_lchcache, ekin2wl(ekin), asVect(indir) );
  } else {
    return m_pimpl->m_scmodel.obj()->crossSection(ekin,indir);
  }
}
