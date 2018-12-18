#ifndef NCrystal_LCUtils_hh
#define NCrystal_LCUtils_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

#include "NCGaussMos.hh"
#include "NCVector.hh"
#include <vector>
#include <cstring>
#include <stdint.h>//cstdint hdr only in C++11

namespace NCrystal {

  //Utilities used to implement layered crystals, intended for internal usage in the LCBragg class.

  class RandomBase;
  class PlaneProvider;

  struct LCPlaneSet {
    //Set of normals with same angle to lcaxis and d-spacing.

    LCPlaneSet(double dspacing, double alpha, double truncangle, double fsquared = 0.0);
    //Add contribution to LCPlaneSet:
    void addFsq(double fsquared) { fsq += fsquared; }
    double twodsp;// 2*dspacing
    double inv_twodsp;// 1/twodsp
    double cosalpha;
    double sinalpha;
    double cosalphaminus;//cos(max[0,alpha-truncangle])
    double cosalphaplus;//cos(alpha+truncangle)
    double fsq;//Total FSquared of all normals in planeset.
    bool isOnAxis() const { return !sinalpha; }

    ~LCPlaneSet(){}
  };

  class LCROI {
  public:
    //Region of interest. In the non-degenerate case this is a phi-rotation
    //range (subset of [0,pi]) and associated PlaneSet for which a given neutron
    //has non-zero (and supposedly smoothly varying) cross-sections. This is
    //given in the standard frame, used by the LCStdFrame class.

    double rotmin;
    double rotmax;
    const LCPlaneSet * planeset;//ptr must be valid for lifetime of LCROI object!
    double normal_sign;//1.0 for normal(s), -1.0 for anti-normal(s).

    //Standard non-degenerate cases:
    LCROI(double rmin, double rmax,const LCPlaneSet* ps, double normsign);
    ~LCROI(){}
    bool contains(double t) const;
    double length() const;

    //degenerate cases, when either neutron or normal is on lcaxis:
    LCROI(const LCPlaneSet* ps, double normsign);
    bool isDegenerate() const;
    bool normalIsOnAxis() const;
    bool neutronIsOnAxis() const;
  };

  class LCROIFinder {
    //Class able to search for ROI's for a given neutron, based on mosaicity
    //bands in the provided planesets. If found by findROIs, they will be
    //appended to roilist. It is slightly faster if subsequent calls to findROIs
    //are ordered by d-spacing (i.e. it is optimised for the case where
    //planesets are kept sorted by d-spacing and findROIs are called for all of
    //them in turn. Note that the search happens in the standard frame, used by
    //the LCStdFrame class.

  public:
    //NB: Results are symmetric between neutron_angle_to_lcaxis -> pi - neutron_angle_to_lcaxis.
    LCROIFinder(double wavelength, double cos_neutron_angle_to_lcaxis,
                double cos_truncation_angle, double sin_truncation_angle);
    ~LCROIFinder();
    void findROIs(const LCPlaneSet * planeset,std::vector<LCROI>& roilist);
  private:
    const double m_wl;
    const double m_c3;
    const double m_s3;
    const double m_cta;
    const double m_sta;
    double m_prev2d;
    double m_c2;
    double m_c23;
    double m_s2;
    double m_s23;
    bool m_s2approximated;
  };

  class LCStdFrame {
    //Helper class use to evaluate cross-sections and scatterings in layered
    //crystals, using an idealised frame where lcaxis=(0,0,1) and the neutron
    //direction is (s3,0,c3). Internally, the GaussMos helper class is used, but
    //wrapped appropriately. Note that this class does NOT apply any truncation
    //cut, which is assumed to take place in the calling code.
  public:
    //Constructor takes same parameters as GaussMos.
    LCStdFrame(double mosaicity, bool mosaicity_is_fhwm = true, double prec = 1e-3, double ntrunc = 0.0 );
    ~LCStdFrame();

    //For reference, we provide const access to underlying GaussMos object
    const GaussMos& gaussMos() const { return m_gm; }

    struct NeutronPars {
      NeutronPars(double wl, double c3, double s3);
      ~NeutronPars(){}
      const double wl;
      const double c3;
      const double s3;
    };
    struct NormalPars {
      NormalPars( const LCPlaneSet* ps, double normal_sign );
      ~NormalPars(){}
      const LCPlaneSet* planeset;
      const double sign;
    };

    //When the normal is on-axis (e.g. parallel to lcaxis), the following two special methods can be used:
    double calcXS_OnAxis( const NeutronPars&, const NormalPars&) const;
    void genScat_OnAxis(RandomBase * rand, const NeutronPars&, const NormalPars&, Vector& outdir ) const;

    //For off-axis normals, one must consider the particular rotation of the
    //crystallite. We define phi so that at phi=0, the normal is given by
    //(sinthn,0,costhn), where thn is the angle between the normals in the
    //planeset and the lcaxis. For phi!=0, this normal is then rotated around
    //lcaxis=(0,0,1), yielding the normal (sinthn*cosphi,sinthn*sinphi,costhn).
    static Vector normalInStdFrame( const NormalPars&, double cosphi, double sinphi );

    //The cross-section for off-axis normals can either be found for a
    //particular crystallite orientation...:
    double calcXS( const NeutronPars&, const NormalPars&, double cosphi ) const;
    //... or it can be integrated over a range of crystallite orientations (the
    //collect object will be passed to the internal AdaptiveSimpsons
    //integrator):
    double calcXSIntegral( const NeutronPars&, const NormalPars&, double phimin, double phimax ) const;

    //Scatterings (sinphi parameter required as well as cosphi, to allow calling
    //code to select phi values in all of [-pi,pi] range):
    void genScat( RandomBase * rand, const NeutronPars&, const NormalPars&,
                  double cosphi, double sinphi, Vector& outdir ) const;

  private:
    const GaussMos m_gm;
  };

  class LCHelper {
    //Class which can provide cross-sections and scatterings for planes with
    //normals not parallel to the lcaxis. Prec and ntrunc parameters will be
    //passed on directly to the internal GaussMos object.
  public:
    LCHelper( Vector lcaxis_crystalframe,
              Vector lcaxis_labframe,
              double mosaicity_fwhm,
              double unitcell_volume_times_natoms,
              PlaneProvider * pp,
              double prec = 1e-3,
              double ntrunc = 0.0 );

    //Usage happens via Cache objects (allowing users of the class to decide
    //upon caching strategies themselves). One should not share Cache objects
    //between different LCHelper instances, except if the cache is first reset
    //via a call to Cache::reset().
    class Cache;
    bool isValid(Cache&, double wavelength, double c3 ) const;
    bool isValid(Cache&, double wavelength, const Vector& indir) const;
    void ensureValid(Cache&, double wavelength, const Vector& indir) const;

    //Valid caches can be used to get cross-sections or generate scatterings:
    double crossSection( Cache&, double wavelength, const Vector& indir ) const;
    void genScatter( Cache&, RandomBase*, double wavelength, const Vector& indir, Vector& outdir ) const;

    //Access without cache.
    void genScatterNoCache( RandomBase*, double wavelength, const Vector& indir, Vector& outdir ) const;
    double crossSectionNoCache( double wavelength, const Vector& indir ) const;

    double braggThreshold() const;//max wavelength, beyond which all cross-sections will be 0.


    //Mechanics:
    ~LCHelper();
#if __cplusplus >= 201103L
    LCHelper(const LCHelper&) = delete;
    void operator=(const LCHelper&) = delete;
#endif

    const GaussMos& gaussMos() { return m_lcstdframe.gaussMos(); }

  private:
    friend class Cache;
    Vector m_lcaxislab;
    std::vector<LCPlaneSet> m_planes;//sorted by dspacing, largest first.
    LCStdFrame m_lcstdframe;
    double m_xsfact;
    void forceUpdateCache( Cache&, uint64_t discr_wl, uint64_t discr_c3 ) const;
    struct Overlay {
      static const unsigned ndata = 8;
      Overlay();
      ~Overlay();
      void clear();
      void prepareNullArray();
      float * data;
      double nonCommulVal(unsigned i) const;
#if __cplusplus >= 201103L
      Overlay(Overlay&& o);
      Overlay& operator=(Overlay&& o);
      Overlay(const Overlay& o) = delete;
      Overlay& operator=(const Overlay& o) = delete;
#else
      Overlay(const Overlay& o);
      Overlay& operator=(const Overlay& o);
#endif
    };
  static void genPhiVal(RandomBase* rand, const LCROI& roi, const Overlay& overlay, double& phi, double& overlay_at_phi);

  public:
    class Cache {
    public:
      Cache();//default constructs invalid cache object
      ~Cache();
      void reset();//invalidates cache, but retains dynamically acquired memory for later reuse.
    private:
      friend class LCHelper;
      std::pair<uint64_t,uint64_t> m_signature;//discretised (wavelength,c3)
      double m_wl;//<-- Neutron wavelength. Actually dediscretized m_signature.first
      double m_c3;//<-- dot(indir,lcaxis). Actually dediscretized m_signature.second.
      double m_s3;//sqrt(1-m_c3*m_c3)
      std::vector<LCROI> m_roilist;
      std::vector<double> m_roixs_commul;//for selecting
      std::vector<Overlay> m_roi_overlays;//for selecting
    };
  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {

  inline LCROI::LCROI(double rmin, double rmax,const LCPlaneSet* ps, double normsign)
    : rotmin(rmin), rotmax(rmax), planeset(ps), normal_sign(normsign)
  {
    nc_assert(ps);
    nc_assert(normsign==1.||normsign==-1.);
    nc_assert(rmax>rmin);
    nc_assert(rmin>=0.0&&rmax<=kPi);
  }

  inline LCROI::LCROI(const LCPlaneSet* ps, double normsign )
    : rotmin(ps->isOnAxis()?0:kPi), rotmax(ps->isOnAxis()?0:kPi), planeset(ps), normal_sign(normsign)
  {
    nc_assert(ps);
    nc_assert(normsign==1.||normsign==-1.);
    nc_assert(rotmin==rotmax);
  }
  inline bool LCROI::normalIsOnAxis() const { nc_assert(planeset->isOnAxis()==(rotmax==0.0)); return rotmax==0.0; }
  inline bool LCROI::neutronIsOnAxis() const { nc_assert(!normalIsOnAxis());/*if both on axis, classify as normal-is-on-axis*/ return rotmin==kPi; }
  inline bool LCROI::isDegenerate() const { return rotmin==rotmax; }
  inline double LCROI::length() const { nc_assert(!isDegenerate()); return rotmax-rotmin; }
  inline bool LCROI::contains(double t) const { nc_assert(!isDegenerate()); return t>=rotmin && t<= rotmax; }


  inline LCStdFrame::NeutronPars::NeutronPars(double thewl, double thec3, double thes3) : wl(thewl), c3(thec3), s3(thes3) {}
  inline LCStdFrame::NormalPars::NormalPars( const LCPlaneSet* ps, double normal_sign )
    : planeset(ps),
      sign(normal_sign)
  {
    nc_assert( normal_sign==1. || normal_sign==-1. );
  }

  inline LCHelper::Cache::Cache() : m_signature(std::numeric_limits<uint64_t>::max(),
                                                std::numeric_limits<uint64_t>::max()),
                                    m_wl(-99), m_c3(-99), m_s3(-99)
  {
    //Starts in same state as after calling Cache::reset()
  }
  inline LCHelper::Cache::~Cache(){}

  inline LCHelper::Overlay::Overlay() : data(0) {}
  inline LCHelper::Overlay::~Overlay() { delete[] data; }
  inline void LCHelper::Overlay::clear() { delete[] data; data = 0; }
  inline void LCHelper::Overlay::prepareNullArray() { if (!data) { data = new float[ndata]; } std::memset(data,0,sizeof(float)*ndata); }
  inline double LCHelper::Overlay::nonCommulVal(unsigned i) const { nc_assert(i<ndata); return i ? data[i]-(double)data[i-1] : (double)data[i]; }
#if __cplusplus >= 201103L
  inline LCHelper::Overlay::Overlay(LCHelper::Overlay&& o) { std::swap(data,o.data); }
  inline LCHelper::Overlay& LCHelper::Overlay::operator=(LCHelper::Overlay&& o) { std::swap(data,o.data); return *this; }
#else
  inline LCHelper::Overlay::Overlay(const LCHelper::Overlay& o) : data(0) { *this = o; }
  inline LCHelper::Overlay& LCHelper::Overlay::operator=(const LCHelper::Overlay& o)
  {
    if (o.data) {
      if (!data)
        data = new float[ndata];
      std::memcpy(data,o.data,ndata*sizeof(float));
    } else {
      clear();
    }
    return *this;
  }
#endif

}

#endif

