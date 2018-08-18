#ifndef NCrystal_Info_hh
#define NCrystal_Info_hh

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

#include "NCrystal/NCDefs.hh"
#include <vector>
#include <stdint.h>//cstdint hdr only in C++11


//////////////////////////////////////////////////////////////////////
// Data class containing information (high level or derived) about  //
// a given crystal. Instances of the class are typically generated  //
// by dedicated factories, based on interpretation of crystallogra- //
// phic data files. Physics models in the form of for example       //
// NCScatter or NCAbsorption instances, are then initialised from   //
// these Info objects, thus providing a useful separation between   //
// data sources and algorithms working on the data.                 //
//////////////////////////////////////////////////////////////////////

namespace NCrystal {

  //Data structures:

  struct NCRYSTAL_API StructureInfo {
    unsigned spacegroup;//From 1-230 if provided, 0 if information not available
    double lattice_a;//angstrom
    double lattice_b;//angstrom
    double lattice_c;//angstrom
    double alpha;//degree
    double beta;//degree
    double gamma;//degree
    double volume;//Aa^3
    unsigned n_atoms;//Number of atoms per unit cell
  };

  struct NCRYSTAL_API HKLInfo {
    double dspacing;//angstrom
    double fsquared;//barn
    int h;
    int k;
    int l;
    unsigned multiplicity;

    //If the HKLInfo source knows the plane normals, they will be provided in
    //the following list as unit vectors. Only half of the normals should be
    //included in this list, since if n is a normal, so is -n. If demi_normals
    //is not empty, it will be true that multiplicity == 2*demi_normals.size().
    struct Normal {
      Normal(double a1, double a2, double a3) : x(a1), y(a2), z(a3) {}
      double x, y, z;
    };
    std::vector<Normal> demi_normals;//TODO for NC2: vector->pointer saves 16B

    //If eqv_hkl is not a null pointer, it contains the corresponding Miller
    //indices of the demi_normals as three 2-byte integers (short). Thus,
    //eqv_hkl has demi_normal.size()*3 entries:
    short * eqv_hkl;
    HKLInfo();
    ~HKLInfo();
    HKLInfo(const HKLInfo &o);
    HKLInfo& operator=(const HKLInfo &o);
#if __cplusplus >= 201103L
    HKLInfo(const HKLInfo &&o);
    HKLInfo& operator=(HKLInfo &&o);
#endif
  };

  typedef std::vector<HKLInfo> HKLList;

  struct NCRYSTAL_API AtomInfo {
    //TODO for NC2: More parameters capable of handling non-natural atoms
    AtomInfo() : atomic_number(0), number_per_unit_cell(0),debye_temp(0.),mean_square_displacement(0.) {}
    unsigned atomic_number;
    unsigned number_per_unit_cell;
    //per-element debye temperature (0.0 if not available, see hasPerElementDebyeTemperature() below):
    double debye_temp;
    //Atomic coordinates (vector must be empty or have number_per_unit_cell
    //entries):
    struct Pos { Pos(double a, double b, double c) : x(a),y(b),z(c) {}; double x, y, z; };
    std::vector<Pos> positions;
    //Mean-square-displacements in angstrom (0.0 if not available):
    double mean_square_displacement;
  };

  typedef std::vector<AtomInfo> AtomList;

  struct NCRYSTAL_API XSectProvider {
    //Provide non-Bragg scattering cross sections.
    //Accept lambda in angstrom and return x-sect in barn
    virtual double xsectScatNonBragg(const double& lambda) const = 0;
    virtual ~XSectProvider(){}
  };

  class NCRYSTAL_API Info : public RCBase {
  public:

    /////////////////////////////////////////
    // Information about crystal structure //
    /////////////////////////////////////////

    bool hasStructureInfo() const;
    const StructureInfo& getStructureInfo() const;


    /////////////////////////////////////////////
    // Information about cross-sections [barn] //
    /////////////////////////////////////////////

    //absorption cross-section (at 2200m/s):
    bool hasXSectAbsorption() const;
    double getXSectAbsorption() const;

    //saturated scattering cross-section (high E limit):
    bool hasXSectFree() const;
    double getXSectFree() const;

    /////////////////////////////////////////////////////////////////////////////////
    // Provides calculation of "background" (non-Bragg diffraction) cross-sections //
    /////////////////////////////////////////////////////////////////////////////////

    bool providesNonBraggXSects() const;
    double xsectScatNonBragg(double lambda) const;

    ///////////////////////////
    // Temperature [kelvin]  //
    ///////////////////////////

    bool hasTemperature() const;
    double getTemperature() const;

    /////////////////////////////////
    // Debye temperature [kelvin]  //
    /////////////////////////////////

    //Global Debye temperature:
    bool hasDebyeTemperature() const;
    double getDebyeTemperature() const;

    //Whether AtomInfo objects have per-element Debye temperatures available:
    bool hasPerElementDebyeTemperature() const;

    //////////////////////
    // Atom Information //
    //////////////////////

    bool hasAtomInfo() const;
    AtomList::const_iterator atomInfoBegin() const;
    AtomList::const_iterator atomInfoEnd() const;

    //Whether AtomInfo objects have atomic coordinates available:
    bool hasAtomPositions() const;

    //Whether AtomInfo objects have mean-square-displacements available:
    bool hasAtomMSD() const;

    //See also hasPerElementDebyeTemperature() above.

    /////////////////////
    // HKL Information //
    /////////////////////

    bool hasHKLInfo() const;
    unsigned nHKL() const;
    HKLList::const_iterator hklBegin() const;//first (==end if empty)
    HKLList::const_iterator hklLast() const;//last (==end if empty)
    HKLList::const_iterator hklEnd() const;
    //The limits:
    double hklDLower() const;
    double hklDUpper() const;
    //The largest/smallest (both returns inf if nHKL=0):
    double hklDMinVal() const;
    double hklDMaxVal() const;

    //////////////////////////////
    // Expanded HKL Information //
    //////////////////////////////

    //Whether HKLInfo objects have demi_normals available:
    bool hasHKLDemiNormals() const;

    //Whether HKLInfo objects have eqv_hkl available:
    bool hasExpandedHKLInfo() const;

    //Search eqv_hkl lists for specific (h,k,l) value. Returns hklEnd() if not found:
    HKLList::const_iterator searchExpandedHKL(short h, short k, short l) const;

    /////////////////////
    // Density [g/cm^3] //
    /////////////////////

    bool hasDensity() const;
    double getDensity() const;

    //////////////////////////////
    // Internals follow here... //
    //////////////////////////////

  public:
    //Methods used by factories when setting up an Info object:
    Info();
    void addAtom(const AtomInfo& ai) {ensureNoLock(); m_atomlist.push_back(ai); };
    void enableHKLInfo(double dlower, double dupper);
    void addHKL(const HKLInfo& hi) { ensureNoLock(); m_hkllist.push_back(hi); }
    void setStructInfo(const StructureInfo& si) { ensureNoLock(); nc_assert_always(si.spacegroup!=999999); m_structinfo = si; }
    void setXSectFree(double x) { ensureNoLock(); m_xsect_free = x; }
    void setXSectAbsorption(double x) { ensureNoLock(); m_xsect_absorption = x; }
    void setTemperature(double t) { ensureNoLock(); m_temp = t; }
    void setDebyeTemperature(double dt) { ensureNoLock(); m_debyetemp = dt; }
    void setDensity(double d) { ensureNoLock(); m_density = d; }
    void setXSectProvider(XSectProvider*xp) { ensureNoLock(); m_xsectprovider = xp; }//assumes ownership
    void objectDone();//Finish up (sorts hkl list (by dspacing first), and atom info list (by Z first)). This locks the instance.
    bool isLocked() const { return m_lock; }

    //Unique ID for this Info instance, useful for downstream caching purposes:
    uint64_t getUniqueID() const { nc_assert(m_lock); return m_uniqueid; }

  private:
    void ensureNoLock();
    uint64_t m_uniqueid;
    StructureInfo m_structinfo;
    AtomList m_atomlist;
    HKLList m_hkllist;//sorted by dspacing first
    double m_hkl_dlower;
    double m_hkl_dupper;
    double m_density;
    double m_xsect_free;
    double m_xsect_absorption;
    double m_temp;
    double m_debyetemp;
    XSectProvider * m_xsectprovider;
    bool m_lock;
  protected:
    virtual ~Info();
  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {
  inline bool Info::hasStructureInfo() const { return m_structinfo.spacegroup!=999999; }
  inline const StructureInfo& Info::getStructureInfo() const  { nc_assert(hasStructureInfo()); return m_structinfo; }
  inline bool Info::hasXSectAbsorption() const { return m_xsect_absorption >= 0.0; }
  inline bool Info::hasXSectFree() const { return m_xsect_free >= 0.0; }
  inline double Info::getXSectAbsorption() const { nc_assert(hasXSectAbsorption()); return m_xsect_absorption; }
  inline double Info::getXSectFree() const { nc_assert(hasXSectFree()); return m_xsect_free; }
  inline bool Info::providesNonBraggXSects() const { return m_xsectprovider!=0; }
  inline double Info::xsectScatNonBragg(double lambda) const  { return m_xsectprovider->xsectScatNonBragg(lambda); }
  inline bool Info::hasTemperature() const { return m_temp > 0.0; }
  inline bool Info::hasDebyeTemperature() const { return m_debyetemp > 0.0; }
  inline double Info::getTemperature() const { nc_assert(hasTemperature()); return m_temp; }
  inline double Info::getDebyeTemperature() const { nc_assert(hasDebyeTemperature()); return m_debyetemp; }
  inline bool Info::hasPerElementDebyeTemperature() const { return hasAtomInfo() && m_atomlist.front().debye_temp > 0.0; }
  inline bool Info::hasAtomPositions() const { return hasAtomInfo() && !m_atomlist.front().positions.empty(); }
  inline bool Info::hasAtomMSD() const { return hasAtomInfo() && m_atomlist.front().mean_square_displacement>0.0; }
  inline bool Info::hasAtomInfo() const  { return !m_atomlist.empty(); }
  inline AtomList::const_iterator Info::atomInfoBegin() const { nc_assert(hasAtomInfo()); return m_atomlist.begin(); }
  inline AtomList::const_iterator Info::atomInfoEnd() const { nc_assert(hasAtomInfo()); return m_atomlist.end(); }
  inline bool Info::hasHKLInfo() const { return m_hkl_dupper>=m_hkl_dlower; }
  inline bool Info::hasExpandedHKLInfo() const { return hasHKLInfo() && !m_hkllist.empty() && m_hkllist.front().eqv_hkl; }
  inline bool Info::hasHKLDemiNormals() const { return hasHKLInfo() && !m_hkllist.empty() && ! m_hkllist.front().demi_normals.empty(); }
  inline unsigned Info::nHKL() const { nc_assert(hasHKLInfo()); return m_hkllist.size(); }
  inline HKLList::const_iterator Info::hklBegin() const { nc_assert(hasHKLInfo()); return m_hkllist.begin(); }
  inline HKLList::const_iterator Info::hklLast() const { nc_assert(hasHKLInfo()); std::size_t nhkl = m_hkllist.size(); return nhkl ? m_hkllist.begin()+(nhkl-1) : m_hkllist.end(); }
  inline HKLList::const_iterator Info::hklEnd() const { nc_assert(hasHKLInfo()); return m_hkllist.end(); }
  inline double Info::hklDLower() const { nc_assert(hasHKLInfo()); return m_hkl_dlower; }
  inline double Info::hklDUpper() const { nc_assert(hasHKLInfo()); return m_hkl_dupper; }
  inline bool Info::hasDensity() const { return m_density > 0.0; }
  inline double Info::getDensity() const { nc_assert(hasDensity()); return m_density; }
  inline HKLInfo::HKLInfo() : dspacing(0.), fsquared(0.), h(0), k(0), l(0),  multiplicity(0), eqv_hkl(0) {}
  inline HKLInfo::~HKLInfo() { delete[] eqv_hkl; }
  inline HKLInfo::HKLInfo(const HKLInfo &o) : eqv_hkl(0) { *this = o; }
}

#endif
