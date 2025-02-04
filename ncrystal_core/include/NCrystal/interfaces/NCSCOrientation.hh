#ifndef NCrystal_SCOrientation_hh
#define NCrystal_SCOrientation_hh

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

#include "NCrystal/core/NCVariant.hh"
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  class NCRYSTAL_API SCOrientation {
  public:

    //Class which defines the orientation of a given single-crystal in the
    //laboratory frame, and which is needed to construct an SCBragg
    //instance.

    static constexpr const double default_tolerance = 1e-4;

    //Construct unspecified orientation, which must subsequently be completed
    //with calls to both setPrimaryDirection and setSecondaryDirection:
    SCOrientation() = default;

    //Provide vectors in the crystal frame and specify their desired direction
    //in the laboratory frame. This is well-formed as long as the opening angle
    //between primary and secondary vectors is similar in the two frames within
    //the stated tolerance in radians, because only the component of each
    //secondary direction ortogonal to the corresponding primary direction will
    //be used. The vectors do not need to be normalised:
    void setPrimaryDirection( const CrystalAxis&, const LabAxis& );
    void setSecondaryDirection( const CrystalAxis&, const LabAxis&,
                                double tolerance = default_tolerance );

    //Often, however, it is more convenient and useful to specify the directions
    //in the crystal frame by simply indicating a point in hkl space, and
    //letting the corresponding direction be the normal of the hkl-plane in
    //question, defined as:
    //
    // a_reci*h + b_reci*k + c_reci*l
    //
    //with (a_reci, b_reci, c_reci) being the reciprocal lattice indices. To
    //define directions in this fashion, use the following methods:
    void setPrimaryDirection( const HKLPoint&, const LabAxis& );
    void setSecondaryDirection( const HKLPoint&, const LabAxis&,
                                double tolerance = default_tolerance );

    //Alternatively, set with OrientDir objects which already encodes the
    //CrystalAxis vs. HKLPoint concept:
    void setPrimaryDirection( const OrientDir& );
    void setSecondaryDirection( const OrientDir&, double tolerance = default_tolerance );

    //Complete when both primary and secondary directions have been set:
    bool isComplete() const;

    //Access state of complete object (throws LogicError if not isComplete()):
    struct Data {
      OrientDir dir1, dir2;
      double dirtol;
    };
    Data getData() const;

    //Miscellaneous:
    void clear();
    void stream(std::ostream&) const;
    bool operator==( const SCOrientation& ) const;

  private:
    Optional<OrientDir> m_dir1data;
    Optional<std::pair<OrientDir,double>> m_dir2data;
  };

  std::ostream& operator<<( std::ostream&, const SCOrientation& );

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline std::ostream& operator<<(std::ostream& os, const SCOrientation& sco )
  {
    sco.stream(os);
    return os;
  }

  inline void SCOrientation::setPrimaryDirection( const CrystalAxis& c, const LabAxis& l )
  {
    setPrimaryDirection( OrientDir{c,l} );
  }

  inline void SCOrientation::setSecondaryDirection( const CrystalAxis& c, const LabAxis& l, double t )
  {
    setSecondaryDirection( OrientDir{c,l}, t );
  }

  inline void SCOrientation::setPrimaryDirection( const HKLPoint& c, const LabAxis& l )
  {
    setPrimaryDirection( OrientDir{c,l} );
  }

  inline void SCOrientation::setSecondaryDirection( const HKLPoint& c, const LabAxis& l, double t )
  {
    setSecondaryDirection( OrientDir{c,l}, t );
  }

  inline bool SCOrientation::isComplete() const
  {
    return m_dir1data.has_value() && m_dir2data.has_value();
  }

  inline void SCOrientation::clear()
  {
    m_dir1data = NullOpt;
    m_dir2data = NullOpt;
  }

  inline bool SCOrientation::operator==( const SCOrientation& o ) const
  {
    return m_dir1data == o.m_dir1data && m_dir2data == o.m_dir2data;
  }

  inline SCOrientation::Data SCOrientation::getData() const
  {
    if (!isComplete())
      NCRYSTAL_THROW(LogicError,"Incomplete SCOrientation object - must set both primary and secondary directions.");
    return { m_dir1data.value(), m_dir2data.value().first, m_dir2data.value().second };
  }

}

#endif
