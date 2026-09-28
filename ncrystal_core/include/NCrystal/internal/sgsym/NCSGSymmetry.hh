#ifndef NCrystal_SGSymmetry_hh
#define NCrystal_SGSymmetry_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"

namespace NCRYSTAL_NAMESPACE {

  //Unit cell parameters (lengths in Aa, angles in degrees):
  struct CellParameters { double a, b, c, alpha, beta, gamma; };

  //Constraints on the unit cell parameters imposed by the symmetry of a space
  //group setting (i.e. the cell metric must be invariant under all rotation
  //parts of the operations), which are derived from the operations. The
  //length a and the angle alpha are never constrained to other parameters.
  //
  //TODO: Replace checkAndCompleteLattice[Angles]() in NCLatticeUtils.hh with
  //this.
  struct SGCellConstraints {
    enum class LengthConstraint { Free, EqualToA };
    enum class AngleConstraint { Free, Is90, Is120, EqualToAlpha };
    LengthConstraint b, c;
    AngleConstraint alpha, beta, gamma;

    //Check fully specified parameters, with exact equality required for
    //constrained values. Also requires positive lengths, angles in (0,180),
    //and alpha<120 if all angles must be equal. Throws BadInput if invalid
    //(including for parameters which are 0, cf. complete(..) below):
    void check( const CellParameters& ) const;

    //Fill in any parameters given as 0 whose values are implied by the
    //constraints (e.g. b and c for cubic space groups, or gamma=120 for
    //hexagonal axes), and then check(..) the result:
    void complete( CellParameters& ) const;
  };

  //Output like "a, b=a, c, alpha=90, beta, gamma=90":
  std::ostream& operator<<( std::ostream&, const SGCellConstraints& );

  //Immutable symmetry information of one of the 530 space group settings,
  //derived from its Hall symbol. Instances are created on first use and live
  //until the end of the program, so references to them remain valid. Obtain
  //them via SGSymmetry::get(spacegroup).
  //
  //The operations are listed in a fixed canonical order: first the coset
  //representatives (one operation per distinct rotation part, with the
  //smallest translation), in sorted order but with the identity first,
  //followed by the same representatives combined with each additional
  //centring vector (in sorted order).

  class SGSymmetry final : private NoCopyMove {
  public:
    //Access (NB: SpaceGroup is taken by value, which unlike member functions
    //or const reference parameters does not trigger false positives from
    //gcc's -Wdangling-reference when called with a temporary SpaceGroup):
    static const SGSymmetry& get( SpaceGroup );

    SpaceGroup spaceGroup() const noexcept { return m_sg; }

    //All operations (incl. centring), coset representatives, and centring
    //vectors (as pure translations, starting with the zero translation):
    const std::vector<SymOp>& operations() const noexcept { return m_ops; }
    const std::vector<SymOp>& representatives() const noexcept;
    const std::vector<SymOp>& centringVectors() const noexcept;

    unsigned order() const noexcept;//number of operations
    char latticeSymbol() const noexcept;//'P','A','B','C','I','R' or 'F'
    bool isCentrosymmetric() const noexcept;

    //Constraints on the unit cell parameters:
    const SGCellConstraints& cellConstraints() const noexcept;

    //Only for internal usage (use SGSymmetry::get(..) instead):
    struct internal_t {};
    SGSymmetry( internal_t, SpaceGroup );
  private:
    SpaceGroup m_sg;
    char m_lattice;
    bool m_centrosymmetric;
    SGCellConstraints m_cellconstraints;
    std::vector<SymOp> m_ops, m_reps, m_centring;
  };

  //Expansion of a site (fractional coordinates) into all its
  //symmetry-equivalent positions. Positions are compared per coordinate,
  //taking periodicity into account. Images of the site closer than 5e-4 in
  //all coordinates are considered identical (i.e. the site is on a special
  //position), in which case the site is first symmetrised (replaced by its
  //average over the operations mapping it onto itself), which for example
  //turns 0.3333 into 1/3. Images which are not identical, but closer than
  //1e-2 in all coordinates, indicate a site close to a special position but
  //not on it (e.g. due to insufficient precision like 0.333 instead of 1/3),
  //in which case a BadInput exception is thrown. These thresholds are based
  //on an analysis of rounding noise and genuinely distinct positions in CIF
  //files from the Crystallography Open Database.

  struct SGSiteOrbit {
    std::vector<Vector> positions;//wrapped into [0,1), symmetrised site first
    unsigned siteSymmetryOrder;//= group order / positions.size()
  };

  SGSiteOrbit expandSiteToOrbit( const SGSymmetry&, const Vector& site );

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline const std::vector<SymOp>& SGSymmetry::representatives() const noexcept
  {
    return m_reps;
  }
  inline const std::vector<SymOp>& SGSymmetry::centringVectors() const noexcept
  {
    return m_centring;
  }
  inline unsigned SGSymmetry::order() const noexcept
  {
    return static_cast<unsigned>( m_ops.size() );
  }
  inline char SGSymmetry::latticeSymbol() const noexcept
  {
    return m_lattice;
  }
  inline bool SGSymmetry::isCentrosymmetric() const noexcept
  {
    return m_centrosymmetric;
  }
  inline const SGCellConstraints& SGSymmetry::cellConstraints() const noexcept
  {
    return m_cellconstraints;
  }
}

#endif
