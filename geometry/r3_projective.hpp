
#pragma once

#include "geometry/r3_element.hpp"

namespace principia {
namespace geometry {

// An |R3Projective| is an element of ℝP², the real projective plane.  |Scalar|
// must be a 1-dimensional vector space over ℝ, typically represented by
// |Quantity| or |double|.
template<typename Scalar>
class R3Projective {
 public:
  // Constructs a point from a set of homogeneous coordinates.  The
  // |coordinates| may not all be zero.
  explicit R3Projective(R3Element<Scalar> const& coordinates);

  // Returns true if and only if the point is at infinity.
  bool at_infinity() const;

  // Returns the slope of the line going from the origin to this point.  Useful
  // to distinguish the points at infinity.  Note that the returned value may be
  // an infinity.
  Scalar slope() const;

  // Returns the Euclidean (non-homegeneous) coordinates of the point.  May be
  // infinities simultaneously.
  Scalar const x() const;
  Scalar const y() const;

 private:
  R3Element<Scalar> coordinates_;

  friend bool operator==(R3Projective<Scalar> const& left,
                         R3Projective<Scalar> const& right);
  friend bool operator!=(R3Projective<Scalar> const& left,
                         R3Projective<Scalar> const& right);
};

template<typename Scalar>
bool operator==(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right);
template<typename Scalar>
bool operator!=(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/projective_body.hpp"
