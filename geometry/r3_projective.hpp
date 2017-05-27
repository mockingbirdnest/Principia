
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
  explicit R3Projective(R3Element<Scalar, Scalar, double> const& coordinates);

  // Returns true if and only if the point is at infinity.
  bool is_at_infinity() const;

  // If the point is at infinity, returns a value that uniquely identifies it
  // among all the points at infinity.
  double point_at_infinity() const;

  // Returns the Euclidean (non-homegeneous) coordinates of the point.  May be
  // infinities simultaneously.
  Scalar const x() const;
  Scalar const y() const;

 private:
  R3Element<Scalar, Scalar, double> coordinates_;

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

#include "geometry/r3_projective_body.hpp"
