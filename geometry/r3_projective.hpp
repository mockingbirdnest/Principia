
#pragma once

#include "geometry/r3_element.hpp"

namespace principia {
namespace geometry {
namespace internal_r3_projective {

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

  // If the point is at infinity, returns a value that identifies it among all
  // the points at infinity.  Note that this identification is not unique, so
  // the comparison of the values returned by this function is not an
  // equivalence on the projective plane.
  double point_at_infinity() const;

  // Returns the Euclidean (non-homegeneous) coordinates of the point.  May be
  // infinities.
  Scalar const x() const;
  Scalar const y() const;

 private:
  R3Element<Scalar, Scalar, double> coordinates_;

  template<typename S>
  friend bool operator==(R3Projective<S> const& left,
                         R3Projective<S> const& right);
  template<typename S>
  friend bool operator!=(R3Projective<S> const& left,
                         R3Projective<S> const& right);
  template<typename Scalar>
  friend std::string DebugString(R3Projective<Scalar> const& r3_projective);
};

// TODO(phl): Improve these operators using |std::fma| so that they define a
// proper equavalence.
template<typename Scalar>
bool operator==(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right);
template<typename Scalar>
bool operator!=(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right);

template<typename Scalar>
std::string DebugString(R3Projective<Scalar> const& r3_projective);

template<typename Scalar>
std::ostream& operator<<(std::ostream& os,
                         R3Projective<Scalar> const& r3_projective);

}  // namespace internal_r3_projective
}  // namespace geometry
}  // namespace principia

#include "geometry/r3_projective_body.hpp"
