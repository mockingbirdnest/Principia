
#pragma once

#include <string>

namespace principia {
namespace geometry {
namespace internal_rp2_point {

// An |RP2Point| is an element of ℝP², the real projective plane.  |Scalar|
// must be a 1-dimensional vector space over ℝ, typically represented by
// |Quantity| or |double|.
template<typename Scalar>
class RP2Point {
 public:
  // Constructs a point from a set of homogeneous coordinates.  This class
  // handles zeros properly and returns infinities of the correct sign.
  explicit RP2Point(Scalar const& x, Scalar const& y, double z);

  // Returns true if and only if the point is at infinity.
  bool is_at_infinity() const;

  // Returns the Euclidean (inhomogeneous) coordinates of the point.  May be
  // (simultaneously) infinities.  The sign of infinities and zeroes is in the
  // proper quadrant.
  Scalar const x() const;
  Scalar const y() const;

 private:
  Scalar x_;
  Scalar y_;
  double z_;

  template<typename S>
  friend bool operator==(RP2Point<S> const& left,
                         RP2Point<S> const& right);
  template<typename S>
  friend bool operator!=(RP2Point<S> const& left,
                         RP2Point<S> const& right);
  template<typename Scalar>
  friend std::string DebugString(RP2Point<Scalar> const& rp2_point);
};

// These operators are implemented using exact multiplication and define a
// proper equivalence.
template<typename Scalar>
bool operator==(RP2Point<Scalar> const& left,
                RP2Point<Scalar> const& right);
template<typename Scalar>
bool operator!=(RP2Point<Scalar> const& left,
                RP2Point<Scalar> const& right);

template<typename Scalar>
std::string DebugString(RP2Point<Scalar> const& rp2_point);

template<typename Scalar>
std::ostream& operator<<(std::ostream& os,
                         RP2Point<Scalar> const& rp2_point);

}  // namespace internal_rp2_point

using internal_rp2_point::RP2Point;

}  // namespace geometry
}  // namespace principia

#include "geometry/rp2_point_body.hpp"
