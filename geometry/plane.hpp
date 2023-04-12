#pragma once

#include <array>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"

namespace principia {
namespace geometry {
namespace _plane {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;

// A subspace of codimension 1 which happens to be a subspace of dimension 2.
// No notion of orientation, neither "clockwise" nor "side of the the plane".
template<typename Frame>
class Plane {
 public:
  template<typename Scalar>
  static Plane OrthogonalTo(Bivector<Scalar, Frame> const& binormal);

  template<typename Scalar>
  static Plane OrthogonalTo(Vector<Scalar, Frame> const& normal);

  std::array<Bivector<double, Frame>, 2> UnitBinormals() const;
  std::array<Vector<double, Frame>, 2> UnitNormals() const;

 private:
  template<typename Scalar>
  explicit Plane(R3Element<Scalar> const& r3_element);

  R3Element<double> unit_;

  template<typename S, typename F>
  friend Vector<S, F> Projection(Bivector<S, F> const& bivector,
                                 Plane<F> const& plane);

  template<typename S, typename F>
  friend Vector<S, F> Projection(Vector<S, F> const& vector,
                                 Plane<F> const& plane);
};

template<typename Scalar, typename Frame>
Vector<Scalar, Frame> Projection(Bivector<Scalar, Frame> const& bivector,
                                 Plane<Frame> const& plane);

template<typename Scalar, typename Frame>
Vector<Scalar, Frame> Projection(Vector<Scalar, Frame> const& vector,
                                 Plane<Frame> const& plane);

}  // namespace internal

using internal::Plane;
using internal::Projection;

}  // namespace _plane
}  // namespace geometry
}  // namespace principia

#include "geometry/plane_body.hpp"
