#pragma once

#include "geometry/plane.hpp"

namespace principia {
namespace geometry {
namespace _plane {
namespace internal {

template<typename Frame>
template<typename Scalar>
Plane<Frame> Plane<Frame>::OrthogonalTo(
    Bivector<Scalar, Frame> const& binormal) {
  return Plane(binormal.coordinates());
}

template<typename Frame>
template<typename Scalar>
Plane<Frame> Plane<Frame>::OrthogonalTo(Vector<Scalar, Frame> const& normal) {
  return Plane(normal.coordinates());
}

template<typename Frame>
std::array<Bivector<double, Frame>, 2> Plane<Frame>::UnitBinormals() const {
  return {Bivector<double, Frame>(unit_), -Bivector<double, Frame>(unit_)};
}

template<typename Frame>
std::array<Vector<double, Frame>, 2> Plane<Frame>::UnitNormals() const {
  return {Vector<double, Frame>(unit_), -Vector<double, Frame>(unit_)};
}

template<typename Frame>
template<typename Scalar>
Plane<Frame>::Plane(R3Element<Scalar> const& r3_element)
    : unit_(Normalize(r3_element)) {}

template<typename Scalar, typename Frame>
Vector<Scalar, Frame> Projection(Bivector<Scalar, Frame> const& bivector,
                                 Plane<Frame> const& plane) {
  auto const some_normal = plane.UnitNormals()[0];
  Trivector<Scalar, Frame> const projection_on_plane =
      Wedge(some_normal, bivector);
  return bivector - some_normal * projection_on_plane;
}

template<typename Scalar, typename Frame>
Vector<Scalar, Frame> Projection(Vector<Scalar, Frame> const& vector,
                                 Plane<Frame> const& plane) {
  auto const some_binormal = plane.UnitBinormals()[1];
  Trivector<Scalar, Frame> const projection_on_plane =
      Wedge(some_binormal, vector);
  return vector - some_binormal * projection_on_plane;
}

}  // namespace internal
}  // namespace _plane
}  // namespace geometry
}  // namespace principia
