#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "base/array.hpp"
#include "geometry/rp2_point.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "geometry/sphere.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace _perspective {
namespace internal {

using namespace principia::base::_array;
using namespace principia::geometry::_rp2_point;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::geometry::_sphere;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Frame>
using Segment = std::pair<Position<Frame>, Position<Frame>>;

template<typename Frame>
using Segments = std::vector<Segment<Frame>>;

// A perspective using the pinhole camera model.  It project a point of
// `FromFrame` to an element of ℝP².  `ToFrame` is the frame of the camera.  In
// that frame the camera is located at the origin and looking at the positive
// z-axis.  The x- and y- axis of the camera correspond to those of ℝP².
template<typename FromFrame, typename ToFrame>
class Perspective final {
 public:
  Perspective(Similarity<ToFrame, FromFrame> const& from_camera,
              Length const& focal);
  Perspective(Similarity<FromFrame, ToFrame> const& to_camera,
              Length const& focal);

  Length const& focal() const;

  // Returns the ℝP² element resulting from the projection of `point`.  This
  // is properly defined for all points other than the camera origin.
  RP2Point<Length, ToFrame> operator()(Position<FromFrame> const& point) const;

  // Returns the part of `segment` that is behind the focal plane as seen from
  // the camera.  Returns nullopt if `segment` is entirely in front of the focal
  // plane.
  std::optional<Segment<FromFrame>>
  SegmentBehindFocalPlane(Segment<FromFrame> const& segment) const;

  // The square of the tangent of the angular distance between the given points
  // as seen from the camera.
  double Tan²AngularDistance(Position<FromFrame> const& p1,
                             Position<FromFrame> const& p2) const;

  Square<Length> SquaredDistanceFromCamera(Position<FromFrame> const& p) const;

  // Returns true iff the `point` is hidden by the `sphere` in this perspective.
  bool IsHiddenBySphere(Position<FromFrame> const& point,
                        Sphere<FromFrame> const& sphere) const;

  // Returns sin² α where α is the half angle under which the `sphere` is seen.
  double SphereSin²HalfAngle(Sphere<FromFrame> const& sphere) const;

  // Returns the (sub)segments of `segment` that are visible in this perspective
  // after taking into account the hiding by `sphere`.  The returned vector has
  // 0, 1, or 2 elements.
  BoundedArray<Segment<FromFrame>, 2> VisibleSegments(
      Segment<FromFrame> const& segment,
      Sphere<FromFrame> const& sphere) const;

  // Same as above, but for hiding with multiple spheres.  For N spheres it may
  // return N + 1 segments.
  Segments<FromFrame> VisibleSegments(
      Segment<FromFrame> const& segment,
      std::vector<Sphere<FromFrame>> const& spheres) const;

 private:
  Similarity<ToFrame, FromFrame> const from_camera_;
  Similarity<FromFrame, ToFrame> const to_camera_;
  Position<FromFrame> const camera_;
  Length const focal_;

  template<typename From, typename To>
  friend std::ostream& operator<<(std::ostream& out,
                                  Perspective<From, To> const& perspective);
};

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    Perspective<FromFrame, ToFrame> const& perspective);

}  // namespace internal

using internal::Perspective;
using internal::Segment;
using internal::Segments;

}  // namespace _perspective
}  // namespace geometry
}  // namespace principia

#include "geometry/perspective_body.hpp"
