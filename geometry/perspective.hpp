
#pragma once

#include <experimental/optional>
#include <utility>
#include <vector>

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/rp2_point.hpp"
#include "geometry/sphere.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

template<typename Vector>
using Segment = std::pair<Point<Vector>, Point<Vector>>;

template<typename Vector>
using Segments = std::vector<Segment<Vector>>;

// A perspective using the pinhole camera model.  It project a point of
// |FromFrame| to an element of ℝP².  |ToFrame| is the frame of the camera.  In
// that frame the camera is located at the origin and looking at the positive
// z-axis.  The x- and y- axis of the camera correspond to those of ℝP².
template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
class Perspective final {
 public:
  Perspective(
      AffineMap<ToFrame, FromFrame, Scalar, LinearMap> const& from_camera,
      Scalar const& focal);
  Perspective(
      AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& to_camera,
      Scalar const& focal);

  Scalar const& focal() const;

  // Returns the ℝP² element resulting from the projection of |point|.  This
  // is properly defined for all points other than the camera origin.
  RP2Point<Scalar, ToFrame> operator()(
      Point<Vector<Scalar, FromFrame>> const& point) const;

  // Returns the part of |segment| that is behind the focal plane as seen from
  // the camera.  Returns nullopt if |segment| is entirely in front of the focal
  // plane.
  std::experimental::optional<Segment<Vector<Scalar, FromFrame>>>
  SegmentBehindFocalPlane(
      Segment<Vector<Scalar, FromFrame>> const& segment) const;

  // Returns true iff the |point| is hidden by the |sphere| in this perspective.
  bool IsHiddenBySphere(Point<Vector<Scalar, FromFrame>> const& point,
                        Sphere<Scalar, FromFrame> const& sphere) const;

  // Returns sin² ɑ where ɑ is the half angle under which the |sphere| is seen.
  double SphereSin²HalfAngle(Sphere<Scalar, FromFrame> const& sphere) const;

  // Returns the (sub)segments of |segment| that are visible in this perspective
  // after taking into account the hiding by |sphere|.  The returned vector has
  // 0, 1, or 2 elements.
  Segments<Vector<Scalar, FromFrame>> VisibleSegments(
      Segment<Vector<Scalar, FromFrame>> const& segment,
      Sphere<Scalar, FromFrame> const& sphere) const;

  // Same as above, but for hiding with multiple spheres.  For N spheres it may
  // return N + 1 segments.
  Segments<Vector<Scalar, FromFrame>> VisibleSegments(
      Segment<Vector<Scalar, FromFrame>> const& segment,
      std::vector<Sphere<Scalar, FromFrame>> const& spheres) const;

 private:
  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> const from_camera_;
  AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const to_camera_;
  Point<Vector<Scalar, FromFrame>> const camera_;
  Scalar const focal_;
};

}  // namespace internal_perspective

using internal_perspective::Perspective;
using internal_perspective::Segment;

}  // namespace geometry
}  // namespace principia

#include "geometry/perspective_body.hpp"
