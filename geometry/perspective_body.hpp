
#pragma once

#include "geometry/perspective.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

using quantities::Product;

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::Perspective(
    AffineMap<ToFrame, FromFrame, Scalar, LinearMap> const& from_camera,
    Scalar const& focal)
    : from_camera_(from_camera),
      to_camera_(from_camera.Inverse()),
      camera_(from_camera_(ToFrame::origin)),
      focal_(focal) {}

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::Perspective(
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& to_camera,
    Scalar const& focal)
    : from_camera_(to_camera.Inverse()),
      to_camera_(to_camera),
      camera_(from_camera_(ToFrame::origin)),
      focal_(focal) {}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
RP2Point<Scalar, ToFrame> Perspective<FromFrame, ToFrame, Scalar, LinearMap>::
operator()(Point<Vector<Scalar, FromFrame>> const& point) const {
  Point<Vector<Scalar, ToFrame>> const point_in_camera = to_camera_(point);
  Vector<Scalar, ToFrame> const displacement_in_camera =
      point_in_camera - ToFrame::origin;
  R3Element<Scalar> const coordinates_in_camera =
      displacement_in_camera.coordinates();
  // This is the actual pinhole camera projection.
  return RP2Point<Scalar, ToFrame>(coordinates_in_camera.x,
                                   coordinates_in_camera.y,
                                   coordinates_in_camera.z / focal_);
}

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
bool Perspective<FromFrame, ToFrame, Scalar, LinearMap>::IsHiddenBySphere(
    Point<Vector<Scalar, FromFrame>> const& point,
    Sphere<Scalar, FromFrame> const& sphere) const {
  using Displacement = Vector<Scalar, FromFrame>;

  Displacement const camera_to_centre = sphere.centre() - camera_;
  Displacement const camera_to_point = point - camera_;
  Displacement const centre_to_point = point - sphere.centre();

  Product<Scalar, Scalar> const& r² = sphere.radius²();
  Product<Scalar, Scalar> const camera_to_centre² =
      InnerProduct(camera_to_centre, camera_to_centre);
  Product<Scalar, Scalar> const camera_to_point² =
      InnerProduct(camera_to_point, camera_to_point);
  Product<Scalar, Scalar> const centre_to_point² =
      InnerProduct(centre_to_point, centre_to_point);

  // If the point lies in the sphere then surely it is hidden.
  bool const is_in_sphere = centre_to_point² < r²;
  if (is_in_sphere) {
    return true;
  }

  // Squared distance between the camera and the horizon, i.e., the circle
  // where the cone from the camera tangents the sphere.  Plain Πυθαγόρας.
  Product<Scalar, Scalar> const camera_to_horizon² = camera_to_centre² - r²;

  // This implicitly gives the cosine of the angle between the centre and the
  // point as seen from the camera.
  Product<Scalar, Scalar> const inner_product =
      InnerProduct(camera_to_point, camera_to_centre);

  // This effectively compares the square cosines of (1) the angle between the
  // centre and the point as seen from the camera and (2) the angle of the cone.
  // If the point does not lie in the cone then surely it is visible.
  bool const is_in_cone =
      inner_product * inner_product > camera_to_horizon² * camera_to_point²;
  if (!is_in_cone) {
    return false;
  }

  // This effectively compares (1) the distance from the camera to the plane of
  // the horizon (the plane where the cone tangents the sphere) and (2) the
  // distance from the camera to the projection of the point on the axis camera-
  // centre.
  bool const is_in_front_of_horizon = inner_product < camera_to_horizon²;
  return !is_in_front_of_horizon;
}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
