
#pragma once

#include <vector>

#include "geometry/perspective.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

using geometry::InnerProduct;
using numerics::SolveQuadraticEquation;
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

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
std::experimental::optional<RP2Point<Scalar, ToFrame>>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::operator()(
    Point<Vector<Scalar, FromFrame>> const& point) const {
  Point<Vector<Scalar, ToFrame>> const point_in_camera = to_camera_(point);
  Vector<Scalar, ToFrame> const displacement_in_camera =
      point_in_camera - ToFrame::origin;
  R3Element<Scalar> const coordinates_in_camera =
      displacement_in_camera.coordinates();
  if (coordinates_in_camera.z < Scalar()) {
    // Do not project points that are behind the camera.
    return std::experimental::nullopt;
  } else {
    // This is the actual pinhole camera projection.
    return RP2Point<Scalar, ToFrame>(coordinates_in_camera.x,
                                     coordinates_in_camera.y,
                                     coordinates_in_camera.z / focal_);
  }
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

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
std::vector<Segment<Vector<Scalar, FromFrame>>>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::VisibleSegments(
    Segment<Vector<Scalar, FromFrame>> const& segment,
    Sphere<Scalar, FromFrame> const& sphere) const {
  // K is the position of the camera, A and B the extremities of the segment,
  // C the centre of the sphere.
  Point<Vector<Scalar, FromFrame>> const& K = camera_;
  Point<Vector<Scalar, FromFrame>> const& A = segment.first;
  Point<Vector<Scalar, FromFrame>> const& B = segment.second;
  Point<Vector<Scalar, FromFrame>> const& C = sphere.centre();

  // H is the projection of C on the plane KAB.  It is such that:
  //   KH = ɑ * KA + β * KB
  // where ɑ and β are computed by solving a linear system.
  Vector<Scalar, FromFrame> const KA = A - K;
  Vector<Scalar, FromFrame> const KB = B - K;
  Vector<Scalar, FromFrame> const KC = C - K;
  auto const KA² = InnerProduct(KA, KA);
  auto const KB² = InnerProduct(KB, KB);
  auto const KAKB = InnerProduct(KA, KB);
  auto const KAKC = InnerProduct(KA, KC);
  auto const KBKC = InnerProduct(KB, KC);

  auto const determinant = KA² * KB² - KAKB * KAKB;
  double const ɑ = (KB² * KAKC - KAKB * KBKC) / determinant;
  double const β = (KA² * KBKC - KAKB * KAKC) / determinant;
  Vector<Scalar, FromFrame> const KH = ɑ * KA + β * KB;

  // The basic check: if H is outside the sphere, there is no intersection and
  // thus no hiding.
  Vector<Scalar, FromFrame> const CH = KH - KC;
  auto const CH² = InnerProduct(CH, CH);
  if (CH² >= sphere.radius²()) {
    return {segment};
  }

  // If ɑ + β < 1, H is between K and the line AB.
  // If ɑ + β > 1, H is behind the line AB when seen from K.
  // If ɑ + β == 1, H is on the line AB.
  if (ɑ + β >= 1) {
    // I is the intersection of the sphere with the line AB.  It is such that:
    //   KI = KA + λ * AB
    // where λ is computed by solving a quadratic equation.
    Vector<Scalar, FromFrame> const AB = B - A;
    Vector<Scalar, FromFrame> const AC = C - A;
    auto const AB² = InnerProduct(AB, AB);
    auto const AC² = InnerProduct(AC, AC);
    auto const ABAC = InnerProduct(AB, AC);

    std::set<double> λs = SolveQuadraticEquation(/*origin=*/0.0,
                                                 /*a0=*/AC² - sphere.radius²(),
                                                 /*a1=*/-2.0 * ABAC,
                                                 /*a2=*/AB²);
    if (λs.size() != 2) {
      // The sphere doesn't intersect the line AB or is tangent to it.  There is
      // no hiding.
      return {segment};
    }
    double const λ1 = *λs.begin();
    double const λ2 = *λs.rbegin();
    if (λ1 <= 0.0 && λ2 >= 1.0) {
      // The sphere swallows the segment.
      return {};
    }
    if (λ1 > 0.0 && λ2 < 1.0) {
      // The sphere hides the middle of the segment.
      return {{A, A + λ1 * AB }, {A + λ2 * AB, B }};
    }
    if (λ1 <= 0.0) {
      // The sphere hides the beginning of the segment.
      return {{A + λ2 * AB, B }};
    }
    {
      CHECK_GE(λ2, 1.0);
      // The sphere hides the end of the segment.
      return {{A, A + λ1 * AB }};
    }
  }

  return {};
}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
