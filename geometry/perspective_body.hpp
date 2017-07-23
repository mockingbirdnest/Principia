
#pragma once

#include "geometry/perspective.hpp"

#include <algorithm>
#include <deque>
#include <vector>

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
      focal_(focal) {
  LOG(INFO)<<"nc:"<<camera_;
}

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
      focal_(focal) {
  LOG(INFO)<<"nc:"<<camera_;
}

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
Segments<Vector<Scalar, FromFrame>>
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
  // where ɑ and β are computed by solving the linear system:
  //   KA·KH = KA·KC = ɑ * KA² + β * KA·KB
  //   KB·KH = KB·KC = ɑ * KA·KB + β * KB²
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

  // The basic check: if H is outside or on the sphere, there is no intersection
  // and thus no hiding.
  Vector<Scalar, FromFrame> const CH = KH - KC;
  auto const CH² = InnerProduct(CH, CH);
  if (CH² >= sphere.radius²()) {
    return {segment};
  }

  // The sphere intersects the plane KAB.  r² = R² - CH² is the square of the
  // radius of the circle formed by this intersection.
  auto const r² = sphere.radius²() - CH²;

  // If ɑ is negative, H is on the other side of the line KB with respect to A.
  // If it is large enough in absolute value the circle doesn't intersect the
  // wedge formed by KA and KB and there is no hiding.  Same for β with respect
  // to the line KA.
  if ((ɑ < 0 && ɑ * ɑ >= r² * KB² / determinant) ||
      (β < 0 && β * β >= r² * KA² / determinant)) {
    return {segment};
  }

  // P is a point of the plane KAB where a line going through K is tangent to
  // the circle .  It is such that:
  //   PH = γ * KA + δ * KB
  // where γ and δ are computed by solving the system:
  //   PH² = r²
  //   PH·KH = r²
  // We know that there are two such points because the sphere intersects the
  // plane.
  auto const KAKH = InnerProduct(KA, KH);
  auto const KBKH = InnerProduct(KB, KH);
  auto const a0 = r² * (r² * KA² - KAKH * KAKH);
  auto const a1 = 2.0 * r² * (KAKB * KAKH - KA² * KBKH);
  auto const a2 =
      KB² * KAKH * KAKH - 2.0 * KAKB * KAKH * KBKH + KA² * KBKH * KBKH;
  std::vector<double> const δs =
      SolveQuadraticEquation(/*origin=*/0.0, a0, a1, a2);
  CHECK_EQ(2, δs.size()) << "a0:" << a0 << " a1:" << a1 << " a2:" << a2
                         << "\nK:" << K << " A:" << A << " B:" << B
                         << " C:" << C;

  // The λs define points Q where the line AB intersects the cone+sphere system,
  // according to the formula:
  //   KQ = KA + λ * AB
  // There can be between 0 and 4 values of λ.
  std::vector<double> λs;
  λs.reserve(4);

  // For each solution of the above quadratic equation, compute the value of λ,
  // if any.
  Vector<Scalar, FromFrame> const AB = B - A;
  for (double const δ : δs) {
    double const γ = (r² - δ * KBKH) / KAKH;
    Vector<Scalar, FromFrame> const PH = γ * KA + δ * KB;
    // Here we have:
    //   KP = (ɑ - γ) * KA + (β - δ) * KB
    // Thus, the value of ɑ + β - γ - δ determines where P lies with respect to
    // the line AB. If it is behind as seen from K, there is no intersection.
    if (ɑ + β - γ - δ < 1) {
      auto const KAPH = InnerProduct(KA, PH);
      auto const ABPH = InnerProduct(AB, PH);
      double const λ = -KAPH / ABPH;
      auto const KQ = KA + λ * AB;
      // It is possible for Q to lie behind the camera, in which case:
      //   KQ·KC <= 0
      // and there is no intersection.
      if (InnerProduct(KQ, KC) > Product<Scalar, Scalar>{}) {
        λs.push_back(λ);
      }
    }
  }

  // Q is the intersection of the sphere with the line AB.  It is such that:
  //   KQ = KA + μ * AB
  // where μ is computed by solving a quadratic equation:
  //   CQ² = R² = (KQ - KC)² = (μ * AB - AC)²
  Vector<Scalar, FromFrame> const AC = C - A;
  auto const AB² = InnerProduct(AB, AB);
  auto const AC² = InnerProduct(AC, AC);
  auto const ABAC = InnerProduct(AB, AC);

  std::vector<double> const μs =
      SolveQuadraticEquation(/*origin=*/0.0,
                             /*a0=*/AC² - sphere.radius²(),
                             /*a1=*/-2.0 * ABAC,
                             /*a2=*/AB²);

  // Merge and sort all the intersections.
  std::copy(μs.begin(), μs.end(), std::back_inserter(λs));
  std::sort(λs.begin(), λs.end());

  // Now we have all the possible intersections of the cone+sphere with the line
  // AB.  Determine which ones fall in the segment AB and compute the final
  // result.
  if (λs.size() < 2) {
    // The cone+sphere doesn't intersect the line AB or is tangent to it.  There
    // is no hiding.
    return {segment};
  }
  double const λ_min = *λs.begin();
  double const λ_max = *λs.rbegin();
  if (λ_min >= 1.0 || λ_max <= 0.0) {
    // All the intersections are outside of the segment AB.
    return {segment};
  }
  if (λ_min <= 0.0 && λ_max >= 1.0) {
    // The cone+sphere swallows the segment.
    return {};
  }
  if (λ_min > 0.0 && λ_max < 1.0) {
    // The cone+sphere hides the middle of the segment.
    return {{A, A + λ_min * AB}, {A + λ_max * AB, B}};
  }
  if (λ_min <= 0.0) {
    // The cone+sphere hides the beginning of the segment.
    DCHECK_GT(λ_max, 0.0);
    return {{A + λ_max * AB, B}};
  }
  {
    DCHECK_GT(λ_max, 1.0);
    DCHECK_LE(λ_min, 1.0);
    // The cone+sphere hides the end of the segment.
    return {{A, A + λ_min * AB}};
  }
}

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
         template<typename, typename> class LinearMap>
Segments<Vector<Scalar, FromFrame>>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::VisibleSegments(
    Segment<Vector<Scalar, FromFrame>> const& segment,
    std::vector<Sphere<Scalar, FromFrame>> const& spheres) const {
  std::vector<Segment<Vector<Scalar, FromFrame>>> old_segments = {segment};
  std::vector<Segment<Vector<Scalar, FromFrame>>> new_segments;
  for (auto const& sphere : spheres) {
    for (auto const& old_segment : old_segments) {
      auto new_segments_for_sphere = VisibleSegments(old_segment, sphere);
      std::move(new_segments_for_sphere.begin(),
                new_segments_for_sphere.end(),
                std::back_inserter(new_segments));
    }
    old_segments = std::move(new_segments);
    new_segments.clear();
  }
  return old_segments;
}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
