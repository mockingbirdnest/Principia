
#pragma once

#include "geometry/perspective.hpp"

#include <algorithm>
#include <deque>
#include <limits>
#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

using geometry::InnerProduct;
using numerics::SolveQuadraticEquation;
using quantities::Product;
using quantities::Square;

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
Scalar const& Perspective<FromFrame, ToFrame, Scalar, LinearMap>::focal()
    const {
  return focal_;
}

template<typename FromFrame,
         typename ToFrame,
         typename Scalar,
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
std::experimental::optional<Segment<Vector<Scalar, FromFrame>>>
Perspective<FromFrame, ToFrame, Scalar, LinearMap>::SegmentBehindFocalPlane(
    Segment<Vector<Scalar, FromFrame>> const& segment) const {
  Vector<double, FromFrame> const z =
      from_camera_.linear_map()(Vector<double, ToFrame>({0.0, 0.0, 1.0}));
  Scalar const z1 = InnerProduct(segment.first - camera_, z);
  Scalar const z2 = InnerProduct(segment.second - camera_, z);
  bool const first_is_visible = z1 >= focal_;
  bool const second_is_visible = z2 >= focal_;
  if (first_is_visible && second_is_visible) {
    return segment;
  } else if (!first_is_visible && !second_is_visible) {
    return std::experimental::nullopt;
  } else {
    // λ determines where the segment intersects the focal plane.
    double const λ = (focal_ - z2) / (z1 - z2);
    auto const intercept =
        Barycentre<Point<Vector<Scalar, FromFrame>>, double>(segment,
                                                             {λ, 1.0 - λ});
    if (first_is_visible) {
      return
          std::experimental::make_optional<Segment<Vector<Scalar, FromFrame>>>(
              {segment.first, intercept});
    } else {
      CHECK(second_is_visible);
      return 
          std::experimental::make_optional<Segment<Vector<Scalar, FromFrame>>>(
              {intercept, segment.second});
    }
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
double Perspective<FromFrame, ToFrame, Scalar, LinearMap>::SphereSin²HalfAngle(
    Sphere<Scalar, FromFrame> const& sphere) const {
  // See VisibleSegments for the notation.
  Point<Vector<Scalar, FromFrame>> const& K = camera_;
  Point<Vector<Scalar, FromFrame>> const& C = sphere.centre();
  Vector<Scalar, FromFrame> const KC = C - K;
  auto const KC² = InnerProduct(KC, KC);
  // Return 1.0 if the camera is within the sphere.
  return std::min(1.0, sphere.radius²() / KC²);
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

  Vector<Scalar, FromFrame> const KA = A - K;
  Vector<Scalar, FromFrame> const KB = B - K;
  Vector<Scalar, FromFrame> const KC = C - K;

  // Bail out if the camera is inside the sphere: the segment is completely
  // hidden.  This is not a common situation, but it would cause no end of
  // trouble below, so we might as well handle it first.
  auto const KC² = InnerProduct(KC, KC);
  if (KC² <= sphere.radius²()) {
    return {};
  }

  // Consider the plane that contains K and is orthogonal to KC.  If the segment
  // AB is entirely in the half-space that doesn't contain C, there is no
  // intersection and no hiding.
  auto const KAKC = InnerProduct(KA, KC);
  auto const KBKC = InnerProduct(KB, KC);
  if (KAKC <= Square<Scalar>{} && KBKC <= Square<Scalar>{}) {
    return {segment};
  }

  // H is the projection of C on the plane KAB.  It is such that:
  //   KH = ɑ * KA + β * KB
  // where ɑ and β are computed by solving the linear system:
  //   KA·KH = KA·KC = ɑ * KA² + β * KA·KB
  //   KB·KH = KB·KC = ɑ * KA·KB + β * KB²
  auto const KA² = InnerProduct(KA, KA);
  auto const KB² = InnerProduct(KB, KB);
  auto const KAKB = InnerProduct(KA, KB);

  auto const determinant = KA² * KB² - KAKB * KAKB;
  double const ɑ = (KB² * KAKC - KAKB * KBKC) / determinant;
  double const β = (KA² * KBKC - KAKB * KAKC) / determinant;
  Vector<Scalar, FromFrame> const KH = ɑ * KA + β * KB;

  // The basic check: if H is outside or on the sphere, the sphere doesn't
  // intersect the plane KAB and therefore there is no hiding.
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

  // S is the intersection of the line AB with the line orthogonal to KC that
  // contains K.  It is such that:
  //   KS = KA + σ * AB
  // Note that σ may be infinite if AB ⊥ KH.  In this case AB.KH is +0.0 and the
  // infinity has the sign of the numerator.  The same applies to τ below.
  Vector<Scalar, FromFrame> const AB = B - A;
  auto const ABKH = InnerProduct(AB, KH);
  double const σ = -KAKH / ABKH;

  // Q is the intersection of the cone with the line AB.  It is such that:
  //   KQ = KA + λ * AB
  // T is the intersection of the line AB with the line orthogonal to KC that
  // contains P.  It is such that:
  //   KT = KA + τ * AB
  // Note that the two choices of P yield the same point T and the same value of
  // τ.  The only part of the cone that is relevant for intersections is the
  // part that is on the other side of T with respect to S.
  // For each solution δ of the above quadratic equation, this loop computes PH
  // and obtains the values of λ and τ.
  std::vector<double> λs;
  λs.reserve(4);
  bool λ_lt_σ = false;
  bool λ_gt_σ = false;
  double infinity = 0.0;  // Might be NaN in cases where it's not used.
  auto const KH² = InnerProduct(KH, KH);
  for (double const δ : δs) {
    double const γ = (r² - δ * KBKH) / KAKH;
    Vector<Scalar, FromFrame> const PH = γ * KA + δ * KB;
    auto const KAPH = InnerProduct(KA, PH);
    auto const KHPH = InnerProduct(KH, PH);
    auto const ABPH = InnerProduct(AB, PH);
    double const τ = (KH² - KHPH - KAKH) / ABKH;
    infinity = std::numeric_limits<double>::infinity() * (τ - σ);
    // Keep track of the location of Q with respect to S.
    double const λ = -KAPH / ABPH;
    λ_lt_σ |= λ < σ;
    λ_gt_σ |= λ > σ;
    // Only keep this λ if it intersects the cone in the "interesting" part.
    // We need to distinguish the case where A, B are in the same order as S, T
    // in the line AB (σ <= τ) from the case where they are in opposite orders
    // (τ <= σ).
    // If AB ⊥ KH, both σ and τ are infinite.  They have the same sign if AB is
    // behind T and opposite signs if AB is between S and T.  Hence the
    // equalities in the comparisons below.
    if ((σ <= τ && τ < λ) || (τ <= σ && λ < τ)) {
      λs.push_back(λ);
    }
  }

  // If we have intersections on both sides of S, this is the hyperbolic case.
  // Add a point at infinity.  This is correct even if λs is empty, i.e. if the
  // intersection with the cone+sphere system is in the sphere.
  if (λ_lt_σ && λ_gt_σ) {
    CHECK_GE(1, λs.size());
    λs.push_back(infinity);
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

  // The line AB doesn't have an interesting intersection with the cone and
  // doesn't intersect the sphere.  There is no hiding.
  if (λs.empty()) {
    return {segment};
  }

  // Now we have all the interesting intersections of the cone+sphere with the
  // line AB.  Determine which ones fall in the segment AB and compute the final
  // result.
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
