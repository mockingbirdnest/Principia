#include "ksp_plugin/rendering_frame.hpp"

#include "geometry/rotation.hpp"
#include "ksp_plugin/celestial.hpp"
#include "quantities/quantities.hpp"

using principia::quantities::ArcTan;
using principia::quantities::Time;
using principia::geometry::Displacement;
using principia::geometry::InnerProduct;
using principia::geometry::Rotation;
using principia::geometry::Wedge;

namespace principia {
namespace ksp_plugin {

namespace {

// Returns an iterator for the first entry in |trajectory| with a time greater
// than or equal to |t|.
// TODO(phl): This is O(N), so we might want to expose a more efficient version.
// But then it's likely that we'll just rewrite this class anyway.
Trajectory<Barycentre>::NativeIterator LowerBound(
    Instant const& t, Trajectory<Barycentre> const& trajectory) {
  for (Trajectory<Barycentre>::NativeIterator it = trajectory.first();
       !it.at_end();
       ++it) {
    if (it.time() >= t) {
      return it;
    }
  }
  LOG(FATAL) << t << " not found in trajectory";
}

}  // namespace

BodyCentredNonRotatingFrame::BodyCentredNonRotatingFrame(
    Celestial<Barycentre> const& body) : body_(body) {}

std::unique_ptr<Trajectory<Barycentre>>
BodyCentredNonRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentre> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentre>> result =
      std::make_unique<Trajectory<Barycentre>>(actual_trajectory.body());
  Trajectory<Barycentre>::NativeIterator actual_timeline =
      actual_trajectory.first();
  Trajectory<Barycentre>::NativeIterator it_in_reference =
      LowerBound(actual_timeline.time(), body_.history());
  DegreesOfFreedom<Barycentre> const& current_reference_state =
      body_.prolongation().last().degrees_of_freedom();
  for (; !actual_timeline.at_end(); ++actual_timeline) {
    Instant const& t = actual_timeline.time();
    DegreesOfFreedom<Barycentre> const& actual_state =
        actual_timeline.degrees_of_freedom();
    while (it_in_reference.time() < t) {
      ++it_in_reference;
    }
    if (it_in_reference.time() == t) {
      DegreesOfFreedom<Barycentre> const& reference_state =
          it_in_reference.degrees_of_freedom();
      // TODO(egg): We should have a vector space structure on
      // |DegreesOfFreedom<Fries>|.
      result->Append(t,
                     {actual_state.position - reference_state.position +
                          current_reference_state.position,
                      actual_state.velocity - reference_state.velocity});
    }
  }
  return std::move(result);
}

BarycentricRotatingFrame::BarycentricRotatingFrame(
    Celestial<Barycentre> const& primary,
    Celestial<Barycentre> const& secondary)
    : primary_(primary),
      secondary_(secondary) {}

std::unique_ptr<Trajectory<Barycentre>>
BarycentricRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentre> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentre>> result =
      std::make_unique<Trajectory<Barycentre>>(actual_trajectory.body());
  Trajectory<Barycentre>::NativeIterator actual_timeline =
      actual_trajectory.first();
  Trajectory<Barycentre>::NativeIterator it_in_primary =
      LowerBound(actual_timeline.time(), primary_.history());
  Trajectory<Barycentre>::NativeIterator it_in_secondary =
      LowerBound(actual_timeline.time(), secondary_.history());
  DegreesOfFreedom<Barycentre> const& current_primary_state =
      primary_.prolongation().last().degrees_of_freedom();
  DegreesOfFreedom<Barycentre> const& current_secondary_state =
      secondary_.prolongation().last().degrees_of_freedom();
  Position<Barycentre> const current_barycentre =
      geometry::Barycentre<Displacement<Barycentre>, Mass>(
          {current_primary_state.position, current_secondary_state.position},
          {primary_.body().mass(), secondary_.body().mass()});
  Displacement<Barycentre> const to =
        current_primary_state.position - current_barycentre;
  for (; !actual_timeline.at_end(); ++actual_timeline) {
    Instant const& t = actual_timeline.time();
    DegreesOfFreedom<Barycentre> const& actual_state =
        actual_timeline.degrees_of_freedom();
    while (it_in_primary.time() < t) {
      ++it_in_primary;
      ++it_in_secondary;
    }
    if (it_in_primary.time() == t) {
      DegreesOfFreedom<Barycentre> const& primary_state =
          it_in_primary.degrees_of_freedom();
      DegreesOfFreedom<Barycentre> const& secondary_state =
          it_in_secondary.degrees_of_freedom();
      Position<Barycentre> const barycentre =
          geometry::Barycentre<Displacement<Barycentre>, Mass>(
              {primary_state.position, secondary_state.position},
              {primary_.body().mass(), secondary_.body().mass()});

      Displacement<Barycentre> const from = primary_state.position - barycentre;
      auto const wedge = Wedge(from, to);
      auto const inverse_product_of_norms = 1 / (from.Norm() * to.Norm());
      Rotation<Barycentre, Barycentre> const rotate =
          Rotation<Barycentre, Barycentre>(
              ArcTan(wedge.Norm() * inverse_product_of_norms,
                     InnerProduct(from, to) * inverse_product_of_norms),
              wedge);
      // TODO(egg): We should have a vector space structure on
      // |DegreesOfFreedom<Fries>|.
      result->Append(t,
                     {rotate(actual_state.position - barycentre) +
                          current_barycentre,
                      // This would not be trivial to compute, but we don't use
                      // it...
                      actual_state.velocity});
    }
  }
  return std::move(result);
}

}  // namespace ksp_plugin
}  // namespace principia
