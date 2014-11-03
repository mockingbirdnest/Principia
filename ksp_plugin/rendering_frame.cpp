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

BodyCentredNonRotatingFrame::BodyCentredNonRotatingFrame(
    Celestial<Barycentre> const& body) : body_(body) {}

std::unique_ptr<Trajectory<Barycentre>>
BodyCentredNonRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentre> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentre>> result =
      std::make_unique<Trajectory<Barycentre>>(actual_trajectory.body());
  Trajectory<Barycentre>::NativeIterator actual_timeline =
      actual_trajectory.first();
  Trajectory<Barycentre>::NativeIterator reference_body_timeline =
      body_.history().first();
  Trajectory<Barycentre>::Timeline::const_iterator it_in_reference =
      reference_body_timeline.lower_bound(actual_timeline.begin()->first);
  DegreesOfFreedom<Barycentre> const& current_reference_state =
      body_.prolongation().last().degrees_of_freedom();
  CHECK(it_in_reference != reference_body_timeline.end());
  for (auto const& pair : actual_timeline) {
    Instant const& t = pair.first;
    DegreesOfFreedom<Barycentre> const& actual_state = pair.second;
    while (it_in_reference->first < t) {
      ++it_in_reference;
    }
    if (it_in_reference->first == t) {
      DegreesOfFreedom<Barycentre> const& reference_state =
          it_in_reference->second;
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
  Trajectory<Barycentre>::NativeIterator const& actual_timeline =
      actual_trajectory.first();
  Trajectory<Barycentre>::NativeIterator const& primary_timeline =
      primary_.history().first();
  Trajectory<Barycentre>::NativeIterator const& secondary_timeline =
      secondary_.history().first();
  Trajectory<Barycentre>::Timeline::const_iterator it_in_primary =
      primary_timeline.lower_bound(actual_timeline.begin()->first);
  Trajectory<Barycentre>::Timeline::const_iterator it_in_secondary =
      secondary_timeline.lower_bound(actual_timeline.begin()->first);
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
  for (auto const& pair : actual_timeline) {
    Instant const& t = pair.first;
    DegreesOfFreedom<Barycentre> const& actual_state = pair.second;
    while (it_in_primary->first < t) {
      ++it_in_primary;
      ++it_in_secondary;
    }
    if (it_in_primary->first == t) {
      DegreesOfFreedom<Barycentre> const& primary_state = it_in_primary->second;
      DegreesOfFreedom<Barycentre> const& secondary_state =
          it_in_secondary->second;
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
