#include "ksp_plugin/rendering_frame.hpp"

#include "ksp_plugin/celestial.hpp"
#include "quantities/quantities.hpp"

using principia::quantities::Time;

namespace principia {
namespace ksp_plugin {

BodyCentredNonRotatingFrame::BodyCentredNonRotatingFrame(
    Celestial<Barycentre> const& body) : body_(body) {};

std::unique_ptr<Trajectory<Barycentre>> const
BodyCentredNonRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentre> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentre>> result =
      std::make_unique<Trajectory<Barycentre>>(actual_trajectory.body());
  Trajectory<Barycentre>::Timeline const& actual_timeline =
      actual_trajectory.timeline();
  Trajectory<Barycentre>::Timeline const& reference_body_timeline =
      body_.history().timeline();
  Trajectory<Barycentre>::Timeline::const_iterator it_in_reference =
      reference_body_timeline.lower_bound(actual_timeline.begin()->first);
  DegreesOfFreedom<Barycentre> const& current_reference_state =
      reference_body_timeline.rbegin()->second;
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
                      actual_state.velocity - reference_state.velocity +
                          current_reference_state.velocity});
    }
  }
  return std::move(result);
}

}  // namespace ksp_plugin
}  // namespace principia
