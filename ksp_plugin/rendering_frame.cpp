#include "ksp_plugin/rendering_frame.hpp"

#include "geometry/rotation.hpp"
#include "ksp_plugin/celestial.hpp"
#include "physics/transforms.hpp"
#include "quantities/quantities.hpp"

using principia::quantities::Angle;
using principia::quantities::ArcTan;
using principia::quantities::Time;
using principia::geometry::Displacement;
using principia::geometry::InnerProduct;
using principia::geometry::Rotation;
using principia::geometry::Wedge;
using principia::physics::BodyCentredNonRotatingTransformingIterator;

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
  return trajectory.first();
}

}  // namespace

BodyCentredNonRotatingFrame::BodyCentredNonRotatingFrame(
    Celestial<Barycentre> const& body) : body_(body) {}

std::unique_ptr<Trajectory<Barycentre>>
BodyCentredNonRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentre> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentre>> result =
      std::make_unique<Trajectory<Barycentre>>(actual_trajectory.body());
  // TODO(phl): Should tag the two frames differently.
  auto it =
      BodyCentredNonRotatingTransformingIterator<Barycentre, Barycentre>(
          body_.prolongation(), &actual_trajectory);
  for (; !it.at_end(); ++it) {
    result->Append(it.time(), it.degrees_of_freedom());
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
  Trajectory<Barycentre>::NativeIterator actual_it = actual_trajectory.first();
  Trajectory<Barycentre>::NativeIterator primary_it =
      LowerBound(actual_it.time(), primary_.history());
  Trajectory<Barycentre>::NativeIterator secondary_it =
      LowerBound(actual_it.time(), secondary_.history());
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
  for (; !actual_it.at_end(); ++actual_it) {
    Instant const& t = actual_it.time();
    DegreesOfFreedom<Barycentre> const& actual_state =
        actual_it.degrees_of_freedom();
    while (primary_it.time() < t) {
      ++primary_it;
      ++secondary_it;
    }
    if (primary_it.time() == t) {
      DegreesOfFreedom<Barycentre> const& primary_state =
          primary_it.degrees_of_freedom();
      DegreesOfFreedom<Barycentre> const& secondary_state =
          secondary_it.degrees_of_freedom();
      Position<Barycentre> const barycentre =
          geometry::Barycentre<Displacement<Barycentre>, Mass>(
              {primary_state.position, secondary_state.position},
              {primary_.body().mass(), secondary_.body().mass()});

      Displacement<Barycentre> const from = primary_state.position - barycentre;
      auto const wedge = Wedge(from, to);
      Angle const angle = ArcTan(wedge.Norm(), InnerProduct(from, to));
      Rotation<Barycentre, Barycentre> const rotation =
          Rotation<Barycentre, Barycentre>(angle, wedge);
      VLOG(1) << "Rotation :\n"
              << "from     : " << from << "\n"
              << "to       : " << to << "\n"
              << "wedge    : " << wedge << "\n"
              << "angle    : " << angle  << "\n"
              << "rotation : " << rotation;
      // TODO(egg): We should have a vector space structure on
      // |DegreesOfFreedom<Fries>|.
      result->Append(t,
                     {rotation(actual_state.position - barycentre) +
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
