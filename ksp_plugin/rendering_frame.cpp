#include "ksp_plugin/rendering_frame.hpp"

#include "ksp_plugin/celestial.hpp"


namespace principia {
namespace ksp_plugin {

BodyCentredNonRotating::BodyCentredNonRotating(
    Celestial<Barycentre> const& body) : body_(body) {};

Trajectory<Barycentre> const BodyCentredNonRotating::ApparentTrajectory(
    Trajectory<Barycentre> const& actual_trajectory) const override {
  Trajectory<Barycentre> result;

}

}  // namespace ksp_plugin
}  // namespace principia
