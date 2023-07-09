#include "ksp_plugin/flight_plan_optimizer.hpp"

#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/reference_frame.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimizer {
namespace internal {

using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_reference_frame;

struct OptimizationVector {
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
};

FlightPlanOptimizer::FlightPlanOptimizer(
    not_null<FlightPlan*> const flight_plan)
    : flight_plan_(flight_plan) {}

void FlightPlanOptimizer::Optimize(int const index, Celestial const& goal) {
  NavigationManœuvre const manœuvre = flight_plan_->GetManœuvre(index);
  OptimizationVector const start_argument{
      .initial_time = manœuvre.burn().timing.initial_time.value(),
      .Δv = manœuvre.burn().intensity.Δv.value()};
}

}  // namespace internal
}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
