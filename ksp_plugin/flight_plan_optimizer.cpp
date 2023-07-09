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
using namespace principia::physics::_apsides;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_reference_frame;

constexpr int max_apsides = 20;

struct OptimizationVector {
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
};

FlightPlanOptimizer::FlightPlanOptimizer(
    not_null<FlightPlan*> const flight_plan)
    : flight_plan_(flight_plan) {}

void FlightPlanOptimizer::Optimize(int const index,
                                   Celestial const& celestial) {
  NavigationManœuvre const manœuvre = flight_plan_->GetManœuvre(index);
  OptimizationVector const start_argument{
      .initial_time = manœuvre.burn().timing.initial_time.value(),
      .Δv = manœuvre.burn().intensity.Δv.value()};
}

Length FlightPlanOptimizer::DistanceToCelestial(Celestial const& celestial) {
  auto const& celestial_trajectory = celestial.trajectory();
  auto const& vessel_trajectory = flight_plan_->GetAllSegments();
  DiscreteTrajectory<Barycentric> apoapsides;
  DiscreteTrajectory<Barycentric> periapsides;
  ComputeApsides(celestial_trajectory,
                  vessel_trajectory,
                  vessel_trajectory.begin(),
                  vessel_trajectory.end(),
                  max_apsides,
                  apoapsides,
                  periapsides);
  Length distance = Infinity<Length>;
  for (const auto& [time, degrees_of_freedom] : periapsides) {
    distance = std::min(distance,
                        (degrees_of_freedom.position() -
                         celestial_trajectory.EvaluatePosition(time))
                            .Norm());
  }
  return distance;
}

}  // namespace internal
}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
