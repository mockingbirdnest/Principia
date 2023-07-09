#include "ksp_plugin/flight_plan_optimizer.hpp"

#include "physics/apsides.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimizer {
namespace internal {

using namespace principia::physics::_apsides;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

constexpr Time absolute_δt = 1 * Second;
constexpr Speed absolute_δv = 1 * Metre / Second;
constexpr int max_apsides = 20;

FlightPlanOptimizer::FlightPlanOptimizer(
    not_null<FlightPlan*> const flight_plan)
    : flight_plan_(flight_plan) {}

void FlightPlanOptimizer::Optimize(int const index,
                                   Celestial const& celestial) {
  NavigationManœuvre const manœuvre = flight_plan_->GetManœuvre(index);
  Argument const start_argument{
      .initial_time = manœuvre.burn().timing.initial_time.value(),
      .Δv = manœuvre.burn().intensity.Δv.value()};
}

Length FlightPlanOptimizer::EvaluateDistanceToCelestial(
    Celestial const& celestial,
    Instant const& begin_time,
    FlightPlan const& flight_plan) {
  auto const& celestial_trajectory = celestial.trajectory();
  auto const& vessel_trajectory = flight_plan.GetAllSegments();
  DiscreteTrajectory<Barycentric> apoapsides;
  DiscreteTrajectory<Barycentric> periapsides;
  ComputeApsides(celestial_trajectory,
                 vessel_trajectory,
                 vessel_trajectory.find(begin_time),
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

FlightPlanOptimizer::LengthGradient
FlightPlanOptimizer::Evaluate𝛁DistanceToCelestial(Celestial const& celestial,
                                                  Argument const& argument,
                                                  int const index,
                                                  FlightPlan& flight_plan) {
  Argument argument_δt = argument;
  argument_δt.initial_time += absolute_δt;

  auto const argument_Δv = argument.Δv.coordinates();
  Argument argument_δx = argument;
  argument_δx.Δv = Velocity<Frenet<Navigation>>(
      {argument_Δv.x + absolute_δv, argument_Δv.y, argument_Δv.z});
  Argument argument_δy = argument;
  argument_δy.Δv = Velocity<Frenet<Navigation>>(
      {argument_Δv.x, argument_Δv.y + absolute_δv, argument_Δv.z});
  Argument argument_δz = argument;
  argument_δz.Δv = Velocity<Frenet<Navigation>>(
      {argument_Δv.x, argument_Δv.y, argument_Δv.z + absolute_δv});

  auto const distance_δt = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δt, index, flight_plan);
  auto const distance_δx = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δx, index, flight_plan);
  auto const distance_δy = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δy, index, flight_plan);
  auto const distance_δz = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δz, index, flight_plan);
}

Length FlightPlanOptimizer::EvaluateDistanceToCelestialWithReplacement(
    Celestial const& celestial,
    Argument const& argument,
    int const index,
    FlightPlan& flight_plan) {
  NavigationManœuvre::Burn burn = flight_plan.GetManœuvre(index).burn();
  burn.intensity = {.Δv = argument.Δv};
  burn.timing = {.initial_time = argument.initial_time};
  flight_plan.Replace(burn, index);
  return EvaluateDistanceToCelestial(
      celestial, argument.initial_time, flight_plan);
}

}  // namespace internal
}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
