#include "ksp_plugin/flight_plan_optimizer.hpp"

#include "physics/apsides.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimizer {
namespace internal {

using namespace principia::physics::_apsides;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

constexpr Time absolute_δt = 1 * Milli(Second);
constexpr Speed absolute_δv = 1 * Milli(Metre) / Second;
constexpr Acceleration time_homogeneization_factor = 1 * Metre / Pow<2>(Second);
constexpr int max_apsides = 20;

FlightPlanOptimizer::FlightPlanOptimizer(
    not_null<FlightPlan*> const flight_plan)
    : flight_plan_(flight_plan) {}

absl::Status FlightPlanOptimizer::Optimize(int const index,
                                           Celestial const& celestial,
                                           Speed const& Δv_tolerance) {
  // The following is a copy, and is not affected by changes to the
  // |flight_plan_|.
  NavigationManœuvre const manœuvre = flight_plan_->GetManœuvre(index);
  LOG(ERROR)<<flight_plan_->GetManœuvre(index).initial_time()<<" "
<<flight_plan_->GetManœuvre(index).Δv();

  //LOG(ERROR)<<EvaluateDistanceToCelestial(celestial, manœuvre.initial_time(), *flight_plan_);

  auto const f = [this, &celestial, index, &manœuvre](
                     HomogeneousArgument const& homogeneous_argument) {
    LOG(ERROR)<<"Function";
    return EvaluateDistanceToCelestialWithReplacement(
        celestial,
        Dehomogeneize(homogeneous_argument),
        manœuvre,
        index,
        *flight_plan_);
  };
  auto const grad_f = [this, &celestial, index, &manœuvre](
                          HomogeneousArgument const& homogeneous_argument) {
    LOG(ERROR)<<"Gradient";
    return Evaluate𝛁DistanceToCelestialWithReplacement(
        celestial,
        Dehomogeneize(homogeneous_argument),
        manœuvre,
        index,
        *flight_plan_);
  };

  auto const solution =
      BroydenFletcherGoldfarbShanno<Length, HomogeneousArgument>(
          Homogeneize(start_argument_), f, grad_f, Δv_tolerance);
  if (solution.has_value()) {
    LOG(ERROR)<<solution.value();
    return ReplaceBurn(
        Dehomogeneize(solution.value()), manœuvre, index, *flight_plan_);
  } else {
    return absl::NotFoundError("No better burn");
  }
}

FlightPlanOptimizer::HomogeneousArgument FlightPlanOptimizer::Homogeneize(
    Argument const& argument) {
  auto const& ΔΔv_coordinates = argument.ΔΔv.coordinates();
  return HomogeneousArgument(
      {argument.Δinitial_time * time_homogeneization_factor,
       ΔΔv_coordinates.x,
       ΔΔv_coordinates.y,
       ΔΔv_coordinates.z});
}

FlightPlanOptimizer::Argument FlightPlanOptimizer::Dehomogeneize(
    HomogeneousArgument const& homogeneous_argument) {
  return Argument{
      .Δinitial_time = homogeneous_argument[0] / time_homogeneization_factor,
      .ΔΔv = Velocity<Frenet<Navigation>>({homogeneous_argument[1],
                                           homogeneous_argument[2],
                                           homogeneous_argument[3]})};
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
                 vessel_trajectory.lower_bound(begin_time),
                 vessel_trajectory.end(),
                 max_apsides,
                 apoapsides,
                 periapsides);
  Length distance = Infinity<Length>;
  Instant t;
  for (const auto& [time, degrees_of_freedom] : periapsides) {
    //LOG(ERROR)<<time<<" "<<(degrees_of_freedom.position() -
    //                     celestial_trajectory.EvaluatePosition(time))
    //                        .Norm();
    if ((degrees_of_freedom.position() -
                         celestial_trajectory.EvaluatePosition(time))
                            .Norm() < distance) {
      t = time;
    }
    distance = std::min(distance,
                        (degrees_of_freedom.position() -
                         celestial_trajectory.EvaluatePosition(time))
                            .Norm());
  }
  LOG(ERROR)<<t<<" "<<distance;
  return distance;
}

FlightPlanOptimizer::LengthGradient
FlightPlanOptimizer::Evaluate𝛁DistanceToCelestialWithReplacement(
    Celestial const& celestial,
    Argument const& argument,
    NavigationManœuvre const& manœuvre,
    int const index,
    FlightPlan& flight_plan) {
  auto const distance = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument, manœuvre, index, flight_plan);

  Argument argument_δt = argument;
  argument_δt.Δinitial_time += absolute_δt;

  auto const argument_ΔΔv = argument.ΔΔv.coordinates();
  Argument argument_δx = argument;
  argument_δx.ΔΔv = Velocity<Frenet<Navigation>>(
      {argument_ΔΔv.x + absolute_δv, argument_ΔΔv.y, argument_ΔΔv.z});
  Argument argument_δy = argument;
  argument_δy.ΔΔv = Velocity<Frenet<Navigation>>(
      {argument_ΔΔv.x, argument_ΔΔv.y + absolute_δv, argument_ΔΔv.z});
  Argument argument_δz = argument;
  argument_δz.ΔΔv = Velocity<Frenet<Navigation>>(
      {argument_ΔΔv.x, argument_ΔΔv.y, argument_ΔΔv.z + absolute_δv});

  auto const distance_δt = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δt, manœuvre, index, flight_plan);
  auto const distance_δx = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δx, manœuvre, index, flight_plan);
  auto const distance_δy = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δy, manœuvre, index, flight_plan);
  auto const distance_δz = EvaluateDistanceToCelestialWithReplacement(
      celestial, argument_δz, manœuvre, index, flight_plan);

  LOG(ERROR)<< (distance_δt - distance) / absolute_δt<<" "<<
      (distance_δx - distance) / absolute_δv<<" "<<
      (distance_δy - distance) / absolute_δv<<" "<<
      (distance_δz - distance) / absolute_δv;

  return LengthGradient({
      (distance_δt - distance) / (absolute_δt * time_homogeneization_factor),
      (distance_δx - distance) / absolute_δv,
      (distance_δy - distance) / absolute_δv,
      (distance_δz - distance) / absolute_δv});
}

Length FlightPlanOptimizer::EvaluateDistanceToCelestialWithReplacement(
    Celestial const& celestial,
    Argument const& argument,
    NavigationManœuvre const& manœuvre,
    int const index,
    FlightPlan& flight_plan) {
  Length distance;
  if (ReplaceBurn(argument, manœuvre, index, flight_plan).ok()) {
    LOG(ERROR)<<argument.Δinitial_time<<" "<<argument.ΔΔv;
    distance = EvaluateDistanceToCelestial(
        celestial, manœuvre.initial_time(), flight_plan);
  } else {
    // If the updated burn cannot replace the existing one (e.g., because it
    // overlaps with the next burn) return an infinite length to move the
    // optimizer away from this place.
    distance = Infinity<Length>;
  }
  CHECK_OK(flight_plan.Replace(manœuvre.burn(), index));
  return distance;
}

absl::Status FlightPlanOptimizer::ReplaceBurn(
    Argument const& argument,
    NavigationManœuvre const& manœuvre,
    int const index,
    FlightPlan& flight_plan) {
  NavigationManœuvre::Burn burn = manœuvre.burn();
  burn.intensity = {.Δv = manœuvre.Δv() + argument.ΔΔv};
  burn.timing = {.initial_time =
                     manœuvre.initial_time() + argument.Δinitial_time};
  auto s = flight_plan.Replace(burn, index);
//  LOG(ERROR)<<flight_plan.GetManœuvre(index).initial_time()<<" "
//<<flight_plan.GetManœuvre(index).Δv();
  return s;
}

}  // namespace internal
}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
