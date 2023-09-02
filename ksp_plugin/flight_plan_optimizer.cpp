#include "ksp_plugin/flight_plan_optimizer.hpp"

#include <algorithm>

#include "geometry/barycentre_calculator.hpp"
#include "physics/apsides.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimizer {
namespace internal {

using std::placeholders::_1;
using std::placeholders::_2;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::physics::_apsides;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

// Keep this consistent with |max_steps_| in |flight_planner.cs|.
constexpr std::int64_t max_steps_in_flight_plan = 1 << 20;

// Conversion factors between |Argument| and |HomogeneousArgument|.
constexpr Time time_homogeneization_factor = 1 * Second;
constexpr Speed speed_homogeneization_factor = 1 * Metre / Second;

// The displacement in |HomogeneousArgument| used to compute the derivatives.
constexpr double δ_homogeneous_argument = 1e-3;

// By how much we extend the flight plan when it is too short to find
// periapsides.
constexpr double flight_plan_extension_factor = 1.05;

constexpr int max_apsides = 20;

class FlightPlanOptimizer::MetricForCelestialCentre
    : public FlightPlanOptimizer::Metric {
 public:
  MetricForCelestialCentre(not_null<FlightPlanOptimizer*> optimizer,
                           NavigationManœuvre manœuvre,
                           int index,
                           not_null<Celestial const*> celestial);

  double Evaluate(
      HomogeneousArgument const& homogeneous_argument) const override;
  Gradient<double, HomogeneousArgument> EvaluateGradient(
      HomogeneousArgument const& homogeneous_argument) const override;
  double EvaluateGateauxDerivative(
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& homogeneous_argument_direction)
      const override;

 private:
  // Has no effect because this metric doesn't mix multiple quantities.
  static constexpr Length scale_ = 1 * Metre;

  not_null<Celestial const*> const celestial_;
};

class FlightPlanOptimizer::MetricForCelestialDistance
    : public FlightPlanOptimizer::Metric {
 public:
  MetricForCelestialDistance(not_null<FlightPlanOptimizer*> optimizer,
                             NavigationManœuvre manœuvre,
                             int index,
                             not_null<Celestial const*> const celestial,
                             Length const& target_distance);

  double Evaluate(
      HomogeneousArgument const& homogeneous_argument) const override;
  Gradient<double, HomogeneousArgument> EvaluateGradient(
      HomogeneousArgument const& homogeneous_argument) const override;
  double EvaluateGateauxDerivative(
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& homogeneous_argument_direction)
      const override;

 private:
  // Has no effect because this metric doesn't mix multiple quantities.
  static constexpr Square<Length> scale_ = 1 * Pow<2>(Metre);

  not_null<Celestial const*> const celestial_;
  Length const target_distance_;
};

FlightPlanOptimizer::MetricForCelestialCentre::MetricForCelestialCentre(
    not_null<FlightPlanOptimizer*> const optimizer,
    NavigationManœuvre manœuvre,
    int const index,
    not_null<Celestial const*> const celestial)
    : Metric(optimizer, std::move(manœuvre), index),
      celestial_(celestial) {}

double FlightPlanOptimizer::MetricForCelestialCentre::Evaluate(
    HomogeneousArgument const& homogeneous_argument) const {
  return optimizer().EvaluateDistanceToCelestialWithReplacement(
      *celestial_, homogeneous_argument, manœuvre(), index()) / scale_;
}

Gradient<double, FlightPlanOptimizer::HomogeneousArgument>
FlightPlanOptimizer::MetricForCelestialCentre::EvaluateGradient(
    HomogeneousArgument const& homogeneous_argument) const {
  return optimizer().Evaluate𝛁DistanceToCelestialWithReplacement(
      *celestial_, homogeneous_argument, manœuvre(), index()) / scale_;
}

double FlightPlanOptimizer::MetricForCelestialCentre::EvaluateGateauxDerivative(
    HomogeneousArgument const& homogeneous_argument,
    Difference<HomogeneousArgument> const& homogeneous_argument_direction)
    const {
  return optimizer()
      .EvaluateGateauxDerivativeOfDistanceToCelestialWithReplacement(
          *celestial_,
          homogeneous_argument,
          homogeneous_argument_direction,
          manœuvre(),
          index()) / scale_;
}

FlightPlanOptimizer::MetricForCelestialDistance::MetricForCelestialDistance(
    not_null<FlightPlanOptimizer*> const optimizer,
    NavigationManœuvre manœuvre,
    int const index,
    not_null<Celestial const*> const celestial,
    Length const& target_distance)
    : Metric(optimizer, manœuvre, index),
      celestial_(celestial),
      target_distance_(target_distance) {}

double FlightPlanOptimizer::MetricForCelestialDistance::Evaluate(
    HomogeneousArgument const& homogeneous_argument) const {
  Length const actual_distance =
      optimizer().EvaluateDistanceToCelestialWithReplacement(
          *celestial_, homogeneous_argument, manœuvre(), index());
  return Pow<2>(actual_distance - target_distance_) / scale_;
}

Gradient<double, FlightPlanOptimizer::HomogeneousArgument>
FlightPlanOptimizer::MetricForCelestialDistance::EvaluateGradient(
    HomogeneousArgument const& homogeneous_argument) const {
  Length const actual_distance =
      optimizer().EvaluateDistanceToCelestialWithReplacement(
          *celestial_, homogeneous_argument, manœuvre(), index());
  Gradient<Length, FlightPlanOptimizer::HomogeneousArgument> const
      actual_gradient = optimizer().Evaluate𝛁DistanceToCelestialWithReplacement(
          *celestial_, homogeneous_argument, manœuvre(), index());
  return 2 * (actual_distance - target_distance_) * actual_gradient / scale_;
}

double
FlightPlanOptimizer::MetricForCelestialDistance::EvaluateGateauxDerivative(
    HomogeneousArgument const& homogeneous_argument,
    Difference<HomogeneousArgument> const& homogeneous_argument_direction)
    const {
  Length const actual_distance =
      optimizer().EvaluateDistanceToCelestialWithReplacement(
          *celestial_, homogeneous_argument, manœuvre(), index());
  Length const actual_gateaux_derivative =
      optimizer().EvaluateGateauxDerivativeOfDistanceToCelestialWithReplacement(
          *celestial_,
          homogeneous_argument,
          homogeneous_argument_direction,
          manœuvre(),
          index());
  return 2 * (actual_distance - target_distance_) * actual_gateaux_derivative /
         scale_;
}

FlightPlanOptimizer::HomogeneousArgument FlightPlanOptimizer::Homogeneize(
    Argument const& argument) {
  auto const& ΔΔv_coordinates = argument.ΔΔv.coordinates();
  return HomogeneousArgument(
      {argument.Δinitial_time / time_homogeneization_factor,
       ΔΔv_coordinates.x / speed_homogeneization_factor,
       ΔΔv_coordinates.y / speed_homogeneization_factor,
       ΔΔv_coordinates.z / speed_homogeneization_factor});
}

FlightPlanOptimizer::Argument FlightPlanOptimizer::Dehomogeneize(
    HomogeneousArgument const& homogeneous_argument) {
  return Argument{
      .Δinitial_time = homogeneous_argument[0] * time_homogeneization_factor,
      .ΔΔv = Velocity<Frenet<Navigation>>(
          {homogeneous_argument[1] * speed_homogeneization_factor,
           homogeneous_argument[2] * speed_homogeneization_factor,
           homogeneous_argument[3] * speed_homogeneization_factor})};
}

FlightPlanOptimizer::Metric::Metric(
    not_null<FlightPlanOptimizer*> const optimizer,
    NavigationManœuvre manœuvre,
    int const index)
    : optimizer_(optimizer),
      manœuvre_(std::move(manœuvre)),
      index_(index) {}

FlightPlanOptimizer& FlightPlanOptimizer::Metric::optimizer() const {
  return *optimizer_;
}

NavigationManœuvre const& FlightPlanOptimizer::Metric::manœuvre() const {
  return manœuvre_;
}

int FlightPlanOptimizer::Metric::index() const {
  return index_;
}

FlightPlanOptimizer::MetricFactory FlightPlanOptimizer::ForCelestialCentre(
    not_null<Celestial const*> const celestial) {
  return [celestial](not_null<FlightPlanOptimizer*> const optimizer,
                     NavigationManœuvre manœuvre,
                     int const index) {
    return make_not_null_unique<MetricForCelestialCentre>(
        optimizer, manœuvre, index, celestial);
  };
}

FlightPlanOptimizer::MetricFactory FlightPlanOptimizer::ForCelestialDistance(
    not_null<Celestial const*> const celestial,
    Length const& target_distance) {
  return [celestial, target_distance](
             not_null<FlightPlanOptimizer*> const optimizer,
             NavigationManœuvre manœuvre,
             int const index) {
    return make_not_null_unique<MetricForCelestialDistance>(
        optimizer, manœuvre, index, celestial, target_distance);
  };
}

FlightPlanOptimizer::FlightPlanOptimizer(
    not_null<FlightPlan*> const flight_plan,
    MetricFactory metric_factory,
    ProgressCallback progress_callback)
    : flight_plan_(flight_plan),
      metric_factory_(std::move(metric_factory)),
      progress_callback_(std::move(progress_callback)) {}

absl::Status FlightPlanOptimizer::Optimize(int const index,
                                           Speed const& Δv_tolerance) {
  // We are going to repeatedly tweak the |flight_plan_|, no point in running
  // the orbit analysers.
  flight_plan_->EnableAnalysis(/*enabled=*/false);

  // Don't reuse the computations from the previous optimization.
  cache_.clear();

  // The following is a copy, and is not affected by changes to the
  // |flight_plan_|.  It is moved into the metric.
  NavigationManœuvre manœuvre = flight_plan_->GetManœuvre(index);
  auto const metric = metric_factory_(this, std::move(manœuvre), index);

  auto const status_or_solution =
      BroydenFletcherGoldfarbShanno<double, HomogeneousArgument>(
          Homogeneize(start_argument_),
          std::bind(&Metric::Evaluate, metric.get(), _1),
          std::bind(&Metric::EvaluateGradient, metric.get(), _1),
          std::bind(&Metric::EvaluateGateauxDerivative, metric.get(), _1, _2),
          Δv_tolerance / speed_homogeneization_factor);
  if (status_or_solution.ok()) {
    auto const& solution = status_or_solution.value();
    auto const replace_status =
        ReplaceBurn(Dehomogeneize(solution), manœuvre, index);
    flight_plan_->EnableAnalysis(/*enabled=*/true);
    return replace_status;
  } else {
    flight_plan_->EnableAnalysis(/*enabled=*/true);
    return status_or_solution.status();
  }
}

Length FlightPlanOptimizer::EvaluateDistanceToCelestial(
    Celestial const& celestial,
    Instant const& begin_time,
    bool const extend_if_needed) const {
  auto const& celestial_trajectory = celestial.trajectory();
  auto const& vessel_trajectory = flight_plan_->GetAllSegments();

  Length distance_at_closest_periapsis;
  for (;;) {
    DiscreteTrajectory<Barycentric> apoapsides;
    DiscreteTrajectory<Barycentric> periapsides;
    ComputeApsides(celestial_trajectory,
                   vessel_trajectory,
                   vessel_trajectory.lower_bound(begin_time),
                   vessel_trajectory.end(),
                   max_apsides,
                   apoapsides,
                   periapsides);
    distance_at_closest_periapsis = Infinity<Length>;
    for (const auto& [time, degrees_of_freedom] : periapsides) {
      distance_at_closest_periapsis =
          std::min(distance_at_closest_periapsis,
                   (degrees_of_freedom.position() -
                    celestial_trajectory.EvaluatePosition(time))
                       .Norm());
    }

    // Evaluate the distance at the end of the trajectory.  If it is smaller
    // than all the periapsides, increase the length of the flight plan until it
    // isn't.
    auto const distance_at_end =
        (vessel_trajectory.rbegin()->degrees_of_freedom.position() -
         celestial_trajectory.EvaluatePosition(celestial_trajectory.t_max()))
            .Norm();
    if (distance_at_end >= distance_at_closest_periapsis) {
      break;
    } else if (!extend_if_needed) {
      return distance_at_end;
    }

    // Try to nudge the desired final time.  This may not succeed, in which case
    // we give up.
    auto const previous_actual_final_time = flight_plan_->actual_final_time();
    auto const new_desired_final_time = Barycentre<Instant, double>(
        {flight_plan_->initial_time(), flight_plan_->desired_final_time()},
        {1 - flight_plan_extension_factor, flight_plan_extension_factor});
    flight_plan_->SetDesiredFinalTime(new_desired_final_time).IgnoreError();
    if (flight_plan_->actual_final_time() <= previous_actual_final_time) {
      break;
    }
  }

  return distance_at_closest_periapsis;
}

Length FlightPlanOptimizer::EvaluateDistanceToCelestialWithReplacement(
    Celestial const& celestial,
    HomogeneousArgument const& homogeneous_argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  if (auto const it = cache_.find(homogeneous_argument); it != cache_.end()) {
    return it->second;
  }

  Length distance;
  Argument const argument = Dehomogeneize(homogeneous_argument);
  auto const replace_status = ReplaceBurn(argument, manœuvre, index);
  if (progress_callback_ != nullptr) {
    progress_callback_(*flight_plan_);
  }
  if (replace_status.ok()) {
    distance = EvaluateDistanceToCelestial(celestial,
                                           manœuvre.initial_time(),
                                           /*extend_if_needed=*/true);
  } else {
    // If the burn could not be replaced, e.g., because the integrator reached
    // its maximal number of steps, evaluate the distance as best as we can,
    // without trying to be smart and extend the flight plan.  This is somewhat
    // iffy, but better than the alternative of returning an infinity, which
    // introduces discontinuities.
    distance = EvaluateDistanceToCelestial(celestial,
                                           manœuvre.initial_time(),
                                           /*extend_if_needed=*/false);
  }
  flight_plan_->Replace(manœuvre.burn(), index).IgnoreError();
  cache_.emplace(homogeneous_argument, distance);
  return distance;
}

FlightPlanOptimizer::LengthGradient
FlightPlanOptimizer::Evaluate𝛁DistanceToCelestialWithReplacement(
    Celestial const& celestial,
    HomogeneousArgument const& homogeneous_argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  auto const distance = EvaluateDistanceToCelestialWithReplacement(
      celestial, homogeneous_argument, manœuvre, index);

  LengthGradient gradient;
  for (int i = 0; i < HomogeneousArgument::dimension; ++i) {
    HomogeneousArgument homogeneous_argument_δi = homogeneous_argument;
    homogeneous_argument_δi[i] += δ_homogeneous_argument;
    auto const distance_δi =
        EvaluateDistanceToCelestialWithReplacement(celestial,
                                                   homogeneous_argument_δi,
                                                   manœuvre,
                                                   index);
    gradient[i] = (distance_δi - distance) / δ_homogeneous_argument;
  }
  return gradient;
}

Length FlightPlanOptimizer::
EvaluateGateauxDerivativeOfDistanceToCelestialWithReplacement(
    Celestial const& celestial,
    HomogeneousArgument const& homogeneous_argument,
    Difference<HomogeneousArgument> const& direction_homogeneous_argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  auto const distance = EvaluateDistanceToCelestialWithReplacement(
      celestial, homogeneous_argument, manœuvre, index);
  double const h = δ_homogeneous_argument /
                   direction_homogeneous_argument.Norm();
  auto const homogeneous_argument_h =
      homogeneous_argument + h * direction_homogeneous_argument;
  auto const distance_δh = EvaluateDistanceToCelestialWithReplacement(
      celestial, homogeneous_argument_h, manœuvre, index);
  return (distance_δh - distance) / h;
}

absl::Status FlightPlanOptimizer::ReplaceBurn(
    Argument const& argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  NavigationManœuvre::Burn burn = manœuvre.burn();
  burn.intensity = {.Δv = manœuvre.Δv() + argument.ΔΔv};
  burn.timing = {.initial_time =
                     manœuvre.initial_time() + argument.Δinitial_time};
  for (;;) {
    auto status = flight_plan_->Replace(burn, index);
    if (true ||
        status.code() != integrators::_ordinary_differential_equations::
                             termination_condition::ReachedMaximalStepCount) {
      return status;
    }

    // If the flight plan integrator reached the maximum step count, increase it
    // if possible.  Keep doing this until we get a different error (or we
    // succeed).
    do {
      auto adaptive_step_parameters = flight_plan_->adaptive_step_parameters();
      auto generalized_adaptive_step_parameters =
          flight_plan_->generalized_adaptive_step_parameters();
      if (adaptive_step_parameters.max_steps() == max_steps_in_flight_plan ||
          generalized_adaptive_step_parameters.max_steps() ==
              max_steps_in_flight_plan) {
        return status;
      }
      adaptive_step_parameters.set_max_steps(
          adaptive_step_parameters.max_steps() * 4);
      generalized_adaptive_step_parameters.set_max_steps(
          generalized_adaptive_step_parameters.max_steps() * 4);
      status = flight_plan_->SetAdaptiveStepParameters(
          adaptive_step_parameters, generalized_adaptive_step_parameters);
    } while (status.code() ==
             integrators::_ordinary_differential_equations::
                 termination_condition::ReachedMaximalStepCount);
  }
}

}  // namespace internal
}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
