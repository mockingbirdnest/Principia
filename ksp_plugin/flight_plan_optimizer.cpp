#include "ksp_plugin/flight_plan_optimizer.hpp"

#include <algorithm>

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
using namespace principia::physics::_apsides;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

constexpr Speed δ_homogeneous_argument = 1 * Milli(Metre) / Second;
constexpr Acceleration time_homogeneization_factor = 1 * Metre / Pow<2>(Second);
constexpr int max_apsides = 20;

void FlightPlanOptimizer::Metric::Initialize(
    not_null<FlightPlanOptimizer*> const optimizer,
    NavigationManœuvre manœuvre,
    int const index) {
  optimizer_ = optimizer;
  manœuvre_ = std::move(manœuvre);
  index_ = index;
}

FlightPlanOptimizer& FlightPlanOptimizer::Metric::optimizer() const {
  return *check_not_null(optimizer_);
}

NavigationManœuvre const& FlightPlanOptimizer::Metric::manœuvre() const {
  return manœuvre_.value();
}

int FlightPlanOptimizer::Metric::index() const {
  return index_.value();
}

class FlightPlanOptimizer::MetricForCelestialCentre
    : public FlightPlanOptimizer::Metric {
 public:
  explicit MetricForCelestialCentre(not_null<Celestial const*> const celestial);

  virtual double Evaluate(
      HomogeneousArgument const& homogeneous_argument) const override;
  virtual Gradient<double, HomogeneousArgument> EvaluateGradient(
      HomogeneousArgument const& homogeneous_argument) const override;
  virtual double EvaluateGateauxDerivative(
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& homogeneous_argument_direction)
      const override;

 private:
  // Has no effect because this metric doesn't mix multiple quantities.
  static constexpr Length scale_ = 1 * Metre;

  not_null<Celestial const*> const celestial_;
};

FlightPlanOptimizer::MetricForCelestialCentre::MetricForCelestialCentre(
    not_null<Celestial const*> const celestial)
    : celestial_(celestial) {}

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

FlightPlanOptimizer::FlightPlanOptimizer(
    not_null<FlightPlan*> const flight_plan,
    ProgressCallback progress_callback)
    : flight_plan_(flight_plan),
      progress_callback_(std::move(progress_callback)) {}

absl::Status FlightPlanOptimizer::Optimize(int const index,
                                           Celestial const& celestial,
                                           Speed const& Δv_tolerance) {
  // We are going to repeatedly tweak the |flight_plan_|, no point in running
  // the orbit analysers.
  flight_plan_->EnableAnalysis(/*enabled=*/false);

  // The following is a copy, and is not affected by changes to the
  // |flight_plan_|.
  NavigationManœuvre const manœuvre = flight_plan_->GetManœuvre(index);

  auto const f = [this, &cache, &celestial, index, &manœuvre](
                     HomogeneousArgument const& homogeneous_argument) {
    return EvaluateDistanceToCelestialWithReplacement(
        celestial, homogeneous_argument, manœuvre, index, *flight_plan_, cache);
  };
  auto const grad_f = [this, &cache, &celestial, index, &manœuvre](
                          HomogeneousArgument const& homogeneous_argument) {
    return Evaluate𝛁DistanceToCelestialWithReplacement(
        celestial, homogeneous_argument, manœuvre, index, *flight_plan_, cache);
  };
  auto const gateaux_derivative_f =
      [this, &cache, &celestial, index, &manœuvre](
    HomogeneousArgument const& homogeneous_argument,
    Difference<HomogeneousArgument> const&
        direction_homogeneous_argument) {
    return EvaluateGateauxDerivativeOfDistanceToCelestialWithReplacement(
        celestial,
        homogeneous_argument,
        direction_homogeneous_argument,
        manœuvre,
        index,
        *flight_plan_,
        cache);
  };

  auto const status_or_solution =
      BroydenFletcherGoldfarbShanno<Length, HomogeneousArgument>(
          Homogeneize(start_argument_),
          f,
          grad_f,
          gateaux_derivative_f,
          Δv_tolerance);
  if (status_or_solution.ok()) {
    auto const& solution = status_or_solution.value();
    auto const replace_status = ReplaceBurn(
        Dehomogeneize(solution), manœuvre, index, *flight_plan_);
    flight_plan_->EnableAnalysis(/*enabled=*/true);
    return replace_status;
  } else {
    flight_plan_->EnableAnalysis(/*enabled=*/true);
    return status_or_solution.status();
  }
}

absl::Status FlightPlanOptimizer::Optimize(Metric const& metric,
                                           int const index,
                                           Speed const& Δv_tolerance) {
  // We are going to repeatedly tweak the |flight_plan_|, no point in running
  // the orbit analysers.
  flight_plan_->EnableAnalysis(/*enabled=*/false);

  // The following is a copy, and is not affected by changes to the
  // |flight_plan_|.
  NavigationManœuvre const manœuvre = flight_plan_->GetManœuvre(index);

  EvaluationCache cache;

  auto const f = [this, &cache, &celestial, index, &manœuvre, target_distance](
                     HomogeneousArgument const& homogeneous_argument) {
    auto const actual_distance = EvaluateDistanceToCelestialWithReplacement(
        celestial, homogeneous_argument, manœuvre, index, *flight_plan_, cache);
    return Pow<2>(actual_distance - target_distance);
  };
  auto const grad_f =
     [this, &cache, &celestial, index, &manœuvre, target_distance](
         HomogeneousArgument const& homogeneous_argument) {
    auto const actual_distance = EvaluateDistanceToCelestialWithReplacement(
        celestial, homogeneous_argument, manœuvre, index, *flight_plan_, cache);
    auto const actual_gradient = Evaluate𝛁DistanceToCelestialWithReplacement(
        celestial, homogeneous_argument, manœuvre, index, *flight_plan_, cache);
    return 2 * (actual_distance - target_distance) * actual_gradient;
  };
  auto const gateaux_derivative_f =
      [this, &cache, &celestial, index, &manœuvre, target_distance](
          HomogeneousArgument const& homogeneous_argument,
          Difference<HomogeneousArgument> const&
              direction_homogeneous_argument) {
    auto const actual_distance = EvaluateDistanceToCelestialWithReplacement(
        celestial, homogeneous_argument, manœuvre, index, *flight_plan_, cache);
    auto const actual_gateaux_derivative =
        EvaluateGateauxDerivativeOfDistanceToCelestialWithReplacement(
        celestial, homogeneous_argument, direction_homogeneous_argument,
        manœuvre, index, *flight_plan_, cache);
    return 2 * (actual_distance - target_distance) * actual_gateaux_derivative;
  };

  auto const status_or_solution =
      BroydenFletcherGoldfarbShanno<Square<Length>, HomogeneousArgument>(
          Homogeneize(start_argument_),
          std::bind(Metric::Evaluate, metric, _1),
          std::bind(Metric::EvaluateGradient, metric, _1),
          std::bind(Metric::EvaluateGateauxDerivative, metric, _1, _2),
          Δv_tolerance);
  if (status_or_solution.ok()) {
    auto const& solution = status_or_solution.value();
    auto const replace_status = ReplaceBurn(
        Dehomogeneize(solution), manœuvre, index, *flight_plan_);
    flight_plan_->EnableAnalysis(/*enabled=*/true);
    return replace_status;
  } else {
    flight_plan_->EnableAnalysis(/*enabled=*/true);
    return status_or_solution.status();
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
    Instant const& begin_time) const {
  auto const& celestial_trajectory = celestial.trajectory();
  auto const& vessel_trajectory = flight_plan_->GetAllSegments();
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
  for (const auto& [time, degrees_of_freedom] : periapsides) {
    distance = std::min(distance,
                        (degrees_of_freedom.position() -
                         celestial_trajectory.EvaluatePosition(time))
                            .Norm());
  }
  if (progress_callback_ != nullptr) {
    progress_callback_(*flight_plan_);
  }
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
  if (ReplaceBurn(argument, manœuvre, index, *flight_plan_).ok()) {
    distance = EvaluateDistanceToCelestial(celestial, manœuvre.initial_time());
  } else {
    // If the updated burn cannot replace the existing one (e.g., because it
    // overlaps with the next burn) return an infinite length to move the
    // optimizer away from this place.
    distance = Infinity<Length>;
  }
  CHECK_OK(flight_plan_->Replace(manœuvre.burn(), index));
  cache_.emplace(homogeneous_argument, distance);
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
  return flight_plan.Replace(burn, index);
}

}  // namespace internal
}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
