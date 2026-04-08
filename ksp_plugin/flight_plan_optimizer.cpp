#include "ksp_plugin/flight_plan_optimizer.hpp"

#include <functional>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "numerics/angle_reduction.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/numbers.hpp"  // 🧙 For π.
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimizer {
namespace internal {

using std::placeholders::_1;
using std::placeholders::_2;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::numerics::_angle_reduction;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;

// Conversion factors between `Argument` and `HomogeneousArgument`.
constexpr Time time_homogeneization_factor = 1 * Second;
constexpr Speed speed_homogeneization_factor = 1 * Metre / Second;

// The displacement in `HomogeneousArgument` used to compute the derivatives.
constexpr double δ_homogeneous_argument = 1e-3;

// By how much we extend the flight plan when it is too short to find
// periapsides.
constexpr double flight_plan_extension_factor = 1.05;

constexpr int max_apsides = 20;

class FlightPlanOptimizer::LinearCombinationOfMetrics
    : public FlightPlanOptimizer::Metric {
 public:
  LinearCombinationOfMetrics(not_null<FlightPlanOptimizer*> optimizer,
                             NavigationManœuvre const& manœuvre,
                             int index,
                             std::vector<MetricFactory> const& factories,
                             std::vector<double> const& weights);

  double Evaluate(
      HomogeneousArgument const& homogeneous_argument) const override;
  Gradient<double, HomogeneousArgument> EvaluateGradient(
      HomogeneousArgument const& homogeneous_argument) const override;
  double EvaluateGateauxDerivative(
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& homogeneous_argument_direction)
      const override;

 private:
  std::vector<std::unique_ptr<Metric>> metrics_;
  std::vector<double> const weights_;
};

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
                             not_null<Celestial const*> celestial,
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

class FlightPlanOptimizer::MetricForInclination
    : public FlightPlanOptimizer::Metric {
 public:
  MetricForInclination(not_null<FlightPlanOptimizer*> optimizer,
                       NavigationManœuvre manœuvre,
                       int index,
                       not_null<Celestial const*> celestial,
                       not_null<std::unique_ptr<NavigationFrame const>> frame,
                       Angle const& target_inclination);

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
  static constexpr Square<Angle> scale_ = 1 * Pow<2>(Degree);

  not_null<Celestial const*> const celestial_;
  not_null<std::unique_ptr<NavigationFrame const>> const frame_;
  Angle const target_inclination_;
};

class FlightPlanOptimizer::MetricForΔv : public FlightPlanOptimizer::Metric {
 public:
  MetricForΔv(not_null<FlightPlanOptimizer*> optimizer,
              NavigationManœuvre manœuvre,
              int index);

  double Evaluate(
      HomogeneousArgument const& homogeneous_argument) const override;
  Gradient<double, HomogeneousArgument> EvaluateGradient(
      HomogeneousArgument const& homogeneous_argument) const override;
  double EvaluateGateauxDerivative(
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& homogeneous_argument_direction)
      const override;

 private:
  NavigationManœuvre UpdatedManœuvre(
      HomogeneousArgument const& homogeneous_argument) const;

  // Has no effect because this metric doesn't mix multiple quantities.
  static constexpr Square<Speed> scale_ = 1 * Pow<2>(Metre / Second);
};

FlightPlanOptimizer::LinearCombinationOfMetrics::LinearCombinationOfMetrics(
    not_null<FlightPlanOptimizer*> optimizer,
    NavigationManœuvre const& manœuvre,
    int const index,
    std::vector<MetricFactory> const& factories,
    std::vector<double> const& weights)
    : Metric(optimizer, manœuvre, index),
      weights_(weights) {
  CHECK_EQ(factories.size(), weights.size());
  for (auto const& factory : factories) {
    metrics_.push_back(factory(optimizer, manœuvre, index));
  }
}

double FlightPlanOptimizer::LinearCombinationOfMetrics::Evaluate(
    HomogeneousArgument const& homogeneous_argument) const {
  double combined_metric = 0.0;
  for (int i = 0; i < metrics_.size(); ++i) {
    combined_metric +=
        metrics_[i]->Evaluate(homogeneous_argument) * weights_[i];
  }
  return combined_metric;
}

Gradient<double, FlightPlanOptimizer::HomogeneousArgument>
FlightPlanOptimizer::LinearCombinationOfMetrics::EvaluateGradient(
    HomogeneousArgument const& homogeneous_argument) const {
  Gradient<double, HomogeneousArgument> combined_gradient{};
  for (int i = 0; i < metrics_.size(); ++i) {
    combined_gradient +=
        metrics_[i]->EvaluateGradient(homogeneous_argument) * weights_[i];
  }
  return combined_gradient;
}

double FlightPlanOptimizer::LinearCombinationOfMetrics::
EvaluateGateauxDerivative(
    HomogeneousArgument const& homogeneous_argument,
    Difference<HomogeneousArgument> const& homogeneous_argument_direction)
    const {
  double combined_derivative = 0.0;
  for (int i = 0; i < metrics_.size(); ++i) {
    combined_derivative +=
        metrics_[i]->EvaluateGateauxDerivative(homogeneous_argument,
                                               homogeneous_argument_direction) *
        weights_[i];
  }
  return combined_derivative;
}

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
    : Metric(optimizer, std::move(manœuvre), index),
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

FlightPlanOptimizer::MetricForInclination::MetricForInclination(
    not_null<FlightPlanOptimizer*> const optimizer,
    NavigationManœuvre manœuvre,
    int const index,
    not_null<Celestial const*> const celestial,
    not_null<std::unique_ptr<NavigationFrame const>> frame,
    Angle const& target_inclination)
    : Metric(optimizer, std::move(manœuvre), index),
      celestial_(celestial),
      frame_(std::move(frame)),
      target_inclination_(target_inclination) {}

double FlightPlanOptimizer::MetricForInclination::Evaluate(
    HomogeneousArgument const& homogeneous_argument) const {
  return Pow<2>(optimizer().EvaluateRelativeInclinationWithReplacement(
      *celestial_,
      *frame_,
      target_inclination_,
      homogeneous_argument,
      manœuvre(),
      index())) / scale_;
}

Gradient<double, FlightPlanOptimizer::HomogeneousArgument>
FlightPlanOptimizer::MetricForInclination::EvaluateGradient(
    HomogeneousArgument const& homogeneous_argument) const {
  Angle const relative_inclination =
      optimizer().EvaluateRelativeInclinationWithReplacement(
          *celestial_,
          *frame_,
          target_inclination_,
          homogeneous_argument,
          manœuvre(),
          index());
  AngleGradient const gradient =
      optimizer().Evaluate𝛁RelativeInclinationWithReplacement(
          *celestial_,
          *frame_,
          target_inclination_,
          homogeneous_argument,
          manœuvre(),
          index());
  return 2 * relative_inclination * gradient / scale_;
}

double FlightPlanOptimizer::MetricForInclination::EvaluateGateauxDerivative(
    HomogeneousArgument const& homogeneous_argument,
    Difference<HomogeneousArgument> const& homogeneous_argument_direction)
    const {
  Angle const relative_inclination =
      optimizer().EvaluateRelativeInclinationWithReplacement(
          *celestial_,
          *frame_,
          target_inclination_,
          homogeneous_argument,
          manœuvre(),
          index());
  Angle const derivative =
      optimizer().EvaluateGateauxDerivativeOfRelativeInclinationWithReplacement(
          *celestial_,
          *frame_,
          target_inclination_,
          homogeneous_argument,
          homogeneous_argument_direction,
          manœuvre(),
          index());
  return 2 * relative_inclination * derivative / scale_;
}

FlightPlanOptimizer::MetricForΔv::MetricForΔv(
    not_null<FlightPlanOptimizer*> const optimizer,
    NavigationManœuvre manœuvre,
    int const index)
    : Metric(optimizer, std::move(manœuvre), index) {}

double FlightPlanOptimizer::MetricForΔv::Evaluate(
    HomogeneousArgument const& homogeneous_argument) const {
  return UpdatedManœuvre(homogeneous_argument).Δv().Norm²() / scale_;
}

Gradient<double, FlightPlanOptimizer::HomogeneousArgument>
FlightPlanOptimizer::MetricForΔv::EvaluateGradient(
    HomogeneousArgument const& homogeneous_argument) const {
  auto const updated_Δv = UpdatedManœuvre(homogeneous_argument).Δv();
  Vector<double, Frenet<Navigation>> const ΔΔv_component_of_grad =
      2 * speed_homogeneization_factor * updated_Δv / scale_;
  auto const& ΔΔv_component_of_grad_coordinates =
      ΔΔv_component_of_grad.coordinates();
  return HomogeneousArgument({0,
                              ΔΔv_component_of_grad_coordinates.x,
                              ΔΔv_component_of_grad_coordinates.y,
                              ΔΔv_component_of_grad_coordinates.z});
}

double FlightPlanOptimizer::MetricForΔv::EvaluateGateauxDerivative(
    HomogeneousArgument const& homogeneous_argument,
    Difference<HomogeneousArgument> const& homogeneous_argument_direction)
    const {
  auto const updated_Δv = UpdatedManœuvre(homogeneous_argument).Δv();
  auto const argument_direction = Dehomogeneize(homogeneous_argument_direction);
  return 2 * InnerProduct(updated_Δv, argument_direction.ΔΔv) / scale_;
}

NavigationManœuvre FlightPlanOptimizer::MetricForΔv::UpdatedManœuvre(
    HomogeneousArgument const& homogeneous_argument) const {
  auto const& original_manœuvre = manœuvre();
  return NavigationManœuvre(
      original_manœuvre.initial_mass(),
      UpdatedBurn(homogeneous_argument, original_manœuvre));
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

FlightPlanOptimizer::MetricFactory FlightPlanOptimizer::LinearCombination(
    std::vector<MetricFactory> const& factories,
    std::vector<double> const& weights) {
  return [factories, weights](not_null<FlightPlanOptimizer*> const optimizer,
                              NavigationManœuvre manœuvre,
                              int const index) {
    return make_not_null_unique<LinearCombinationOfMetrics>(
        optimizer, std::move(manœuvre), index, factories, weights);
  };
}

FlightPlanOptimizer::MetricFactory FlightPlanOptimizer::ForCelestialCentre(
    not_null<Celestial const*> const celestial) {
  return [celestial](not_null<FlightPlanOptimizer*> const optimizer,
                     NavigationManœuvre manœuvre,
                     int const index) {
    return make_not_null_unique<MetricForCelestialCentre>(
        optimizer, std::move(manœuvre), index, celestial);
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
        optimizer, std::move(manœuvre), index, celestial, target_distance);
  };
}

FlightPlanOptimizer::MetricFactory FlightPlanOptimizer::ForInclination(
    not_null<Celestial const*> const celestial,
    std::function<not_null<std::unique_ptr<NavigationFrame>>()> frame_factory,
    Angle const& target_inclination) {
  return
      [celestial, frame_factory = std::move(frame_factory), target_inclination](
          not_null<FlightPlanOptimizer*> const optimizer,
          NavigationManœuvre manœuvre,
          int const index) {
        return make_not_null_unique<MetricForInclination>(optimizer,
                                                          std::move(manœuvre),
                                                          index,
                                                          celestial,
                                                          frame_factory(),
                                                          target_inclination);
      };
}

FlightPlanOptimizer::MetricFactory FlightPlanOptimizer::ForΔv() {
  return [](not_null<FlightPlanOptimizer*> const optimizer,
            NavigationManœuvre manœuvre,
            int const index) {
    return make_not_null_unique<MetricForΔv>(
        optimizer, std::move(manœuvre), index);
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
  // We are going to repeatedly tweak the `flight_plan_`, no point in running
  // the orbit analysers.
  flight_plan_->EnableAnalysis(/*enabled=*/false);

  // Don't reuse the computations from the previous optimization.
  cache_.clear();

  // The following is a copy, and is not affected by changes to the
  // `flight_plan_`.
  NavigationManœuvre const manœuvre = flight_plan_->GetManœuvre(index);
  auto const metric = metric_factory_(this, manœuvre, index);

  auto const status_or_solution =
      BroydenFletcherGoldfarbShanno<double, HomogeneousArgument>(
          Homogeneize(start_argument_),
          std::bind(&Metric::Evaluate, metric.get(), _1),
          std::bind(&Metric::EvaluateGradient, metric.get(), _1),
          std::bind(&Metric::EvaluateGateauxDerivative, metric.get(), _1, _2),
          Δv_tolerance / speed_homogeneization_factor);
  if (status_or_solution.ok()) {
    auto const& solution = status_or_solution.value();
    return flight_plan_->Replace(UpdatedBurn(solution, manœuvre), index);
  } else {
    return status_or_solution.status();
  }
}

DistinguishedPoints<Barycentric>::value_type
FlightPlanOptimizer::EvaluateClosestPeriapsis(
    Celestial const& celestial,
    Instant const& begin_time,
    bool const extend_if_needed) const {
  auto const& celestial_trajectory = celestial.trajectory();
  auto const& vessel_trajectory = flight_plan_->GetAllSegments();

  Length distance_at_closest_periapsis;
  std::optional<DistinguishedPoints<Barycentric>::value_type> closest_periapsis;
  for (;;) {
    DistinguishedPoints<Barycentric> apoapsides;
    DistinguishedPoints<Barycentric> periapsides;
    ComputeApsides(celestial_trajectory,
                   vessel_trajectory,
                   vessel_trajectory.lower_bound(begin_time),
                   vessel_trajectory.end(),
                   /*t_max=*/InfiniteFuture,
                   max_apsides,
                   apoapsides,
                   periapsides);
    distance_at_closest_periapsis = Infinity<Length>;
    for (auto const& periapsis : periapsides) {
      auto const& [time, degrees_of_freedom] = periapsis;
      Length const distance = (degrees_of_freedom.position() -
                               celestial_trajectory.EvaluatePosition(time))
                                  .Norm();
      if (distance < distance_at_closest_periapsis) {
        distance_at_closest_periapsis = distance;
        closest_periapsis.emplace(periapsis);
      }
    }

    auto make_distinguished_point =
        [](DiscreteTrajectory<Barycentric>::value_type const& v) {
          return DistinguishedPoints<Barycentric>::value_type{
              v.time, v.degrees_of_freedom};
        };

    // Evaluate the distance at the end of the trajectory.  If it is smaller
    // than all the periapsides, increase the length of the flight plan until it
    // isn't.
    auto const& end_point = vessel_trajectory.back();
    auto const distance_at_end =
        (end_point.degrees_of_freedom.position() -
         celestial_trajectory.EvaluatePosition(end_point.time)).Norm();
    if (distance_at_end >= distance_at_closest_periapsis) {
      break;
    } else if (!extend_if_needed) {
      return make_distinguished_point(end_point);
    }

    // Try to nudge the desired final time.  This may not succeed, in which case
    // we give up.
    auto const previous_actual_final_time = flight_plan_->actual_final_time();
    auto const new_desired_final_time = Barycentre(
        {flight_plan_->initial_time(), flight_plan_->desired_final_time()},
        {1 - flight_plan_extension_factor, flight_plan_extension_factor});
    flight_plan_->SetDesiredFinalTime(new_desired_final_time).IgnoreError();
    if (flight_plan_->actual_final_time() <= previous_actual_final_time) {
      return make_distinguished_point(vessel_trajectory.back());
    }
  }

  return closest_periapsis.value();
}

DistinguishedPoints<Barycentric>::value_type
FlightPlanOptimizer::EvaluatePeriapsisWithReplacement(
    Celestial const& celestial,
    HomogeneousArgument const& homogeneous_argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  if (auto const it = cache_.find(homogeneous_argument); it != cache_.end()) {
    return it->second;
  }

  auto const replace_status =
      flight_plan_->Replace(UpdatedBurn(homogeneous_argument, manœuvre), index);
  if (progress_callback_ != nullptr) {
    progress_callback_(*flight_plan_);
  }

  // If the burn could not be replaced, e.g., because the integrator reached its
  // maximal number of steps, evaluate the distance as best as we can, without
  // trying to be smart and extend the flight plan.  This is somewhat iffy, but
  // better than the alternative of returning an infinity, which introduces
  // discontinuities.
  auto const periapsis =
      EvaluateClosestPeriapsis(celestial,
                               manœuvre.initial_time(),
                               /*extend_if_needed=*/replace_status.ok());

  flight_plan_->Replace(manœuvre.burn(), index).IgnoreError();
  cache_.emplace(homogeneous_argument, periapsis);
  return periapsis;
}

Length FlightPlanOptimizer::EvaluateDistanceToCelestialWithReplacement(
    Celestial const& celestial,
    HomogeneousArgument const& homogeneous_argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  auto const [time, degrees_of_freedom] = EvaluatePeriapsisWithReplacement(
      celestial, homogeneous_argument, manœuvre, index);
  return (degrees_of_freedom.position() -
          celestial.trajectory().EvaluatePosition(time)).Norm();
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

Angle FlightPlanOptimizer::EvaluateRelativeInclinationWithReplacement(
    Celestial const& celestial,
    NavigationFrame const& frame,
    Angle const& target_inclination,
    HomogeneousArgument const& homogeneous_argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  auto const [time, barycentric_degrees_of_freedom] =
      EvaluatePeriapsisWithReplacement(
          celestial, homogeneous_argument, manœuvre, index);
  auto const navigation_degrees_of_freedom =
      frame.ToThisFrameAtTime(time)(barycentric_degrees_of_freedom);
  auto const r = navigation_degrees_of_freedom.position() - Navigation::origin;
  auto const v = navigation_degrees_of_freedom.velocity();
  Angle const i =
      AngleBetween(Wedge(r, v), Bivector<double, Navigation>({0, 0, 1}));
  return ReduceAngle<-π, π>(i - target_inclination);
}

FlightPlanOptimizer::AngleGradient
FlightPlanOptimizer::Evaluate𝛁RelativeInclinationWithReplacement(
    Celestial const& celestial,
    NavigationFrame const& frame,
    Angle const& target_inclination,
    HomogeneousArgument const& homogeneous_argument,
    NavigationManœuvre const& manœuvre,
    int const index) {
  auto const angle =
      EvaluateRelativeInclinationWithReplacement(celestial,
                                                 frame,
                                                 target_inclination,
                                                 homogeneous_argument,
                                                 manœuvre,
                                                 index);

  AngleGradient gradient;
  for (int k = 0; k < HomogeneousArgument::dimension; ++k) {
    HomogeneousArgument homogeneous_argument_δk = homogeneous_argument;
    homogeneous_argument_δk[k] += δ_homogeneous_argument;
    auto const angle_δk =
        EvaluateRelativeInclinationWithReplacement(celestial,
                                                   frame,
                                                   target_inclination,
                                                   homogeneous_argument_δk,
                                                   manœuvre,
                                                   index);
    gradient[k] = (angle_δk - angle) / δ_homogeneous_argument;
  }
  return gradient;
}

Angle FlightPlanOptimizer::
    EvaluateGateauxDerivativeOfRelativeInclinationWithReplacement(
        Celestial const& celestial,
        NavigationFrame const& frame,
        Angle const& target_inclination,
        HomogeneousArgument const& homogeneous_argument,
        Difference<HomogeneousArgument> const& direction_homogeneous_argument,
        NavigationManœuvre const& manœuvre,
        int const index) {
  auto const angle =
      EvaluateRelativeInclinationWithReplacement(celestial,
                                                  frame,
                                                  target_inclination,
                                                  homogeneous_argument,
                                                  manœuvre,
                                                  index);
  double const h =
      δ_homogeneous_argument / direction_homogeneous_argument.Norm();
  auto const homogeneous_argument_h =
      homogeneous_argument + h * direction_homogeneous_argument;
  auto const angle_δh =
      EvaluateRelativeInclinationWithReplacement(celestial,
                                                 frame,
                                                 target_inclination,
                                                 homogeneous_argument_h,
                                                 manœuvre,
                                                 index);
  return (angle_δh - angle) / h;
}

NavigationManœuvre::Burn FlightPlanOptimizer::UpdatedBurn(
    HomogeneousArgument const& homogeneous_argument,
    NavigationManœuvre const& manœuvre) {
  auto const argument = Dehomogeneize(homogeneous_argument);
  NavigationManœuvre::Burn burn = manœuvre.burn();
  burn.intensity = {.Δv = manœuvre.Δv() + argument.ΔΔv};
  burn.timing = {.initial_time =
                     manœuvre.initial_time() + argument.Δinitial_time};
  return burn;
}

}  // namespace internal
}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
