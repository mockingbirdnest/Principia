#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/gradient_descent.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/reference_frame.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimizer {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_frames;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_gradient_descent;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_reference_frame;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// A class to optimize a flight to go through or near a celestial.  This class
// is *not* thread-safe.
class FlightPlanOptimizer {
 public:
  // A point in the phase space on which optimization happens.  It represents a
  // change in a burn, with respect to the pre-optimization value of that burn.
  struct Argument {
    Time Δinitial_time;
    Velocity<Frenet<Navigation>> ΔΔv;
  };

  // For the gradient descent algorithm, `Argument` is transformed into a
  // homogeneous array of numbers.
  using HomogeneousArgument = FixedVector<double, 4>;

  // These functions convert between the two representations of `Argument`.
  static HomogeneousArgument Homogeneize(Argument const& argument);
  static Argument Dehomogeneize(
      HomogeneousArgument const& homogeneous_argument);

  // A metric is an algorithm used to evaluate the quality of a solution to the
  // minimization problem.
  class Metric {
   public:
    virtual ~Metric() = default;

    virtual double Evaluate(
        HomogeneousArgument const& homogeneous_argument) const = 0;
    virtual Gradient<double, HomogeneousArgument> EvaluateGradient(
        HomogeneousArgument const& homogeneous_argument) const = 0;
    virtual double EvaluateGateauxDerivative(
        HomogeneousArgument const& homogeneous_argument,
        Difference<HomogeneousArgument> const& homogeneous_argument_direction)
        const = 0;

   protected:
    Metric(not_null<FlightPlanOptimizer*> optimizer,
           NavigationManœuvre manœuvre,
           int index);

    FlightPlanOptimizer& optimizer() const;
    NavigationManœuvre const& manœuvre() const;
    int index() const;

   private:
    not_null<FlightPlanOptimizer*> const optimizer_;
    NavigationManœuvre const manœuvre_;
    int const index_;
  };

  // A metric factory is passed at construction of the optimizer and then used
  // to construct actual metric objects during optimization.
  using MetricFactory = std::function<not_null<std::unique_ptr<Metric>>(
      not_null<FlightPlanOptimizer*> optimizer,
      NavigationManœuvre manœuvre,
      int index)>;

  static MetricFactory LinearCombination(
      std::vector<MetricFactory> const& factories,
      std::vector<double> const& weights);

  static MetricFactory ForCelestialCentre(
      not_null<Celestial const*> celestial);
  static MetricFactory ForCelestialDistance(
      not_null<Celestial const*> celestial,
      Length const& target_distance);
  static MetricFactory ForInclination(
      not_null<Celestial const*> celestial,
      std::function<not_null<std::unique_ptr<NavigationFrame>>()> frame,
      Angle const& target_inclination);
  static MetricFactory ForΔv();

  // Called throughout the optimization to let the client know the tentative
  // state of the flight plan.
  using ProgressCallback = std::function<void(FlightPlan const&)>;

  // Constructs an optimizer for `flight_plan`.  `flight_plan` must outlive this
  // object.
  FlightPlanOptimizer(not_null<FlightPlan*> flight_plan,
                      MetricFactory metric_factory,
                      ProgressCallback progress_callback = nullptr);

  // Optimizes the manœuvre at the given `index` to minimize the metric passed
  // at construction.  The `Δv_tolerance` is used for the initial choice of the
  // step and for deciding when to stop, and must be small enough to not miss
  // interesting features of the trajectory, and large enough to avoid costly
  // startup steps.  Changes the flight plan passed at construction.
  absl::Status Optimize(int index,
                        Speed const& Δv_tolerance);

 private:
  class LinearCombinationOfMetrics;
  class MetricForCelestialCentre;
  class MetricForCelestialDistance;
  class MetricForInclination;
  class MetricForΔv;

  // Function evaluations are very expensive, as they require integrating a
  // flight plan and finding periapsides.  We don't want do to them
  // unnecessarily.  You generally don't want to hash floats, but it's a case
  // where it's kosher.
  using EvaluationCache =
      absl::flat_hash_map<HomogeneousArgument,
                          DiscreteTrajectory<Barycentric>::value_type>;

  using LengthGradient = Gradient<Length, HomogeneousArgument>;
  using AngleGradient = Gradient<Angle, HomogeneousArgument>;

  // Compute the closest periapsis of the `flight_plan` with respect to the
  // `celestial`, occurring after `begin_time`.  If `extend_if_needed` is true,
  // the flight plan is extended until its end is not the point that minimizes
  // the metric.
  DiscreteTrajectory<Barycentric>::value_type EvaluateClosestPeriapsis(
      Celestial const& celestial,
      Instant const& begin_time,
      bool extend_if_needed) const;

  // Replaces the manœuvre at the given `index` based on the `argument`, and
  // computes the closest periapis.  Leaves the `flight_plan` unchanged.
  DiscreteTrajectory<Barycentric>::value_type EvaluatePeriapsisWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  Length EvaluateDistanceToCelestialWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  // Replaces the manœuvre at the given `index` based on the `argument`, and
  // computes the gradient of the closest periapis with respect to the
  // `argument`.  Leaves the `flight_plan` unchanged.
  LengthGradient Evaluate𝛁DistanceToCelestialWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  Length EvaluateGateauxDerivativeOfDistanceToCelestialWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& direction_homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  Angle EvaluateRelativeInclinationWithReplacement(
      Celestial const& celestial,
      NavigationFrame const& frame,
      Angle const& target_inclination,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  AngleGradient Evaluate𝛁RelativeInclinationWithReplacement(
      Celestial const& celestial,
      NavigationFrame const& frame,
      Angle const& target_inclination,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  Angle EvaluateGateauxDerivativeOfRelativeInclinationWithReplacement(
      Celestial const& celestial,
      NavigationFrame const& frame,
      Angle const& target_inclination,
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& direction_homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  // Returns a burn obtained by applying the changes in `homogeneous_argument`
  // to the `manœuvre`.
  static NavigationManœuvre::Burn UpdatedBurn(
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre);

  static constexpr Argument start_argument_{};
  not_null<FlightPlan*> const flight_plan_;
  MetricFactory const metric_factory_;
  ProgressCallback const progress_callback_;

  EvaluationCache cache_;
};

}  // namespace internal

using internal::FlightPlanOptimizer;

}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
