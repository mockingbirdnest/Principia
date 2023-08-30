#pragma once

#include <optional>

#include "absl/container/flat_hash_map.h"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/gradient_descent.hpp"
#include "physics/reference_frame.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include <functional>
#include <absl/status/status.h>

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
using namespace principia::physics::_reference_frame;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// A class to optimize a flight to go through or near a celestial.  This class
// is *not* thread-safe.
class FlightPlanOptimizer {
 public:
  // The |Argument| is relative to the current properties of the burn.
  struct Argument {
    Time Δinitial_time;
    Velocity<Frenet<Navigation>> ΔΔv;
  };

  // The data structure passed to the gradient descent algorithm.
  using HomogeneousArgument = FixedVector<Speed, 4>;

  // REMOVE BEFORE FLIGHT: Comment.
  class Metric {
   public:
    // REMOVE BEFORE FLIGHT: Ugly that we mutate the metric.
    void Initialize(not_null<FlightPlanOptimizer*> optimizer,
                    NavigationManœuvre manœuvre,
                    int index);

    virtual double Evaluate(
        HomogeneousArgument const& homogeneous_argument) const = 0;
    virtual Gradient<double, HomogeneousArgument> EvaluateGradient(
        HomogeneousArgument const& homogeneous_argument) const = 0;
    virtual double EvaluateGateauxDerivative(
        HomogeneousArgument const& homogeneous_argument,
        Difference<HomogeneousArgument> const& homogeneous_argument_direction)
        const = 0;

   protected:
    FlightPlanOptimizer& optimizer() const;
    NavigationManœuvre const& manœuvre() const;
    int index() const;

   private:
    FlightPlanOptimizer* optimizer_ = nullptr;
    std::optional<NavigationManœuvre> manœuvre_;
    std::optional<int> index_;
  };

  static Metric ForCelestialCentre(not_null<Celestial const*> celestial);
  static Metric ForCelestialDistance(not_null<Celestial const*> celestial,
                                     Length const& target_distance);

  // Called throughout the optimization to let the client know the tentative
  // state of the flight plan.
  using ProgressCallback = std::function<void(FlightPlan const&)>;

  // Constructs an optimizer for |flight_plan|.  |flight_plan| must outlive this
  // object.
  FlightPlanOptimizer(not_null<FlightPlan*> flight_plan,
                      ProgressCallback progress_callback = nullptr);

  // Optimizes the manœuvre at the given |index| to go through (or close to)
  // |celestial|.  The |Δv_tolerance| is used for the initial choice of the step
  // and for deciding when to stop, and must be small enough to not miss
  // interesting features of the trajectory, and large enough to avoid costly
  // startup steps.  Changes the flight plan passed at construction.
  // REMOVE BEFORE FLIGHT: Comment.
  absl::Status Optimize(Metric& metric,
                        int index,
                        Speed const& Δv_tolerance);

 private:
  class MetricForCelestialCentre;
  class MetricForCelestialDistance;

  // Function evaluations are very expensive, as they require integrating a
  // flight plan and finding periapsides.  We don't want do to them
  // unnecessarily.  You generally don't want to hash floats, but it's a case
  // where it's kosher.
  using EvaluationCache = absl::flat_hash_map<HomogeneousArgument, Length>;

  using LengthGradient = Gradient<Length, HomogeneousArgument>;

  static HomogeneousArgument Homogeneize(Argument const& argument);
  static Argument Dehomogeneize(
      HomogeneousArgument const& homogeneous_argument);

  // Compute the closest periapsis of the |flight_plan| with respect to the
  // |celestial|, occurring after |begin_time|.
  Length EvaluateDistanceToCelestial(Celestial const& celestial,
                                     Instant const& begin_time) const;

  // Replaces the manœuvre at the given |index| based on the |argument|, and
  // computes the closest periapis.  Leaves the |flight_plan| unchanged.
  Length EvaluateDistanceToCelestialWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index);

  // Replaces the manœuvre at the given |index| based on the |argument|, and
  // computes the gradient of the closest periapis with respect to the
  // |argument|.  Leaves the |flight_plan| unchanged.
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

  // Replaces the burn at the given |index| based on the |argument|.
  static absl::Status ReplaceBurn(Argument const& argument,
                                  NavigationManœuvre const& manœuvre,
                                  int index,
                                  FlightPlan& flight_plan);

  static constexpr Argument start_argument_{};
  not_null<FlightPlan*> const flight_plan_;
  ProgressCallback const progress_callback_;

  EvaluationCache cache_;
};

}  // namespace internal

using internal::FlightPlanOptimizer;

}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
