#pragma once

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

// A class to optimize a flight to go through or near a celestial.
class FlightPlanOptimizer {
 public:
  // The |Argument| is relative to the current properties of the burn.
  struct Argument {
    Time Δinitial_time;
    Velocity<Frenet<Navigation>> ΔΔv;
  };

  class Metric {
   public:
    virtual double Evaluate(Argument const& argument) const = 0;
    virtual Gradient<double, Argument> EvaluateGradient(
        Argument const& argument) const = 0;
    virtual double EvaluateGateauxDerivative(
        Argument const& argument,
        Difference<Argument> const& direction) const = 0;
  };

  Metric ForCelestialCentre(not_null<Celestial const*> celestial);
  Metric ForCelestialDistance(not_null<Celestial const*> celestial,
                              Length const& distance);

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
  absl::Status Optimize(Metric const& metric,
                        int index,
                        Speed const& Δv_tolerance);

 private:
  class MetricForCelestialCentre;

  // The data structure passed to the gradient descent algorithm.
  using HomogeneousArgument = FixedVector<Speed, 4>;

  // Function evaluations are very expensive, as they require integrating a
  // flight plan and finding periapsides.  We don't want do to them
  // unnecessarily.  You generally don't want to hash floats, but it's a case
  // where it's kosher.
  using EvaluationCache = absl::flat_hash_map<HomogeneousArgument, Length>;

  using LengthField = Field<Length, HomogeneousArgument>;
  using LengthGradient = Gradient<Length, HomogeneousArgument>;

  static HomogeneousArgument Homogeneize(Argument const& argument);
  static Argument Dehomogeneize(
      HomogeneousArgument const& homogeneous_argument);

  // Compute the closest periapsis of the |flight_plan| with respect to the
  // |celestial|, occurring after |begin_time|.
  Length EvaluateDistanceToCelestial(Celestial const& celestial,
                                     Instant const& begin_time,
                                     FlightPlan const& flight_plan);

  // Replaces the manœuvre at the given |index| based on the |argument|, and
  // computes the closest periapis.  Leaves the |flight_plan| unchanged.
  Length EvaluateDistanceToCelestialWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index,
      FlightPlan& flight_plan);

  // Replaces the manœuvre at the given |index| based on the |argument|, and
  // computes the gradient of the closest periapis with respect to the
  // |argument|.  Leaves the |flight_plan| unchanged.
  LengthGradient Evaluate𝛁DistanceToCelestialWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index,
      FlightPlan& flight_plan);

  Length EvaluateGateauxDerivativeOfDistanceToCelestialWithReplacement(
      Celestial const& celestial,
      HomogeneousArgument const& homogeneous_argument,
      Difference<HomogeneousArgument> const& direction_homogeneous_argument,
      NavigationManœuvre const& manœuvre,
      int index,
      FlightPlan& flight_plan);

  // Replaces the burn at the given |index| based on the |argument|.
  static absl::Status ReplaceBurn(Argument const& argument,
                                  NavigationManœuvre const& manœuvre,
                                  int index,
                                  FlightPlan& flight_plan);

  static constexpr Argument start_argument_{};
  not_null<FlightPlan*> const flight_plan_;
  ProgressCallback const progress_callback_;

  EvaluationCache cache_;

  friend bool operator==(Argument const& left, Argument const& right);
  template<typename H>
  friend H AbslHashValue(H h, Argument const& argument);
};

}  // namespace internal

using internal::FlightPlanOptimizer;

}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
