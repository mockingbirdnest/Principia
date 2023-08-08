#pragma once

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
  // Constructs an optimizer for |flight_plan|.  |flight_plan| must outlive this
  // object.
  explicit FlightPlanOptimizer(not_null<FlightPlan*> flight_plan);

  // Optimizes the manœuvre at the given |index| to go through (or close to)
  // |celestial|.  The |Δv_tolerance| is used for the initial choice of the step
  // and for deciding when to stop, and must be small enough to not miss
  // interesting features of the trajectory, and large enough to avoid costly
  // startup steps.  Changes the flight plan passed at construction.
  absl::Status Optimize(int index,
                        Celestial const& celestial,
                        Speed const& Δv_tolerance);

  //TODO(phl)comment
  absl::Status Optimize(int index,
                        Celestial const& celestial,
                        Length const& distance,
                        Speed const& Δv_tolerance);

 private:
  // The |Argument| is relative to the current properties of the burn.
  struct Argument {
    Time Δinitial_time;
    Velocity<Frenet<Navigation>> ΔΔv;
  };

  // The data structure passed to the gradient descent algorithm.
  using HomogeneousArgument = FixedVector<Speed, 4>;

  using LengthField = Field<Length, HomogeneousArgument>;
  using LengthGradient = Gradient<Length, HomogeneousArgument>;

  static HomogeneousArgument Homogeneize(Argument const& argument);
  static Argument Dehomogeneize(
      HomogeneousArgument const& homogeneous_argument);

  // Compute the closest periapsis of the |flight_plan| with respect to the
  // |celestial|, occurring after |begin_time|.
  static Length EvaluateDistanceToCelestial(Celestial const& celestial,
                                            Instant const& begin_time,
                                            FlightPlan const& flight_plan);

  // Replaces the manœuvre at the given |index| based on the |argument|, and
  // computes the closest periapis.  Leaves the |flight_plan| unchanged.
  static Length EvaluateDistanceToCelestialWithReplacement(
      Celestial const& celestial,
      Argument const& argument,
      NavigationManœuvre const& manœuvre,
      int index,
      FlightPlan& flight_plan);

  // Replaces the manœuvre at the given |index| based on the |argument|, and
  // computes the gradient of the closest periapis with respect to the
  // |argument|.  Leaves the |flight_plan| unchanged.
  static LengthGradient Evaluate𝛁DistanceToCelestialWithReplacement(
      Celestial const& celestial,
      Argument const& argument,
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
};

}  // namespace internal

using internal::FlightPlanOptimizer;

}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
