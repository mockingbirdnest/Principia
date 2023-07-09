#pragma once

#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "numerics/gradient_descent.hpp"
#include "physics/reference_frame.hpp"
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
using namespace principia::numerics::_gradient_descent;
using namespace principia::physics::_reference_frame;
using namespace principia::quantities::_quantities;

class FlightPlanOptimizer {
 public:
  explicit FlightPlanOptimizer(not_null<FlightPlan*> flight_plan);

  void Optimize(int index, Celestial const& celestial);

 private:
  struct Argument {
    Instant initial_time;
    Velocity<Frenet<Navigation>> Δv;
  };

  using LengthField = Field<Length, Argument>;
  using LengthGradient = Gradient<Length, Argument>;

  static Length EvaluateDistanceToCelestial(Celestial const& celestial,
                                            Instant const& begin_time,
                                            FlightPlan const& flight_plan);
  static Length EvaluateDistanceToCelestialWithReplacement(
      Celestial const& celestial,
      Argument const& argument,
      int index,
      FlightPlan& flight_plan);

  static LengthGradient Evaluate𝛁DistanceToCelestial(Celestial const& celestial,
                                                     Argument const& argument,
                                                     int const index,
                                                     FlightPlan& flight_plan);

  not_null<FlightPlan*> const flight_plan_;
};

}  // namespace internal

using internal::FlightPlanOptimizer;

}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
