#pragma once

#include "base/not_null.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimizer {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_flight_plan;

class FlightPlanOptimizer {
 public:
  explicit FlightPlanOptimizer(not_null<FlightPlan*> flight_plan);

  void Optimize(int index, Celestial const& goal);

 private:
  not_null<FlightPlan*> const flight_plan_;
};

}  // namespace internal

using internal::FlightPlanOptimizer;

}  // namespace _flight_plan_optimizer
}  // namespace ksp_plugin
}  // namespace principia
