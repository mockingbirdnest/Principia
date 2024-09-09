#pragma once

#include <memory>

#include "absl/synchronization/mutex.h"
#include "base/jthread.hpp"
#include "base/not_null.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/flight_plan_optimizer.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimization_driver {
namespace internal {

using namespace principia::base::_jthread;
using namespace principia::base::_not_null;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_flight_plan_optimizer;
using namespace principia::quantities::_named_quantities;

class FlightPlanOptimizationDriver {
 public:
  struct Parameters {
    int index;
    Speed Δv_tolerance;
  };

  FlightPlanOptimizationDriver(
      not_null<std::shared_ptr<FlightPlan>> const& flight_plan,
      FlightPlanOptimizer::MetricFactory metric_factory);

  virtual ~FlightPlanOptimizationDriver();

  // Returns the last flight plan evaluated by the optimizer.
  std::shared_ptr<FlightPlan> last_flight_plan() const;

  // Returns true if the optimization is done.
  bool done() const;

  // Cancels any optimization in progress.
  void Interrupt();

  // Starts an optimization with the given parameters.  Has no effect if an
  // optimization is already happening.
  void RequestOptimization(Parameters const& parameters);

  // Returns the last parameters passed to `RequestOptimization`, or nullopt if
  // `RequestOptimization` was not called.
  std::optional<Parameters> const& last_parameters() const;

  // Waits for the current optimization (if any) to complete.
  void Wait() const;

 private:
  // The progress callback of the optimizer.
  void UpdateLastFlightPlan(FlightPlan const& flight_plan);

  // The flight plan being optimized, asynchronously modified by the optimizer.
  FlightPlan flight_plan_under_optimization_;
  FlightPlanOptimizer flight_plan_optimizer_;

  mutable absl::Mutex lock_;
  jthread optimizer_;
  bool optimizer_idle_ GUARDED_BY(lock_) = true;
  std::optional<Parameters> last_parameters_ GUARDED_BY(lock_);

  // The last flight plan evaluated by the optimizer.
  not_null<std::shared_ptr<FlightPlan>> last_flight_plan_ GUARDED_BY(lock_);
};

}  // namespace internal

using internal::FlightPlanOptimizationDriver;

}  // namespace _flight_plan_optimization_driver
}  // namespace ksp_plugin
}  // namespace principia
