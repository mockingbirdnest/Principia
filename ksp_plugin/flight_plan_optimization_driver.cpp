#include "ksp_plugin/flight_plan_optimization_driver.hpp"

#include <functional>
#include <memory>
#include <utility>

#include "absl/status/status.h"
#include "absl/synchronization/mutex.h"
#include "glog/logging.h"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan_optimization_driver {
namespace internal {

using std::placeholders::_1;

FlightPlanOptimizationDriver::FlightPlanOptimizationDriver(
    not_null<std::shared_ptr<FlightPlan>> const& flight_plan,
    FlightPlanOptimizer::MetricFactory metric_factory)
    : flight_plan_under_optimization_(*flight_plan),  // Copy.
      flight_plan_optimizer_(&flight_plan_under_optimization_,
                             std::move(metric_factory),
                             [this](FlightPlan const& flight_plan) {
                               UpdateLastFlightPlan(flight_plan);
                             }),
      last_flight_plan_(flight_plan) {}

FlightPlanOptimizationDriver::~FlightPlanOptimizationDriver() {
  // Ensure that we do not have a thread still running with references to the
  // members of this class when those are destroyed.
  Interrupt();
}

std::shared_ptr<FlightPlan>
FlightPlanOptimizationDriver::last_flight_plan() const {
  absl::ReaderMutexLock l(&lock_);
  return last_flight_plan_;
}

bool FlightPlanOptimizationDriver::done() const {
  absl::ReaderMutexLock l(&lock_);
  return optimizer_idle_;
}

void FlightPlanOptimizationDriver::Interrupt() {
  optimizer_ = jthread();
  // We are single-threaded here, no need to lock.
  optimizer_idle_ = true;
}

void FlightPlanOptimizationDriver::RequestOptimization(
    Parameters const& parameters) {
  // Only process this request if there is no analysis in progress.
  absl::MutexLock l(&lock_);
  if (optimizer_idle_) {
    last_parameters_ = parameters;
    optimizer_idle_ = false;
    optimizer_ = MakeStoppableThread([this, parameters]() {
      const absl::Status optimization_status = flight_plan_optimizer_.Optimize(
          parameters.index, parameters.Δv_tolerance);

      absl::MutexLock l(&lock_);
      if (optimization_status.ok()) {
        last_flight_plan_ =
            make_not_null_shared<FlightPlan>(flight_plan_under_optimization_);
        last_flight_plan_->EnableAnalysis(/*enabled=*/true);
      } else {
        LOG(WARNING) << "Optimization returned " << optimization_status;
      }
      optimizer_idle_ = true;
    });
  }
}

std::optional<FlightPlanOptimizationDriver::Parameters> const&
FlightPlanOptimizationDriver::last_parameters() const {
  absl::ReaderMutexLock l(&lock_);
  return last_parameters_;
}

void FlightPlanOptimizationDriver::Wait() const {
  absl::ReaderMutexLock l(&lock_);
  lock_.Await(absl::Condition(&optimizer_idle_));
}

void FlightPlanOptimizationDriver::UpdateLastFlightPlan(
    FlightPlan const& flight_plan) {
  absl::MutexLock l(&lock_);
  last_flight_plan_ = make_not_null_shared<FlightPlan>(flight_plan);
  // TODO(phl): This is wasteful, but otherwise, what happens if we interrupt
  // optimization?
  last_flight_plan_->EnableAnalysis(/*enabled=*/true);
}

}  // namespace internal
}  // namespace _flight_plan_optimization_driver
}  // namespace ksp_plugin
}  // namespace principia
