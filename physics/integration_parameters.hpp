#pragma once

#include "base/not_null.hpp"
#include "integrators/integrators.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_integration_parameters {

using base::not_null;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::FixedStepSizeIntegrator;
using quantities::Length;
using quantities::Speed;
using quantities::Time;

template<typename ODE>
class AdaptiveStepParameters final {
 public:
  // The |length_| and |speed_integration_tolerance|s are used to compute the
  // |tolerance_to_error_ratio| for step size control.  The number of steps is
  // limited to |max_steps|.
  AdaptiveStepParameters(AdaptiveStepSizeIntegrator<ODE> const& integrator,
                         std::int64_t max_steps,
                         Length const& length_integration_tolerance,
                         Speed const& speed_integration_tolerance);

  AdaptiveStepSizeIntegrator<ODE> const& integrator() const;
  std::int64_t max_steps() const;
  Length length_integration_tolerance() const;
  Speed speed_integration_tolerance() const;

  void set_max_steps(std::int64_t max_steps);
  void set_length_integration_tolerance(
      Length const& length_integration_tolerance);
  void set_speed_integration_tolerance(
      Speed const& speed_integration_tolerance);

  void WriteToMessage(
      not_null<serialization::AdaptiveStepParameters*> message) const;
  static AdaptiveStepParameters ReadFromMessage(
      serialization::AdaptiveStepParameters const& message);

 private:
  // This will refer to a static object returned by a factory.
  not_null<AdaptiveStepSizeIntegrator<ODE> const*> integrator_;
  std::int64_t max_steps_;
  Length length_integration_tolerance_;
  Speed speed_integration_tolerance_;
};

template<typename ODE>
class FixedStepParameters final {
 public:
  FixedStepParameters(FixedStepSizeIntegrator<ODE> const& integrator,
                      Time const& step);

  FixedStepSizeIntegrator<ODE> const& integrator() const;
  Time const& step() const;

  void WriteToMessage(
      not_null<serialization::FixedStepParameters*> message) const;
  static FixedStepParameters ReadFromMessage(
      serialization::FixedStepParameters const& message);

 private:
  // This will refer to a static object returned by a factory.
  not_null<FixedStepSizeIntegrator<ODE> const*> integrator_;
  Time step_;
};

}  // namespace internal_integration_parameters

using internal_integration_parameters::AdaptiveStepParameters;
using internal_integration_parameters::FixedStepParameters;

}  // namespace physics
}  // namespace principia

#include "physics/integration_parameters_body.hpp"
