#pragma once

#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_equipotential {

using base::not_null;
using geometry::Bivector;
using geometry::InfiniteFuture;
using geometry::Instant;
using geometry::Position;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::ExplicitFirstOrderOrdinaryDifferentialEquation;
using quantities::Acceleration;
using quantities::Difference;
using quantities::Infinity;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Second;

// TODO(phl): Similar class in Ephemeris.  Move to a common place.
template<typename ODE>
class ODEAdaptiveStepParameters final {
 public:
  ODEAdaptiveStepParameters(AdaptiveStepSizeIntegrator<ODE> const& integrator,
                            std::int64_t max_steps,
                            Length const& length_integration_tolerance);

  AdaptiveStepSizeIntegrator<ODE> const& integrator() const;
  std::int64_t max_steps() const;
  Length length_integration_tolerance() const;

 private:
  // This will refer to a static object returned by a factory.
  not_null<AdaptiveStepSizeIntegrator<ODE> const*> integrator_;
  std::int64_t max_steps_;
  Length length_integration_tolerance_;
};

template<typename InertialFrame, typename Frame>
class Equipotential {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  // Equipotential lines are parameterized by this dimentionless quantity.
  using IndependentVariable = double;

  // The first state variable is a point of an equipotential.  The second state
  // variable is the (scale-free) intensity of the braking.
  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                     Position<Frame>,
                                                     double>;
  using AdaptiveParameters = ODEAdaptiveStepParameters<ODE>;

  Equipotential(
      AdaptiveParameters const& adaptive_parameters,
      not_null<DynamicFrame<InertialFrame, Frame> const*> dynamic_frame);

  // Computes an equipotential line going through the given point.
  typename ODE::State ComputeLine(Bivector<double, Frame> const& plane,
                                  Instant const& t,
                                  Position<Frame> const& position) const;

  // Computes an equipotential line for the total energy determined by the
  // |degrees_of_freedom|.
  typename ODE::State ComputeLine(
      Bivector<double, Frame> const& plane,
      Instant const& t,
      DegreesOfFreedom<Frame> const& degrees_of_freedom) const;

 private:
  using IndependentVariableDifference =
      typename ODE::IndependentVariableDifference;
  using State = typename ODE::State;
  using StateVariation = typename ODE::StateVariation;
  using SystemState = typename ODE::SystemState;
  using SystemStateError = typename ODE::SystemStateError;

  static constexpr IndependentVariable const s_initial_ = 0;
  static constexpr IndependentVariable const s_final_ =
      Infinity<IndependentVariable>;
  static constexpr IndependentVariableDifference const initial_s_step_ = 1;
  static constexpr Length const characteristic_length_ = 1 * Metre;

  // TODO(phl): One or both of these values should probably be a parameter.
  static constexpr double β_max_ = 1e6;
  static constexpr double β_tolerance_ = 1;

  absl::Status RightHandSide(Bivector<double, Frame> const& plane,
                             Position<Frame> const& position,
                             Instant const& t,
                             IndependentVariable s,
                             State const& state,
                             StateVariation& state_variation) const;

  double ToleranceToErrorRatio(IndependentVariableDifference current_s_step,
                               SystemStateError const& error) const;

  AdaptiveParameters const& adaptive_parameters_;
  not_null<DynamicFrame<InertialFrame, Frame> const*> const dynamic_frame_;
};

}  // namespace internal_equipotential

using internal_equipotential::Equipotential;

}  // namespace physics
}  // namespace principia

#include "physics/equipotential_body.hpp"
