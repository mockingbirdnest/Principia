#pragma once

#include "physics/equipotential.hpp"

#include <functional>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"
#include "numerics/double_precision.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_equipotential {

using geometry::InfiniteFuture;
using geometry::Normalize;
using geometry::Vector;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaIntegrator;
using integrators::ExplicitFirstOrderOrdinaryDifferentialEquation;
using integrators::IntegrationProblem;
using integrators::methods::DormandPrince1986RK547FC;
using numerics::DoublePrecision;
using quantities::Difference;
using quantities::Speed;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;
using ::std::placeholders::_1;
using ::std::placeholders::_2;


template<typename ODE>
template<typename Frame>
Equipotential<Frame>::ODEAdaptiveStepParameters<ODE>::ODEAdaptiveStepParameters(
    AdaptiveStepSizeIntegrator<ODE> const& integrator,
    std::int64_t max_steps,
    Length const& length_integration_tolerance)
    : integrator_(&integrator),
      max_steps_(max_steps),
      length_integration_tolerance_(length_integration_tolerance) {}

template<typename Frame>
Equipotential<Frame>::Equipotential(
    AdaptiveParameters const& adaptive_parameters,
    Ephemeris<Frame> const& ephemeris)
    : adaptive_parameters_(adaptive_parameters),
      ephemeris_(&ephemeris) {}

template<typename Frame>
std::vector<Position<Frame>> Equipotential<Frame>::ComputeLine(
    Bivector<double, Frame> const& plane,
    Position<Frame> const& position,
    Instant const& t) {
  ODE equation{
      .compute_derivative =
          [&ephemeris, &plane, &t](
              IndependentVariable const& s,
              typename ODE<Frame>::State const& state,
              typename ODE<Frame>::StateVariation& state_variation) {
            Velocity<Frame> velocity;
            auto const status = RightHandSide(
                ephemeris, plane, t, std::get<0>(state).front(), velocity);
            std::get<0>(state_variation)[0] = velocity;
            return status;
          }};
  typename ODE<Frame>::SystemState initial_state(
      {{position}}, s_initial);
  IntegrationProblem<ODE<Frame>> const problem{
      .equation = std::move(equation),
      .initial_state = std::move(initial_state)};

  typename AdaptiveStepSizeIntegrator<ODE<Frame>>::Parameters const
      integrator_parameters(
          /*first_time_step=*/s_initial_step,
          /*safety_factor=*/0.9,
          /*max_steps=*/1000,
          /*last_step_is_exact=*/true);

  std::vector<Position<Frame>> equipotential;

  typename AdaptiveStepSizeIntegrator<ODE<Frame>>::AppendState const append_state =
      [&equipotential](typename ODE<Frame>::SystemState const& system_state) {
        equipotential.push_back(std::get<0>(system_state.y).front().value);
      };
  auto const tolerance_to_error_ratio =
      std::bind(ToleranceToErrorRatio<Frame>,
                _1, _2);


  auto const instance = integrator.NewInstance(
      problem, append_state, tolerance_to_error_ratio, integrator_parameters);
  auto status = instance->Solve(s_final);

  return equipotential;
}

template<typename Frame>
absl::Status Equipotential<Frame>::RightHandSide(
    Bivector<double, Frame> const& plane,
    Position<Frame> const& position,
    Instant const& t,
    IndependentVariable const& s,
    State const& state,
    StateVariation& state_variation) {
  // First state variable.
  auto const& γₛ = std::get<0>(state).front();
  auto const dVǀᵧ₍ₛ₎ =
      ephemeris.ComputeGravitationalAccelerationOnMasslessBody(position, t);
  Velocity<Frame> const γʹ = Normalize(plane * dVǀᵧ₍ₛ₎) * characteristic_speed;

  // Second state variable.
  auto const& γ₀ = position;
  double const βʹ = characteristic_speed * s / (γₛ - γ₀).Norm();

  std::get<0>(state_variation).front() = γʹ;
  std::get<1>(state_variation).front() = βʹ;

  return absl::OkStatus();
}

template<typename Frame>
double Equipotential<Frame>::ToleranceToErrorRatio(
    Difference<IndependentVariable> const& current_s_step,
    typename ODE::SystemStateError const& error) {
  Length const max_length_error = std::get<0>(error).front().Norm();
  double const max_braking_error = Abs(std::get<1>(error).front());
  return std::min(length_integration_tolerance / max_length_error,
                  max_braking_error);
}

}  // namespace internal_equipotential
}  // namespace physics
}  // namespace principia
