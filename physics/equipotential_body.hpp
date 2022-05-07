#pragma once

#include <functional>

#include "absl/status/status.h"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/double_precision.hpp"
#include "physics/equipotential.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_equipotential {

using geometry::InfiniteFuture;
using geometry::Normalize;
using geometry::Vector;
using geometry::Velocity;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::EmbeddedExplicitRungeKuttaIntegrator;
using integrators::ExplicitFirstOrderOrdinaryDifferentialEquation;
using integrators::IntegrationProblem;
using integrators::methods::DormandPrince1986RK547FC;
using numerics::DoublePrecision;
using quantities::Difference;
using quantities::Length;
using quantities::Speed;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;
using ::std::placeholders::_1;
using ::std::placeholders::_2;

using IndependentVariable = Instant;

IndependentVariable const s_initial;
IndependentVariable const s_final = InfiniteFuture;
Difference<IndependentVariable> const s_initial_step = 1 * Second;
Speed const characteristic_speed = 1 * Metre / Second;

Length const length_integration_tolerance = 1 * Metre;

template<typename Frame>
using ODE = ExplicitFirstOrderOrdinaryDifferentialEquation<Position<Frame>>;

template<typename Frame>
absl::Status RightHandSide(Ephemeris<Frame> const& ephemeris,
                           Bivector<double, Frame> const& plane,
                           IndependentVariable const& t,
                           Position<Frame> const& position,
                           Velocity<Frame>& velocity) {
  auto const dVǀᵧ₍ₛ₎ =
      ephemeris.ComputeGravitationalAccelerationOnMasslessBody(position, t);
  velocity = Normalize(plane * dVǀᵧ₍ₛ₎) * characteristic_speed;
  return absl::OkStatus();
}

template<typename Frame>
double ToleranceToErrorRatio(
    Difference<IndependentVariable> const& current_s_step,
    typename ODE<Frame>::SystemStateError const& error) {
  Length const max_length_error = std::get<0>(error).front().Norm();
  return length_integration_tolerance / max_length_error;
}

template<typename Frame>
std::vector<Position<Frame>> ComputeEquipotential(
    Ephemeris<Frame> const& ephemeris,
    Bivector<double, Frame> const& plane,
    Position<Frame> const& position,
    Instant const& t) {
  static_assert(Frame::is_inertial);

  auto const& integrator =
      EmbeddedExplicitRungeKuttaIntegrator<DormandPrince1986RK547FC,
                                           Position<Frame>>();

  ODE<Frame> equation{
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

}  // namespace internal_equipotential
}  // namespace physics
}  // namespace principia
