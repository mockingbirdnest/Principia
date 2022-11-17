#pragma once

#include "physics/equipotential.hpp"

#include <functional>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/gradient_descent.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace physics {
namespace internal_equipotential {

using geometry::Normalize;
using geometry::Displacement;
using geometry::Trivector;  // We don't use this every day.
using geometry::Vector;
using geometry::Wedge;
using integrators::IntegrationProblem;
using numerics::BroydenFletcherGoldfarbShanno;
using numerics::DoublePrecision;
using quantities::Abs;
using quantities::Frequency;
using quantities::Pow;
using quantities::SpecificEnergy;
using quantities::Square;
using quantities::Time;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

// If the potential is below the total energy by this factor, return an empty
// equipotential line.
constexpr double energy_tolerance = 0x1p-24;

// The following function projects a vector on the plane orthogonal to |plane|.
// TODO(phl): Why don't we have projections?
template<typename Scalar, typename Frame>
Vector<Scalar, Frame> ProjectedVector(Bivector<double, Frame> const& plane,
                                      Vector<Scalar, Frame> const& vector) {
  Trivector<Scalar, Frame> const projection_on_plane = Wedge(plane, vector);
  return vector - plane * projection_on_plane;
}

template<typename ODE>
ODEAdaptiveStepParameters<ODE>::ODEAdaptiveStepParameters(
    AdaptiveStepSizeIntegrator<ODE> const& integrator,
    std::int64_t const max_steps,
    Length const& length_integration_tolerance)
    : integrator_(&integrator),
      max_steps_(max_steps),
      length_integration_tolerance_(length_integration_tolerance) {}

template<typename ODE>
AdaptiveStepSizeIntegrator<ODE> const&
ODEAdaptiveStepParameters<ODE>::integrator() const {
  return *integrator_;
}

template<typename ODE>
std::int64_t ODEAdaptiveStepParameters<ODE>::max_steps() const {
  return max_steps_;
}

template<typename ODE>
Length ODEAdaptiveStepParameters<ODE>::length_integration_tolerance() const {
  return length_integration_tolerance_;
}

template<typename InertialFrame, typename Frame>
Equipotential<InertialFrame, Frame>::Equipotential(
    AdaptiveParameters const& adaptive_parameters,
    not_null<DynamicFrame<InertialFrame, Frame> const*> const dynamic_frame)
    : adaptive_parameters_(adaptive_parameters),
      dynamic_frame_(dynamic_frame) {}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLine(
    Bivector<double, Frame> const& plane,
    Instant const& t,
    Position<Frame> const& position) const -> State {
  ODE equation{
      .compute_derivative = std::bind(
          &Equipotential::RightHandSide, this, plane, position, t, _1, _2, _3)};
  SystemState initial_state(s_initial_, {{position}, {/*β=*/0}});
  IntegrationProblem<ODE> const problem{
      .equation = std::move(equation),
      .initial_state = std::move(initial_state)};

  typename AdaptiveStepSizeIntegrator<ODE>::Parameters const
      integrator_parameters(
          /*first_time_step=*/initial_s_step_,
          /*safety_factor=*/0.9,
          /*max_steps=*/adaptive_parameters_.max_steps(),
          /*last_step_is_exact=*/true);

  State equipotential;
  typename AdaptiveStepSizeIntegrator<ODE>::AppendState const append_state =
      [&equipotential](SystemState const& system_state) {
        std::get<0>(equipotential)
            .push_back(std::get<0>(system_state.y).front().value);
        std::get<1>(equipotential)
            .push_back(std::get<1>(system_state.y).front().value);
      };

  auto const tolerance_to_error_ratio =
      std::bind(&Equipotential::ToleranceToErrorRatio, this, _1, _2);

  auto const instance = adaptive_parameters_.integrator().NewInstance(
      problem, append_state, tolerance_to_error_ratio, integrator_parameters);
  auto status = instance->Solve(s_final_);

  return equipotential;
}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLine(
    Bivector<double, Frame> const& plane,
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) const -> State {
  // Compute the total (specific) energy.
  auto const potential_energy =
      dynamic_frame_->GeometricPotential(t, degrees_of_freedom.position());
  auto const kinetic_energy = 0.5 * degrees_of_freedom.velocity().Norm²();
  auto const total_energy = potential_energy + kinetic_energy;

  // The function on which we perform gradient descent is defined to have a
  // minimum at a position where the potential is equal to the total energy.
  auto const f = [this, t, total_energy](Position<Frame> const& position) {
    return Pow<2>(dynamic_frame_->GeometricPotential(t, position) -
                  total_energy);
  };

  auto const grad_f = [this, &plane, t, total_energy](
      Position<Frame> const& position) {
    // To keep the problem bidimensional we eliminate any off-plane component of
    // the gradient.
    return ProjectedVector(
        plane,
        -2 * (dynamic_frame_->GeometricPotential(t, position) - total_energy) *
            dynamic_frame_->RotationFreeGeometricAccelerationAtRest(t,
                                                                    position));
  };

  // Do the gradient descent to find a point on the equipotential having the
  // total energy.
  // NOTE(phl): Unclear if |length_integration_tolerance| is the right thing to
  // use below.
  auto const equipotential_position =
      BroydenFletcherGoldfarbShanno<Square<SpecificEnergy>, Position<Frame>>(
          degrees_of_freedom.position(),
          f,
          grad_f,
          adaptive_parameters_.length_integration_tolerance());

  // The BFGS algorithm will put us at the minimum of f, but that may be a point
  // that has (significantly) less energy that our total energy.  No point in
  // building a line in that case.
  if (dynamic_frame_->GeometricPotential(t, equipotential_position) <
      total_energy - Abs(total_energy) * energy_tolerance) {
    return State{};
  }

  // Compute that equipotential.
  return ComputeLine(plane, t, equipotential_position);
}

template<typename InertialFrame, typename Frame>
absl::Status Equipotential<InertialFrame, Frame>::RightHandSide(
    Bivector<double, Frame> const& plane,
    Position<Frame> const& position,
    Instant const& t,
    IndependentVariable const s,
    State const& state,
    StateVariation& state_variation) const {
  // First state variable.
  auto const& γₛ = std::get<0>(state).front();
  auto const dVǀᵧ₍ₛ₎ =
      dynamic_frame_->RotationFreeGeometricAccelerationAtRest(t, γₛ);
  Displacement<Frame> const γʹ =
      Normalize(plane * dVǀᵧ₍ₛ₎) * characteristic_length_;

  // Second state variable.
  double const β = std::get<1>(state).front();
  auto const& γ₀ = position;
  double const βʹ = s == s_initial_ ? 0
                                    : Pow<2>(characteristic_length_) *
                                          (s - s_initial_) / (γₛ - γ₀).Norm²();

  std::get<0>(state_variation).front() = γʹ;
  std::get<1>(state_variation).front() = βʹ;

  return β > β_max_ ? absl::AbortedError("β reached max") : absl::OkStatus();
}

template<typename InertialFrame, typename Frame>
double Equipotential<InertialFrame, Frame>::ToleranceToErrorRatio(
    IndependentVariableDifference const current_s_step,
    SystemStateError const& error) const {
  Length const max_length_error = std::get<0>(error).front().Norm();
  double const max_braking_error = Abs(std::get<1>(error).front());
  return std::min(
      adaptive_parameters_.length_integration_tolerance() / max_length_error,
      β_tolerance_ / max_braking_error);
}

}  // namespace internal_equipotential
}  // namespace physics
}  // namespace principia
