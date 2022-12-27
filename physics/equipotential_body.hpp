#pragma once

#include "physics/equipotential.hpp"

#include <functional>
#include <optional>
#include <tuple>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/gradient_descent.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

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
using quantities::si::Radian;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

// If the potential is below the total energy by this factor, return an empty
// equipotential line.
constexpr double energy_tolerance = 0x1p-24;

template<typename InertialFrame, typename Frame>
Equipotential<InertialFrame, Frame>::Equipotential(
    AdaptiveParameters const& adaptive_parameters,
    not_null<DynamicFrame<InertialFrame, Frame> const*> const dynamic_frame)
    : adaptive_parameters_(adaptive_parameters),
      dynamic_frame_(dynamic_frame) {}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLine(
    Plane<Frame> const& plane,
    Instant const& t,
    Position<Frame> const& position) const -> States {
  auto const binormal = plane.UnitBinormals().front();
  ODE equation{
      .compute_derivative = std::bind(
          &Equipotential::RightHandSide,
          this, binormal, position, t, _1, _2, _3)};
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

  States equipotential;
  typename AdaptiveStepSizeIntegrator<ODE>::AppendState const append_state =
      [&equipotential](SystemState const& system_state) {
        State state;
        for (auto const& y0 : std::get<0>(system_state.y)) {
          std::get<0>(state).push_back(y0.value);
        }
        for (auto const& y1 : std::get<1>(system_state.y)) {
          std::get<1>(state).push_back(y1.value);
        }
        equipotential.push_back(state);
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
    Plane<Frame> const& plane,
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) const -> States {
  // Compute the total (specific) energy.
  auto const potential_energy =
      dynamic_frame_->GeometricPotential(t, degrees_of_freedom.position());
  auto const kinetic_energy = 0.5 * degrees_of_freedom.velocity().Norm²();
  auto const total_energy = potential_energy + kinetic_energy;

  return ComputeLine(plane, t, degrees_of_freedom.position(), total_energy);
}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLine(
    Plane<Frame> const& plane,
    Instant const& t,
    Position<Frame> const& start_position,
    SpecificEnergy const& total_energy) const -> States {
  auto const lines = ComputeLines(plane, t, {start_position}, total_energy);
  CHECK_EQ(1, lines.size());
  return lines[0];
}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLines(
    Plane<Frame> const& plane,
    Instant const& t,
    std::vector<Position<Frame>> const& start_positions,
    SpecificEnergy const& total_energy) const -> std::vector<States> {
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
    return Projection(
        -2 * (dynamic_frame_->GeometricPotential(t, position) - total_energy) *
            dynamic_frame_->RotationFreeGeometricAccelerationAtRest(t,
                                                                    position),
        plane);
  };

  std::vector<States> lines;
  for (auto const& start_position : start_positions) {
    // Compute the winding number of every line already found with respect to
    // |start_position|.  If any line "turns around" that position, we don't
    // need to compute a new equipotential, it would just duplicate one we
    // already have.
    bool must_compute_line = true;
    for (auto const& line : lines) {
      std::vector<Position<Frame>> positions;
      for (auto const& state : line) {
        auto const& [positions_of_state, _] = state;
        positions.push_back(positions_of_state.front());
      }
      std::int64_t const winding_number =
          WindingNumber(plane, start_position, positions);
      if (winding_number > 0) {
        must_compute_line = false;
        break;
      }
    }
    if (!must_compute_line) {
      continue;
    }

    // Do the gradient descent to find a point on the equipotential having the
    // total energy.
    // NOTE(phl): Unclear if |length_integration_tolerance| is the right thing
    // to use below.
    auto const equipotential_position =
        BroydenFletcherGoldfarbShanno<Square<SpecificEnergy>, Position<Frame>>(
            start_position,
            f,
            grad_f,
            adaptive_parameters_.length_integration_tolerance());
    CHECK(equipotential_position.has_value());

    // The BFGS algorithm will put us at the minimum of f, but that may be a
    // point that has (significantly) less energy that our total energy.  No
    // point in building a line in that case.
    if (dynamic_frame_->GeometricPotential(t, equipotential_position.value()) <
        total_energy - Abs(total_energy) * energy_tolerance) {
      lines.push_back(States{});
      continue;
    }

    // Compute that equipotential.
    lines.push_back(ComputeLine(plane, t, equipotential_position.value()));
  }

  return lines;
}

template<typename InertialFrame, typename Frame>
absl::Status Equipotential<InertialFrame, Frame>::RightHandSide(
    Bivector<double, Frame> const& binormal,
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
      Normalize(binormal * dVǀᵧ₍ₛ₎) * characteristic_length_;

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

template<typename InertialFrame, typename Frame>
std::int64_t Equipotential<InertialFrame, Frame>::WindingNumber(
    Plane<Frame> const& plane,
    Position<Frame> const& position,
    std::vector<Position<Frame>> const& line) const {
  auto const binormal = plane.UnitBinormals().front();
  Angle angle;
  int previous_i = line.size() - 1;
  for (int i = 0; i < line.size(); ++i) {
    auto const& previous_point = line[previous_i];
    auto const& point = line[i];
    angle += OrientedAngleBetween(previous_point - position,
                                  point - position,
                                  binormal);
    previous_i = i;
  }
  return static_cast<std::int64_t>(std::round(Abs(angle) / (2 * π * Radian)));
}

}  // namespace internal_equipotential
}  // namespace physics
}  // namespace principia
