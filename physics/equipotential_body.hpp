#pragma once

#include "physics/equipotential.hpp"

#include <functional>
#include <optional>
#include <set>
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
using integrators::InitialValueProblem;
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
    Position<Frame> const& position) const -> DependentVariables {
  auto const binormal = plane.UnitBinormals().front();
  ODE equation{
      .compute_derivative = std::bind(
          &Equipotential::RightHandSide,
          this, binormal, position, t, _1, _2, _3)};
  State initial_state(s_initial_, {{position}, {/*β=*/0}});
  InitialValueProblem<ODE> const problem{
      .equation = std::move(equation),
      .initial_state = std::move(initial_state)};

  typename AdaptiveStepSizeIntegrator<ODE>::Parameters const
      integrator_parameters(
          /*first_time_step=*/initial_s_step_,
          /*safety_factor=*/0.9,
          /*max_steps=*/adaptive_parameters_.max_steps(),
          /*last_step_is_exact=*/true);

  DependentVariables equipotential;
  typename AdaptiveStepSizeIntegrator<ODE>::AppendState const append_state =
      [&equipotential](State const& state) {
        std::get<0>(equipotential)
            .push_back(std::get<0>(state.y).front().value);
        std::get<1>(equipotential)
            .push_back(std::get<1>(state.y).front().value);
      };

  auto const tolerance_to_error_ratio =
      std::bind(&Equipotential::ToleranceToErrorRatio, this, _1, _2, _3);

  auto const instance = adaptive_parameters_.integrator().NewInstance(
      problem, append_state, tolerance_to_error_ratio, integrator_parameters);
  auto status = instance->Solve(s_final_);

  return equipotential;
}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLine(
    Plane<Frame> const& plane,
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) const
    -> DependentVariables {
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
    SpecificEnergy const& total_energy) const -> DependentVariables {
  auto const lines = ComputeLines(plane, t, {start_position}, total_energy);
  CHECK_EQ(1, lines.size());
  return lines[0];
}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLines(
    Plane<Frame> const& plane,
    Instant const& t,
    std::vector<Position<Frame>> const& start_positions,
    SpecificEnergy const& total_energy) const
    -> std::vector<DependentVariables> {
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

  std::vector<DependentVariables> lines;
  for (auto const& start_position : start_positions) {
    // Compute the winding number of every line already found with respect to
    // |start_position|.  If any line "turns around" that position, we don't
    // need to compute a new equipotential, it would just duplicate one we
    // already have.
    bool must_compute_line = true;
    for (auto const& l : lines) {
      auto const& line = std::get<0>(l);
      std::int64_t const winding_number =
          WindingNumber(plane, start_position, line);
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
      lines.push_back(DependentVariables{});
      continue;
    }
    // Compute that equipotential.
    lines.push_back(ComputeLine(plane, t, equipotential_position.value()));
  }

  return lines;
}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLines(
    Plane<Frame> const& plane,
    Instant const& t,
    std::vector<Position<Frame>> const& peaks,
    std::vector<Well> const& wells,
    std::function<Position<Frame>(Position<Frame>)> towards_infinity,
    SpecificEnergy const& energy) const -> std::vector<DependentVariables> {
  using WellIterator = typename std::vector<Well>::const_iterator;
  LOG(ERROR) << "V=" << energy;

  // A |PeakDelineation| represents:
  // 1. the set of wells that are not yet delineated from a peak by
  //    equipotentials already computed;
  // 2. whether the well is delineated from the “well at infinity”.
  struct PeakDelineation {
    std::set<WellIterator> indistinct_wells;
    bool delineated_from_infinity = false;
  };

  // |peak_delineations[i]| corresponds to |peaks[i]|.
  std::vector<PeakDelineation> peak_delineations(peaks.size());
  for (auto& delineation : peak_delineations) {
    for (auto it = wells.begin(); it != wells.end(); ++it) {
      delineation.indistinct_wells.insert(it);
    }
  }

  std::vector<DependentVariables> lines;
  for (int i = 0; i < peaks.size(); ++i) {
    auto const& delineation = peak_delineations[i];
    Position<Frame> const& peak = peaks[i];
    // Ignore |peak| if it is below |energy|.
    if (dynamic_frame_->GeometricPotential(t, peak) < energy) {
      LOG(ERROR) << "Ignoring peak " << i << " which is below energy";
      continue;
    }
    LOG(ERROR) << "Delineating peak " << i;
    while (!delineation.indistinct_wells.empty() ||
           !delineation.delineated_from_infinity) {
      std::optional<WellIterator> expected_delineated_well;
      bool expect_delineation_from_infinity = false;
      if (!delineation.indistinct_wells.empty()) {
        // Try to delineate |peak| from the first of its |indistinct_wells|.
        LOG(ERROR) << delineation.indistinct_wells.size()
                   << " wells to delineate";
        expected_delineated_well = *delineation.indistinct_wells.begin();
        Well const well = **expected_delineated_well;
        Length const r = (peak - well.position).Norm();
        if (dynamic_frame_->GeometricPotential(
                t,
                Barycentre(std::pair(peak, well.position),
                           std::pair(well.radius, r - well.radius))) >=
            energy) {
          // TODO(phl): This happens when we find the peak at the centre of the
          // Earth.
          LOG(ERROR) << "well " << *expected_delineated_well - wells.begin()
                     << " is weird";
          peak_delineations[i].indistinct_wells.erase(
              *expected_delineated_well);
          continue;
        }
        Length const x = numerics::Brent(
            [&](Length const& x) {
              return dynamic_frame_->GeometricPotential(
                         t,
                         Barycentre(std::pair(peak, well.position),
                                    std::pair(x, r - x))) -
                     energy;
            },
            well.radius,
            r);
        Position<Frame> const equipotential_position =
            Barycentre(std::pair(peak, well.position), std::pair(x, r - x));
        lines.push_back(ComputeLine(plane, t, equipotential_position));
      } else {
        // Try to delineate |peak| from the well at infinity.
        LOG(ERROR) << "Not delineated from infinity";
        expect_delineation_from_infinity = true;
        Position<Frame> const far_away = towards_infinity(peak);
        if (dynamic_frame_->GeometricPotential(t, far_away) >= energy) {
          // TODO(phl): This happens when we find the peak at the centre of the
          // Earth.
          LOG(ERROR) << "far away point is weird";
          peak_delineations[i].delineated_from_infinity = true;
          continue;
        }
        double const x = numerics::Brent(
            [&](double const& x) {
              return dynamic_frame_->GeometricPotential(
                         t,
                         Barycentre(std::pair(peak, far_away),
                                    std::pair(x, 1 - x))) -
                     energy;
            },
            0.0,
            1.0);
        Position<Frame> const equipotential_position =
            Barycentre(std::pair(peak, far_away), std::pair(x, 1 - x));
        lines.push_back(ComputeLine(plane, t, equipotential_position));
      }
      auto const& line = std::get<0>(lines.back());

      // Figure out whether the equipotential introduces new delineations.
      std::set<WellIterator> enclosed_wells;
      for (auto it = wells.begin(); it != wells.end(); ++it) {
        std::int64_t const winding_number =
            WindingNumber(plane, it->position, line);
        if (winding_number > 0) {
          enclosed_wells.insert(it);
        }
      }
      for (int j = 0; j < peaks.size(); ++j) {
        bool const peak_j_enclosed = WindingNumber(plane, peaks[j], line) > 0;
        if (peak_j_enclosed && !peak_delineations[j].delineated_from_infinity) {
          LOG(ERROR) << "line delineates peak " << j << " from infinity";
        }
        peak_delineations[j].delineated_from_infinity |= peak_j_enclosed;
        for (auto it = peak_delineations[j].indistinct_wells.begin();
             it != peak_delineations[j].indistinct_wells.end();) {
          if (enclosed_wells.contains(*it) != peak_j_enclosed) {
            LOG(ERROR) << "line delineates peak " << j << " from well "
                       << *it - wells.begin();
            it = peak_delineations[j].indistinct_wells.erase(it);
          } else {
            ++it;
          }
        }
        if (j == i) {
          if (expected_delineated_well.has_value() &&
              peak_delineations[i].indistinct_wells.contains(
                  *expected_delineated_well)) {
            LOG(ERROR) << "Failed to delineate peak " << i << " from well "
                       << *expected_delineated_well - wells.begin();
            peak_delineations[i].indistinct_wells.erase(
                *expected_delineated_well);
          }
          if (expect_delineation_from_infinity &&
              !peak_delineations[i].delineated_from_infinity) {
            LOG(ERROR) << "Failed to delineate peak " << i << " from infinity";
            peak_delineations[i].delineated_from_infinity = true;
          }
        }
      }
    }
  }

  return lines;
}

template<typename InertialFrame, typename Frame>
absl::Status Equipotential<InertialFrame, Frame>::RightHandSide(
    Bivector<double, Frame> const& binormal,
    Position<Frame> const& position,
    Instant const& t,
    IndependentVariable const s,
    DependentVariables const& values,
    DependentVariableDerivatives& derivatives) const {
  // First state variable.
  auto const& γₛ = std::get<0>(values).front();
  auto const dVǀᵧ₍ₛ₎ =
      dynamic_frame_->RotationFreeGeometricAccelerationAtRest(t, γₛ);
  Displacement<Frame> const γʹ =
      Normalize(binormal * dVǀᵧ₍ₛ₎) * characteristic_length_;

  // Second state variable.
  double const β = std::get<1>(values).front();
  auto const& γ₀ = position;
  double const βʹ = s == s_initial_ ? 0
                                    : Pow<2>(characteristic_length_) *
                                          (s - s_initial_) / (γₛ - γ₀).Norm²();

  std::get<0>(derivatives).front() = γʹ;
  std::get<1>(derivatives).front() = βʹ;

  return β > β_max_ ? absl::AbortedError("β reached max") : absl::OkStatus();
}

template<typename InertialFrame, typename Frame>
double Equipotential<InertialFrame, Frame>::ToleranceToErrorRatio(
    IndependentVariableDifference const current_s_step,
    State const& /*state*/,
    State::Error const& error) const {
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
