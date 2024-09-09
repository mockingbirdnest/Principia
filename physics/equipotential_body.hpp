#pragma once

#include "physics/equipotential.hpp"

#include <functional>
#include <optional>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/gradient_descent.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace physics {
namespace _equipotential {
namespace internal {

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_gradient_descent;
using namespace principia::numerics::_root_finders;
using namespace principia::quantities::_elementary_functions;

// If the potential is below the total energy by this factor, return an empty
// equipotential line.
constexpr double energy_tolerance = 0x1p-24;

template<typename InertialFrame, typename Frame>
Equipotential<InertialFrame, Frame>::Equipotential(
    AdaptiveParameters const& adaptive_parameters,
    not_null<ReferenceFrame<InertialFrame, Frame> const*> const reference_frame,
    Length const& characteristic_length)
    : adaptive_parameters_(adaptive_parameters),
      reference_frame_(reference_frame),
      characteristic_length_(characteristic_length) {}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLine(
    Plane<Frame> const& plane,
    Instant const& t,
    Position<Frame> const& position) const -> Line {
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

  Line equipotential;
  typename AdaptiveStepSizeIntegrator<ODE>::AppendState const append_state =
      [this, &equipotential, binormal, position, t](State const& state) {
        auto const& [double_q, double_β] = state.y;
        DependentVariables values;
        auto& [q, β] = values;
        q = double_q.value;
        β = double_β.value;
        DependentVariableDerivatives derivatives;
        RightHandSide(binormal, position, t, state.s.value, values, derivatives)
            .IgnoreError();
        auto const& [qʹ, βʹ] = derivatives;
        CHECK_OK(equipotential.Append(
            Instant() +
                state.s.value * reinterpret_independent_variable_as_time,
            {q, qʹ / reinterpret_independent_variable_as_time}));
      };

  auto const tolerance_to_error_ratio =
      std::bind(&Equipotential::ToleranceToErrorRatio, this, _1, _2, _3);

  auto const instance = adaptive_parameters_.integrator().NewInstance(
      problem, append_state, tolerance_to_error_ratio, integrator_parameters);
  auto status = instance->Solve(s_final_);

  return equipotential;
}

template<typename InertialFrame, typename Frame>
auto Equipotential<InertialFrame, Frame>::ComputeLines(
    Plane<Frame> const& plane,
    Instant const& t,
    std::vector<Position<Frame>> const& peaks,
    std::vector<Well> const& wells,
    std::function<Position<Frame>(Position<Frame>)> towards_infinity,
    SpecificEnergy const& energy) const -> Lines {
  using WellIterator = typename std::vector<Well>::const_iterator;

  // A `PeakDelineation` represents:
  // 1. the set of wells that are not yet delineated from a peak by
  //    equipotentials already computed;
  // 2. whether the well is delineated from the “well at infinity”.
  struct PeakDelineation {
    std::set<WellIterator> indistinct_wells;
    bool delineated_from_infinity = false;
  };

  // `peak_delineations[i]` corresponds to `peaks[i]`.
  std::vector<PeakDelineation> peak_delineations(peaks.size());
  for (auto& delineation : peak_delineations) {
    for (auto it = wells.begin(); it != wells.end(); ++it) {
      delineation.indistinct_wells.insert(it);
    }
  }

  Lines lines;
  for (int i = 0; i < peaks.size(); ++i) {
    auto const& delineation = peak_delineations[i];
    Position<Frame> const& peak = peaks[i];

    // Ignore `peak` if it is below `energy`.
    if (reference_frame_->GeometricPotential(t, peak) < energy) {
      continue;
    }

    while (!delineation.indistinct_wells.empty() ||
           !delineation.delineated_from_infinity) {
      std::optional<WellIterator> expected_delineated_well;
      bool expect_delineation_from_infinity = false;
      if (!delineation.indistinct_wells.empty()) {
        // Try to delineate `peak` from the first of its `indistinct_wells`.
        expected_delineated_well = *delineation.indistinct_wells.begin();
        Well const well = **expected_delineated_well;
        Length const r = (peak - well.position).Norm();
        if (reference_frame_->GeometricPotential(
                t,
                Barycentre({peak, well.position},
                           {well.radius, r - well.radius})) >=
            energy) {
          // The point at the edge of the well in the direction of the peak is
          // above the energy; this should not happen (the edge of the well
          // should be close enough to the singularity to be below any
          // interesting energy).
          // Give up on separating the peak from the well.
          // TODO(phl): This happens when we find the peak at the centre of the
          // Earth.
          peak_delineations[i].indistinct_wells.erase(
              *expected_delineated_well);
          continue;
        }
        // Look for a point on the equipotential along the line between the peak
        // and the edge of the well.
        Length const x = Brent(
            [&](Length const& x) {
              return reference_frame_->GeometricPotential(
                         t, Barycentre({peak, well.position}, {x, r - x})) -
                     energy;
            },
            well.radius,
            r);
        Position<Frame> const equipotential_position =
            Barycentre({peak, well.position}, {x, r - x});
        lines.push_back(ComputeLine(plane, t, equipotential_position));
      } else {
        // Try to delineate `peak` from the well at infinity; this works as for
        // an actual well, but instead of picking the point on the edge of the
        // well in the direction of the peak we generate a far away point based
        // on the peak (corresponding to a point on the edge of the well at
        // infinity).
        expect_delineation_from_infinity = true;
        Position<Frame> const far_away = towards_infinity(peak);
        if (reference_frame_->GeometricPotential(t, far_away) >= energy) {
          // The far away point is too high in the potential, presumably not far
          // enough.  Give up on separating this peak from infinity.
          peak_delineations[i].delineated_from_infinity = true;
          continue;
        }
        double const x = Brent(
            [&](double const& x) {
              return reference_frame_->GeometricPotential(
                         t, Barycentre({peak, far_away}, {x, 1 - x})) -
                     energy;
            },
            0.0,
            1.0);
        Position<Frame> const equipotential_position =
            Barycentre({peak, far_away}, {x, 1 - x});
        lines.push_back(ComputeLine(plane, t, equipotential_position));
      }
      std::vector<Position<Frame>> positions;
      for (auto const& [s, dof] : lines.back()) {
        positions.push_back(dof.position());
      }

      // Figure out whether the equipotential introduces new delineations.
      std::set<WellIterator> enclosed_wells;
      for (auto it = wells.begin(); it != wells.end(); ++it) {
        std::int64_t const winding_number =
            WindingNumber(plane, it->position, positions);
        if (winding_number > 0) {
          enclosed_wells.insert(it);
        }
      }
      for (int j = 0; j < peaks.size(); ++j) {
        bool const peak_j_enclosed =
            WindingNumber(plane, peaks[j], positions) > 0;
        peak_delineations[j].delineated_from_infinity |= peak_j_enclosed;
        for (auto it = peak_delineations[j].indistinct_wells.begin();
             it != peak_delineations[j].indistinct_wells.end();) {
          if (enclosed_wells.contains(*it) != peak_j_enclosed) {
            it = peak_delineations[j].indistinct_wells.erase(it);
          } else {
            ++it;
          }
        }
        if (j == i) {
          if (expected_delineated_well.has_value() &&
              peak_delineations[i].indistinct_wells.contains(
                  *expected_delineated_well)) {
            peak_delineations[i].indistinct_wells.erase(
                *expected_delineated_well);
          }
          if (expect_delineation_from_infinity &&
              !peak_delineations[i].delineated_from_infinity) {
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
  auto const& [γₛ, β] = values;
  // First state variable.
  auto const dVǀᵧ₍ₛ₎ =
      reference_frame_->RotationFreeGeometricAccelerationAtRest(t, γₛ);
  Displacement<Frame> const γʹ =
      Normalize(binormal * dVǀᵧ₍ₛ₎) * characteristic_length_;

  // Second state variable.
  auto const& γ₀ = position;
  double const βʹ = s == s_initial_ ? 0
                                    : Pow<2>(characteristic_length_) *
                                          (s - s_initial_) / (γₛ - γ₀).Norm²();

  derivatives = {γʹ, βʹ};

  return β > β_max_ ? absl::AbortedError("β reached max") : absl::OkStatus();
}

template<typename InertialFrame, typename Frame>
double Equipotential<InertialFrame, Frame>::ToleranceToErrorRatio(
    IndependentVariableDifference const current_s_step,
    State const& /*state*/,
    typename State::Error const& error) const {
  Length const max_length_error = std::get<0>(error).Norm();
  double const max_braking_error = Abs(std::get<1>(error));
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

}  // namespace internal
}  // namespace _equipotential
}  // namespace physics
}  // namespace principia
