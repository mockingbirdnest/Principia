#pragma once

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include "physics/lagrange_equipotentials.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/plane.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"
#include "numerics/global_optimization.hpp"
#include "numerics/root_finders.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _lagrange_equipotentials {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_plane;
using namespace principia::integrators::_embedded_explicit_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::numerics::_global_optimization;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// The distance between the bodies is 1 m in the rotating-pulsating frame, and
// the origin is at their barycentre.  A 3 m box is large enough to encompass
// the L₄ and L₅ Lagrange points, as well as the whole ridge they sit on,
// regardless of the mass ratio.
constexpr Length box_side = 3 * Metre;
// This is a length in the rotating-pulsating frame, so really one billionth of
// the distance between the bodies.
constexpr Length characteristic_length = 1 * Nano(Metre);
constexpr std::int64_t max_steps = 1000;
constexpr std::int64_t points_per_round = 1000;
constexpr Length local_search_tolerance = 1e-3 * Metre;

template<typename Inertial, typename RotatingPulsating>
LagrangeEquipotentials<Inertial, RotatingPulsating>::LagrangeEquipotentials(
    not_null<Ephemeris<Inertial> const*> const ephemeris)
    : ephemeris_(ephemeris) {}

template<typename Inertial, typename RotatingPulsating>
absl::StatusOr<typename LagrangeEquipotentials<Inertial, RotatingPulsating>::
                   Equipotentials>
LagrangeEquipotentials<Inertial, RotatingPulsating>::ComputeLines(
    Parameters const& parameters) {
  Equipotentials result;
  Instant const& t = parameters.time;

  RotatingPulsatingReferenceFrame<Inertial, RotatingPulsating> const
      reference_frame(ephemeris_, parameters.primaries, parameters.secondaries);
  auto const plane = Plane<RotatingPulsating>::OrthogonalTo(
      Vector<double, RotatingPulsating>({0, 0, 1}));

  auto const potential = [&reference_frame,
                          &t](Position<RotatingPulsating> const& position) {
    return reference_frame.GeometricPotential(t, position);
  };
  auto const gradient = [&reference_frame,
                         &t](Position<RotatingPulsating> const& position) {
    auto const acceleration = reference_frame.GeometricAcceleration(
        t, {position, Velocity<RotatingPulsating>{}});
    // Note the sign: the acceleration goes down the potential, we need the
    // gradient which goes up.
    return -Vector<Acceleration, RotatingPulsating>(
        {acceleration.coordinates()[0],
         acceleration.coordinates()[1],
         Acceleration{}});
  };

  typename MultiLevelSingleLinkage<SpecificEnergy,
                                   Position<RotatingPulsating>,
                                   2>::Box const box = {
      .centre = RotatingPulsating::origin,
      .vertices = {
          Displacement<RotatingPulsating>({box_side, 0 * Metre, 0 * Metre}),
          Displacement<RotatingPulsating>({0 * Metre, box_side, 0 * Metre})}};

  Equipotential<Inertial, RotatingPulsating> const equipotential(
      {EmbeddedExplicitRungeKuttaIntegrator<
           DormandPrince1986RK547FC,
           typename Equipotential<Inertial, RotatingPulsating>::ODE>(),
       max_steps,
       /*length_integration_tolerance=*/characteristic_length},
      &reference_frame,
      characteristic_length);

  // TODO(phl): Make this interruptible.
  auto const arg_maximorum =
      MultiLevelSingleLinkage<SpecificEnergy, Position<RotatingPulsating>, 2>(
          box, potential, gradient)
          .FindGlobalMaxima(
              points_per_round,
              /*number_of_rounds=*/std::nullopt,
              local_search_tolerance);
  SpecificEnergy maximum_maximorum = -Infinity<SpecificEnergy>;
  for (auto const& arg_maximum : arg_maximorum) {
    auto const maximum = potential(arg_maximum);
    if (!IsFinite(maximum)) {
      return absl::OutOfRangeError(absl::StrCat("Singular maximum ",
                                                DebugString(maximum),
                                                " at ",
                                                DebugString(arg_maximum)));
    }
    result.maxima.emplace(maximum, arg_maximum);
    maximum_maximorum = std::max(maximum_maximorum, maximum);
  }

  BarycentreCalculator<Position<Inertial>, GravitationalParameter>
      primary_Inertial_position;
  for (not_null primary : parameters.primaries) {
    primary_Inertial_position.Add(
        ephemeris_->trajectory(primary)->EvaluatePosition(t),
        primary->gravitational_parameter());
  }
  BarycentreCalculator<Position<Inertial>, GravitationalParameter>
      secondary_Inertial_position;
  for (not_null secondary : parameters.secondaries) {
    secondary_Inertial_position.Add(
        ephemeris_->trajectory(secondary)->EvaluatePosition(t),
        secondary->gravitational_parameter());
  }
  Length const r =
      (secondary_Inertial_position.Get() - primary_Inertial_position.Get())
          .Norm();

  Position<RotatingPulsating> const primary_position =
      reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
          primary_Inertial_position.Get());
  Position<RotatingPulsating> const secondary_position =
      reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
          secondary_Inertial_position.Get());

  // For a system, the radius is the maximum distance from the barycentre over
  // the balls of radius min_radius centred at each body.
  Length primary_radius;
  if (parameters.primaries.size() == 1) {
    primary_radius = parameters.primaries.front()->min_radius();
  } else {
    for (not_null const primary : parameters.primaries) {
      primary_radius =
          std::max(primary_radius,
                   (primary_Inertial_position.Get() -
                    ephemeris_->trajectory(primary)->EvaluatePosition(t))
                           .Norm() +
                       primary->min_radius());
    }
  }

  Length secondary_radius;
  if (parameters.secondaries.size() == 1) {
    secondary_radius = parameters.secondaries.front()->min_radius();
  } else {
    for (not_null const secondary : parameters.secondaries) {
      secondary_radius =
          std::max(secondary_radius,
                   (secondary_Inertial_position.Get() -
                    ephemeris_->trajectory(secondary)->EvaluatePosition(t))
                           .Norm() +
                       secondary->min_radius());
    }
  }

  // TODO(egg): Consider additional wells.
  std::vector<typename Equipotential<Inertial, RotatingPulsating>::Well> wells{
      {secondary_position, secondary_radius / r * (1 * Metre)},
      {primary_position, primary_radius / r * (1 * Metre)}};

  // L₁ lies between the primary and the secondary, and on  that segment it is a
  // maximum of the potential.
  auto const potential_on_primary_secondary_segment = [&](double const x) {
    return potential(
        Barycentre({secondary_position, primary_position}, {x, 1 - x}));
  };
  double const arg_approx_l1 = Brent(
      potential_on_primary_secondary_segment,
      0.0,
      1.0,
      std::greater<>{});
  // L₂ lies a short distance away the smaller body in the direction opposite
  // the barycentre.  We look for it between the smaller body and a point as far
  // from the second body as the barycentre, but in the other direction.
  auto const potential_on_secondary_outward_segment = [&](double x) {
    return potential(
        Barycentre({secondary_position,
                    RotatingPulsating::origin +
                        2 * (secondary_position - RotatingPulsating::origin)},
                   {x, 1 - x}));
  };
  double const arg_approx_l2 = Brent(
      potential_on_secondary_outward_segment,
      0.0,
      1.0,
      std::greater<>{});
  SpecificEnergy const approx_l1_energy =
      potential_on_primary_secondary_segment(arg_approx_l1);
  SpecificEnergy const approx_l2_energy =
      potential_on_secondary_outward_segment(arg_approx_l2);

  // This energy is about 42% of the way down from L₄/L₅ to L₃ with equal
  // masses, and tends towards around 35% as the masses diverge.  Drawing this
  // equipotential usually distinguishes a small region around L₄/L₅ from the
  // rest of the ridge that contains them and L₃.
  SpecificEnergy const l45_separator =
      maximum_maximorum - (maximum_maximorum - approx_l1_energy) /
                              (4 * Sqrt(primary_Inertial_position.weight() /
                                        secondary_Inertial_position.weight()));

  auto const towards_infinity = [](Position<RotatingPulsating> q) {
    return RotatingPulsating::origin +
           Normalize(q - RotatingPulsating::origin) * box_side;
  };
  for (int i = 1; i <= parameters.levels; ++i) {
    RETURN_IF_STOPPED;
    SpecificEnergy const energy =
        maximum_maximorum - i * (1.0 / parameters.l1_level *
                                 (maximum_maximorum - approx_l1_energy));
    // TODO(phl): Make this interruptible.
    auto lines = equipotential.ComputeLines(
        plane, t, arg_maximorum, wells, towards_infinity, energy);
    result.lines.emplace(energy, std::move(lines));
  }
  if (parameters.show_l245_level) {
    for (SpecificEnergy const energy : {approx_l2_energy, l45_separator}) {
      auto lines = equipotential.ComputeLines(
          plane, t, arg_maximorum, wells, towards_infinity, energy);
      result.lines.emplace(energy, std::move(lines));
    }
  }
  return result;
}

}  // namespace internal
}  // namespace _lagrange_equipotentials
}  // namespace physics
}  // namespace principia
