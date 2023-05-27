#pragma once

#include "physics/lagrange_equipotentials.hpp"
#include "numerics/global_optimization.hpp"

namespace principia {
namespace physics {
namespace _lagrange_equipotentials {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_plane;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::numerics::_global_optimization;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

template<typename Inertial, typename RotatingPulsating>
absl::StatusOr<typename Equipotential<Inertial, RotatingPulsating>::Lines>
LagrangeEquipotentials<Inertial, RotatingPulsating>::ComputeLines(
    Parameters const& parameters) {
  Instant const& t = parameters.time;

  RotatingPulsatingReferenceFrame<Barycentric, Navigation> const
      reference_frame(ephemeris_, parameters.primaries, parameters.secondaries);
  auto const plane =
      Plane<Navigation>::OrthogonalTo(Vector<double, Navigation>({0, 0, 1}));

  auto const potential = [&reference_frame,
                          &t](Position<Navigation> const& position) {
    return reference_frame.GeometricPotential(t, position);
  };
  auto const gradient = [&reference_frame,
                         &t](Position<Navigation> const& position) {
    auto const acceleration = reference_frame.GeometricAcceleration(
        t, {position, Velocity<Navigation>{}});
    // Note the sign: the acceleration goes down the potential, we need the
    // gradient which goes up.
    return -Vector<Acceleration, Navigation>({acceleration.coordinates()[0],
                                              acceleration.coordinates()[1],
                                              Acceleration{}});
  };

  // The distance between the bodies is 1 m in the rotating-pulsating frame, and
  // the origin is at their barycentre.  A 3 m box is large enough to encompass
  // the L₄ and L₅ Lagrange points, as well as the whole ridge they sit on,
  // regardless of the mass ratio.
  const MultiLevelSingleLinkage<SpecificEnergy, Position<Navigation>, 2>::Box
      box = {.centre = Navigation::origin,
             .vertices = {
                 Displacement<Navigation>({3 * Metre, 0 * Metre, 0 * Metre}),
                 Displacement<Navigation>({0 * Metre, 3 * Metre, 0 * Metre})}};

  // This is a length in the rotating-pulsating frame, so really one billionth
  // of the distance between the bodies.
  constexpr Length characteristic_length = 1 * Nano(Metre);
  Equipotential<Barycentric, Navigation> const equipotential(
      {EmbeddedExplicitRungeKuttaIntegrator<
           DormandPrince1986RK547FC,
           Equipotential<Barycentric, Navigation>::ODE>(),
       /*max_steps=*/1000,
       /*length_integration_tolerance=*/characteristic_length},
      &reference_frame,
      characteristic_length);

  // TODO(phl): Make this interruptible.
  auto const arg_maximorum =
      MultiLevelSingleLinkage<SpecificEnergy, Position<Navigation>, 2>(
          box, potential, gradient)
          .FindGlobalMaxima(
              /*points_per_round=*/1000,
              /*number_of_rounds=*/std::nullopt,
              /*local_search_tolerance=*/1e-3 * Metre);
  SpecificEnergy maximum_maximorum = -Infinity<SpecificEnergy>;
  for (auto const& arg_maximum : arg_maximorum) {
    auto const maximum = potential(arg_maximum);
    if (!IsFinite(maximum)) {
      return absl::OutOfRangeError(absl::StrCat("Improper maximum ",
                                                DebugString(maximum),
                                                " at ",
                                                DebugString(arg_maximum)));
    }
    maximum_maximorum = std::max(maximum_maximorum, maximum);
  }

  BarycentreCalculator<Position<Barycentric>, GravitationalParameter>
      primary_barycentric_position;
  for (not_null primary : parameters.primaries) {
    primary_barycentric_position.Add(
        ephemeris_->trajectory(primary)->EvaluatePosition(t),
        primary->gravitational_parameter());
  }
  BarycentreCalculator<Position<Barycentric>, GravitationalParameter>
      secondary_barycentric_position;
  for (not_null secondary : parameters.secondaries) {
    secondary_barycentric_position.Add(
        ephemeris_->trajectory(secondary)->EvaluatePosition(t),
        secondary->gravitational_parameter());
  }
  Length const r = (secondary_barycentric_position.Get() -
                    primary_barycentric_position.Get())
                       .Norm();

  Position<Navigation> const primary_position =
      reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
          primary_barycentric_position.Get());
  Position<Navigation> const secondary_position =
      reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
          secondary_barycentric_position.Get());

  Length primary_radius;
  if (parameters.primaries.size() == 1) {
    primary_radius = parameters.primaries.front()->min_radius();
  } else {
    for (not_null const primary : parameters.primaries) {
      primary_radius =
          std::max(primary_radius,
                   (primary_barycentric_position.Get() -
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
                   (secondary_barycentric_position.Get() -
                    ephemeris_->trajectory(secondary)->EvaluatePosition(t))
                           .Norm() +
                       secondary->min_radius());
    }
  }

  // TODO(egg): Consider additional wells.
  std::vector<Equipotential<Barycentric, Navigation>::Well> wells{
      {secondary_position, secondary_radius / r * (1 * Metre)},
      {primary_position, primary_radius / r * (1 * Metre)}};

  double const arg_approx_l1 = Brent(
      [&](double const x) {
        return potential(
            Barycentre(std::pair(secondary_position, primary_position),
                       std::pair(x, 1 - x)));
      },
      0.0,
      1.0,
      std::greater<>{});
  double const arg_approx_l2 = Brent(
      [&](double x) {
        // L₂ lies a short distance away the smaller body in the direction
        // opposite the barycentre.  We look for it between the smaller body and
        // a point as far from the second body as the barycentre, but in the
        // other direction.
        return potential(Barycentre(
            std::pair(secondary_position,
                      Navigation::origin +
                          2 * (secondary_position - Navigation::origin)),
            std::pair(x, 1 - x)));
      },
      0.0,
      1.0,
      std::greater<>{});
  SpecificEnergy const approx_l1_energy =
      potential(Barycentre(std::pair(secondary_position, primary_position),
                           std::pair(arg_approx_l1, 1 - arg_approx_l1)));
  SpecificEnergy const approx_l2_energy = potential(Barycentre(
      std::pair(
          secondary_position,
          Navigation::origin + 2 * (secondary_position - Navigation::origin)),
      std::pair(arg_approx_l2, 1 - arg_approx_l2)));

  // This energy is about 42% of the way down from L₄/L₅ to L₃ with equal
  // masses, and tends towards around 35% as the masses diverge.  Drawing this
  // equipotential usually distinguishes a small region around L₄/L₅ from the
  // rest of the ridge that contains them and L₃.
  SpecificEnergy const l45_separator =
      maximum_maximorum -
      (maximum_maximorum - approx_l1_energy) /
          (4 * Sqrt(primary_barycentric_position.weight() /
                    secondary_barycentric_position.weight()));

  typename Equipotential<Inertial, RotatingPulsating>::Lines result;

  auto const towards_infinity = [](Position<Navigation> q) {
    return Navigation::origin + Normalize(q - Navigation::origin) * 3 * Metre;
  };
  for (int i = 1; i <= parameters.levels; ++i) {
    RETURN_IF_STOPPED;
    SpecificEnergy const energy =
        maximum_maximorum - i * (1.0 / parameters.l1_level *
                                 (maximum_maximorum - approx_l1_energy));
    // TODO(phl): Make this interruptible.
    auto lines = equipotential.ComputeLines(
        plane, t, arg_maximorum, wells, towards_infinity, energy);
    result.insert(equipotentials.lines.end(),
                                std::make_move_iterator(lines.begin()),
                                std::make_move_iterator(lines.end()));
  }
  if (parameters.show_l245_level) {
    for (SpecificEnergy const energy : {approx_l2_energy, l45_separator}) {
      auto lines = equipotential.ComputeLines(
          plane, t, arg_maximorum, wells, towards_infinity, energy);
      result.insert(equipotentials.lines.end(),
                                  std::make_move_iterator(lines.begin()),
                                  std::make_move_iterator(lines.end()));
    }
  }
  return result;
}

}  // namespace internal
}  // namespace _lagrange_equipotentials
}  // namespace physics
}  // namespace principia
