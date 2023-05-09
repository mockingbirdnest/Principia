#include "geometric_potential_plotter.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "numerics/global_optimization.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"

namespace principia {
namespace ksp_plugin {
namespace _geometric_potential_plotter {
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

GeometricPotentialPlotter::GeometricPotentialPlotter(
    not_null<Ephemeris<Barycentric>*> const ephemeris)
    : ephemeris_(ephemeris) {}

void GeometricPotentialPlotter::Interrupt() {
  plotter_ = jthread();
  // We are single-threaded here, no need to lock.
  plotter_idle_ = true;
}

void GeometricPotentialPlotter::RequestEquipotentials(
    Parameters const& parameters) {
  last_parameters_ = parameters;
  absl::MutexLock l(&lock_);
  // Only process this request if there is no analysis in progress.
  if (plotter_idle_) {
    plotter_idle_ = false;
    plotter_ = MakeStoppableThread([this, parameters]() {
      auto const status = PlotEquipotentials(parameters);
      if (!status.ok() && !absl::IsCancelled(status)) {
        LOG(ERROR) << "Error while computing " << parameters.primary->name()
                   << "-" << parameters.secondaries.front()->name()
                   << " equipotentials at " << parameters.time << ": "
                   << status;
      }
    });
  }
}

std::optional<GeometricPotentialPlotter::Parameters> const&
GeometricPotentialPlotter::last_parameters() const {
  return last_parameters_;
}

void GeometricPotentialPlotter::RefreshEquipotentials() {
  absl::MutexLock l(&lock_);
  if (next_equipotentials_.has_value()) {
    equipotentials_ = std::move(next_equipotentials_);
    next_equipotentials_.reset();
  }
}

GeometricPotentialPlotter::Equipotentials const*
GeometricPotentialPlotter::equipotentials() const {
  return equipotentials_.has_value() ? &*equipotentials_ : nullptr;
}

absl::Status GeometricPotentialPlotter::PlotEquipotentials(
    Parameters const& parameters) {
  Instant const& t = parameters.time;

  RotatingPulsatingReferenceFrame<Barycentric, Navigation> const
      reference_frame(ephemeris_, parameters.primary, parameters.secondaries);
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
      absl::MutexLock l(&lock_);
      // We don’t reset |next_equipotentials_|, so that if this is transient, we
      // keep the old ones until the problem goes away.
      plotter_idle_ = true;
      return absl::OutOfRangeError(absl::StrCat("Improper maximum ",
                                                DebugString(maximum),
                                                " at ",
                                                DebugString(arg_maximum)));
    }
    maximum_maximorum = std::max(maximum_maximorum, maximum);
  }

  Position<Barycentric> const primary_barycentric_position =
      ephemeris_->trajectory(parameters.primary)->EvaluatePosition(t);
  BarycentreCalculator<Position<Barycentric>, GravitationalParameter>
      secondary_barycentric_position;
  for (not_null secondary : parameters.secondaries) {
    secondary_barycentric_position.Add(
        ephemeris_->trajectory(secondary)->EvaluatePosition(t),
        secondary->gravitational_parameter());
  }
  Length const r =
      (secondary_barycentric_position.Get() - primary_barycentric_position)
          .Norm();

  Position<Navigation> const primary_position =
      reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
          primary_barycentric_position);
  Position<Navigation> const secondary_position =
      reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
          secondary_barycentric_position.Get());

  // TODO(egg): Consider additional wells.
  std::vector<Equipotential<Barycentric, Navigation>::Well> wells{
      {secondary_position,
       parameters.secondaries.front()->min_radius() / r * (1 * Metre)},
      {primary_position, parameters.primary->min_radius() / r * (1 * Metre)}};

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
          (4 * Sqrt(parameters.primary->gravitational_parameter() /
                    secondary_barycentric_position.weight()));

  Equipotentials equipotentials{.lines = {},
                                .parameters = parameters};

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
    equipotentials.lines.insert(equipotentials.lines.end(),
                                std::make_move_iterator(lines.begin()),
                                std::make_move_iterator(lines.end()));
  }
  if (parameters.show_l245_level) {
    for (SpecificEnergy const energy :
         {approx_l2_energy, l45_separator}) {
      auto lines = equipotential.ComputeLines(
          plane, t, arg_maximorum, wells, towards_infinity, energy);
      equipotentials.lines.insert(equipotentials.lines.end(),
                                  std::make_move_iterator(lines.begin()),
                                  std::make_move_iterator(lines.end()));
    }
  }

  absl::MutexLock l(&lock_);
  next_equipotentials_ = std::move(equipotentials);
  plotter_idle_ = true;
  return absl::OkStatus();
}

}  // namespace internal
}  // namespace _geometric_potential_plotter
}  // namespace ksp_plugin
}  // namespace principia
