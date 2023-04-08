// .\Release\x64\benchmarks.exe --benchmark_filter=RigidReferenceFrame --benchmark_repetitions=5  // NOLINT(whitespace/line_length)

#include <memory>
#include <utility>
#include <vector>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "base/not_null.hpp"
#include "base/status_utilities.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "glog/logging.h"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/barycentric_rotating_reference_frame.hpp"
#include "physics/body.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace physics {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::ksp_plugin::_frames;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

using Rendering = Frame<struct RenderingTag>;

template<typename F, template<typename> class T>
void FillLinearTrajectory(Position<F> const& initial,
                          Velocity<F> const& velocity,
                          Instant const& t0,
                          Time const& Δt,
                          int const steps,
                          T<F>& trajectory) {
  for (int i = 0; i < steps; ++i) {
    Time const iΔt = i * Δt;
    Displacement<F> const displacement_i = velocity * iΔt;
    CHECK_OK(trajectory.Append(
        t0 + iΔt, DegreesOfFreedom<F>(initial + displacement_i, velocity)));
  }
}

// This code is derived from Plugin::RenderTrajectory.
std::vector<std::pair<Position<Barycentric>, Position<Barycentric>>>
ApplyReferenceFrame(
    not_null<Body const*> const body,
    not_null<RigidReferenceFrame<Barycentric, Rendering>*> const reference_frame,
    DiscreteTrajectory<Barycentric>::iterator const& begin,
    DiscreteTrajectory<Barycentric>::iterator const& end) {
  std::vector<std::pair<Position<Barycentric>,
                        Position<Barycentric>>> result;

  // Compute the trajectory in the rendering frame.
  DiscreteTrajectory<Rendering> intermediate_trajectory;
  for (auto it = begin; it != end; ++it) {
    auto const& [time, degrees_of_freedom] = *it;
    CHECK_OK(intermediate_trajectory.Append(
        time,
        reference_frame->ToThisFrameAtTime(time)(degrees_of_freedom)));
  }

  // Render the trajectory at current time in |Rendering|.
  Instant const& current_time = intermediate_trajectory.back().time;
  DiscreteTrajectory<Rendering>::iterator initial_it =
      intermediate_trajectory.begin();
  DiscreteTrajectory<Rendering>::iterator const intermediate_end =
      intermediate_trajectory.end();
  auto to_rendering_frame_at_current_time =
      reference_frame->FromThisFrameAtTime(current_time).rigid_transformation();
  if (initial_it != intermediate_end) {
    for (auto final_it = initial_it;
         ++final_it, final_it != intermediate_end;
         initial_it = final_it) {
      result.emplace_back(to_rendering_frame_at_current_time(
                              initial_it->degrees_of_freedom.position()),
                          to_rendering_frame_at_current_time(
                              final_it->degrees_of_freedom.position()));
    }
  }
  return result;
}

void BM_BodyCentredNonRotatingReferenceFrame(benchmark::State& state) {
  Time const Δt = 5 * Minute;
  int const steps = state.range(0);

  SolarSystem<Barycentric> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt",
      /*ignore_frame=*/true);
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<Barycentric>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order5Optimal,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*step=*/45 * Minute));
  CHECK_OK(ephemeris->Prolong(solar_system.epoch() + steps * Δt));

  not_null<MassiveBody const*> const earth =
      solar_system.massive_body(*ephemeris, "Earth");

  MasslessBody probe;
  Position<Barycentric> probe_initial_position =
      Barycentric::origin + Displacement<Barycentric>(
                                     {0.5 * AstronomicalUnit,
                                      -1 * AstronomicalUnit,
                                      0 * AstronomicalUnit});
  Velocity<Barycentric> probe_velocity =
      Velocity<Barycentric>({0 * si::Unit<Speed>,
                             100 * Kilo(Metre) / Second,
                             0 * si::Unit<Speed>});
  DiscreteTrajectory<Barycentric> probe_trajectory;
  FillLinearTrajectory<Barycentric, DiscreteTrajectory>(probe_initial_position,
                                                        probe_velocity,
                                                        solar_system.epoch(),
                                                        Δt,
                                                        steps,
                                                        probe_trajectory);

  BodyCentredNonRotatingReferenceFrame<Barycentric, Rendering>
      reference_frame(ephemeris.get(), earth);
  for (auto _ : state) {
    auto v = ApplyReferenceFrame(&probe,
                                 &reference_frame,
                                 probe_trajectory.begin(),
                                 probe_trajectory.end());
  }
}

void BM_BarycentricRotatingReferenceFrame(benchmark::State& state) {
  Time const Δt = 5 * Minute;
  int const steps = state.range(0);

  SolarSystem<Barycentric> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt",
      /*ignore_frame=*/true);
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<Barycentric>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order5Optimal,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*step=*/45 * Minute));
  CHECK_OK(ephemeris->Prolong(solar_system.epoch() + steps * Δt));

  not_null<MassiveBody const*> const earth =
      solar_system.massive_body(*ephemeris, "Earth");
  not_null<MassiveBody const*> const venus =
      solar_system.massive_body(*ephemeris, "Venus");

  MasslessBody probe;
  Position<Barycentric> probe_initial_position =
      Barycentric::origin + Displacement<Barycentric>({0.5 * AstronomicalUnit,
                                                       -1 * AstronomicalUnit,
                                                       0 * AstronomicalUnit});
  Velocity<Barycentric> probe_velocity =
      Velocity<Barycentric>({0 * si::Unit<Speed>,
                             100 * Kilo(Metre) / Second,
                             0 * si::Unit<Speed>});
  DiscreteTrajectory<Barycentric> probe_trajectory;
  FillLinearTrajectory<Barycentric, DiscreteTrajectory>(probe_initial_position,
                                                        probe_velocity,
                                                        solar_system.epoch(),
                                                        Δt,
                                                        steps,
                                                        probe_trajectory);

  BarycentricRotatingReferenceFrame<Barycentric, Rendering>
      reference_frame(ephemeris.get(), earth, venus);
  for (auto _ : state) {
    auto v = ApplyReferenceFrame(&probe,
                                 &reference_frame,
                                 probe_trajectory.begin(),
                                 probe_trajectory.end());
  }
}

int const iterations = (1000 << 10) + 1;

BENCHMARK(BM_BodyCentredNonRotatingReferenceFrame)
    ->Arg(iterations)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_BarycentricRotatingReferenceFrame)
    ->Arg(iterations)
    ->Unit(benchmark::kMillisecond);

}  // namespace physics
}  // namespace principia
