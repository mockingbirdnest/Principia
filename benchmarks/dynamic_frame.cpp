
// .\Release\x64\benchmarks.exe --benchmark_filter=DynamicFrame --benchmark_repetitions=5  // NOLINT(whitespace/line_length)

#include <experimental/optional>
#include <memory>
#include <utility>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/barycentric_rotating_dynamic_frame.hpp"
#include "physics/body.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/solar_system.hpp"
#include "serialization/geometry.pb.h"

// This must come last because apparently it redefines CDECL.
#include "benchmark/benchmark.h"

namespace principia {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::Exp;
using geometry::Frame;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Time;
using quantities::astronomy::EarthMass;
using quantities::astronomy::JulianYear;
using quantities::si::AstronomicalUnit;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace physics {

namespace {
const Length tolerance = 0.01 * Metre;
}  // namespace

using Rendering = Frame<serialization::Frame::TestTag,
                        serialization::Frame::TEST, false>;

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
    trajectory.Append(t0 + iΔt,
                      DegreesOfFreedom<F>(initial + displacement_i,
                                          velocity));
  }
}

// This code is derived from Plugin::RenderTrajectory.
std::vector<std::pair<Position<ksp_plugin::Barycentric>,
                      Position<ksp_plugin::Barycentric>>>
ApplyDynamicFrame(
    not_null<Body const*> const body,
    not_null<DynamicFrame<ksp_plugin::Barycentric, Rendering>*> const
        dynamic_frame,
    DiscreteTrajectory<ksp_plugin::Barycentric>::Iterator const& begin,
    DiscreteTrajectory<ksp_plugin::Barycentric>::Iterator const& end) {
  std::vector<std::pair<Position<ksp_plugin::Barycentric>,
                        Position<ksp_plugin::Barycentric>>> result;

  // Compute the trajectory in the rendering frame.
  DiscreteTrajectory<Rendering> intermediate_trajectory;
  for (auto it = begin; it != end; ++it) {
    intermediate_trajectory.Append(
        it.time(),
        dynamic_frame->ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
  }

  // Render the trajectory at current time in |Rendering|.
  Instant const& current_time = intermediate_trajectory.last().time();
  DiscreteTrajectory<Rendering>::Iterator initial_it =
      intermediate_trajectory.Begin();
  DiscreteTrajectory<Rendering>::Iterator const intermediate_end =
      intermediate_trajectory.End();
  auto to_rendering_frame_at_current_time =
      dynamic_frame->FromThisFrameAtTime(current_time).rigid_transformation();
  if (initial_it != intermediate_end) {
    for (auto final_it = initial_it;
         ++final_it, final_it != intermediate_end;
         initial_it = final_it) {
      result.emplace_back(to_rendering_frame_at_current_time(
                              initial_it.degrees_of_freedom().position()),
                          to_rendering_frame_at_current_time(
                              final_it.degrees_of_freedom().position()));
    }
  }
  return result;
}

void BM_BodyCentredNonRotatingDynamicFrame(benchmark::State& state) {
  Time const Δt = 5 * Minute;
  int const steps = state.range_x();

  SolarSystem<ksp_plugin::Barycentric> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt");
  auto const ephemeris = solar_system.MakeEphemeris(
      /*fitting_tolerance=*/5 * Milli(Metre),
      Ephemeris<ksp_plugin::Barycentric>::FixedStepParameters(
          McLachlanAtela1992Order5Optimal<Position<ksp_plugin::Barycentric>>(),
          /*step=*/45 * Minute));
  ephemeris->Prolong(solar_system.epoch() + steps * Δt);

  not_null<MassiveBody const*> const earth =
      solar_system.massive_body(*ephemeris, "Earth");

  MasslessBody probe;
  Position<ksp_plugin::Barycentric> probe_initial_position =
      ksp_plugin::Barycentric::origin + Displacement<ksp_plugin::Barycentric>(
                                     {0.5 * AstronomicalUnit,
                                      -1 * AstronomicalUnit,
                                      0 * AstronomicalUnit});
  Velocity<ksp_plugin::Barycentric> probe_velocity =
      Velocity<ksp_plugin::Barycentric>({0 * SIUnit<Speed>(),
                                  100 * Kilo(Metre) / Second,
                                  0 * SIUnit<Speed>()});
  DiscreteTrajectory<ksp_plugin::Barycentric> probe_trajectory;
  FillLinearTrajectory<ksp_plugin::Barycentric, DiscreteTrajectory>(
      probe_initial_position,
      probe_velocity,
      solar_system.epoch(),
      Δt,
      steps,
      probe_trajectory);

  BodyCentredNonRotatingDynamicFrame<ksp_plugin::Barycentric, Rendering>
      dynamic_frame(ephemeris.get(), earth);
  while (state.KeepRunning()) {
    auto v = ApplyDynamicFrame(&probe,
                               &dynamic_frame,
                               probe_trajectory.Begin(),
                               probe_trajectory.End());
  }
}

void BM_BarycentricRotatingDynamicFrame(benchmark::State& state) {
  Time const Δt = 5 * Minute;
  int const steps = state.range_x();

  SolarSystem<ksp_plugin::Barycentric> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2433282_500000000.proto.txt");
  auto const ephemeris = solar_system.MakeEphemeris(
      /*fitting_tolerance=*/5 * Milli(Metre),
      Ephemeris<ksp_plugin::Barycentric>::FixedStepParameters(
          McLachlanAtela1992Order5Optimal<Position<ksp_plugin::Barycentric>>(),
          /*step=*/45 * Minute));
  ephemeris->Prolong(solar_system.epoch() + steps * Δt);

  not_null<MassiveBody const*> const earth =
      solar_system.massive_body(*ephemeris, "Earth");
  not_null<MassiveBody const*> const venus =
      solar_system.massive_body(*ephemeris, "Venus");

  MasslessBody probe;
  Position<ksp_plugin::Barycentric> probe_initial_position =
      ksp_plugin::Barycentric::origin + Displacement<ksp_plugin::Barycentric>(
                                     {0.5 * AstronomicalUnit,
                                      -1 * AstronomicalUnit,
                                      0 * AstronomicalUnit});
  Velocity<ksp_plugin::Barycentric> probe_velocity =
      Velocity<ksp_plugin::Barycentric>({0 * SIUnit<Speed>(),
                                  100 * Kilo(Metre) / Second,
                                  0 * SIUnit<Speed>()});
  DiscreteTrajectory<ksp_plugin::Barycentric> probe_trajectory;
  FillLinearTrajectory<ksp_plugin::Barycentric, DiscreteTrajectory>(
      probe_initial_position,
      probe_velocity,
      solar_system.epoch(),
      Δt,
      steps,
      probe_trajectory);

  BarycentricRotatingDynamicFrame<ksp_plugin::Barycentric, Rendering>
      dynamic_frame(ephemeris.get(), earth, venus);
  while (state.KeepRunning()) {
    auto v = ApplyDynamicFrame(&probe,
                               &dynamic_frame,
                               probe_trajectory.Begin(),
                               probe_trajectory.End());
  }
}

int const iterations = (1000 << 10) + 1;

BENCHMARK(BM_BodyCentredNonRotatingDynamicFrame)->Arg(iterations);
BENCHMARK(BM_BarycentricRotatingDynamicFrame)->Arg(iterations);

}  // namespace physics
}  // namespace principia
