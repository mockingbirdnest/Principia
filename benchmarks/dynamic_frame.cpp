
// .\Release\benchmarks.exe --benchmark_filter=DynamicFrame --benchmark_repetitions=5          // NOLINT(whitespace/line_length)

#include <experimental/optional>
#include <memory>
#include <utility>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "quantities/astronomy.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "physics/barycentric_rotating_dynamic_frame.hpp"
#include "physics/body.hpp"
#include "physics/body_centered_non_rotating_dynamic_frame.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "serialization/geometry.pb.h"

// This must come last because apparently it redefines CDECL.
#include "benchmark/benchmark.h"

namespace principia {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::Exp;
using geometry::Frame;
using geometry::Position;
using geometry::Velocity;
using quantities::AngularFrequency;
using quantities::SIUnit;
using quantities::Time;
using physics::BarycentricRotatingDynamicFrame;
using physics::Body;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::DynamicFrame;
using physics::MassiveBody;
using physics::MasslessBody;
using quantities::astronomy::EarthMass;
using quantities::astronomy::JulianYear;
using quantities::si::AstronomicalUnit;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

namespace benchmarks {

namespace {
const Length kTolerance = 0.01 * Metre;
}  // namespace

using World1 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST1, true>;
using World2 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST2, false>;

template<typename F, template<typename F> class T>
void FillCircularTrajectory(Position<F> const& center,
                            Position<F> const& initial,
                            AngularVelocity<F> const& angular_velocity,
                            Time const& Δt,
                            int const steps,
                            not_null<T<F>*> const trajectory) {
  Displacement<F> const radius = initial - center;
  for (int i = 0; i < steps; ++i) {
    Time const t_i = i * Δt;
    Displacement<F> const displacement_i = Exp(angular_velocity * t_i)(radius);
    Velocity<F> const velocity_i = angular_velocity * displacement_i / Radian;
    trajectory->Append(Instant(t_i),
                       DegreesOfFreedom<F>(initial + displacement_i,
                                           velocity_i));
  }
}

template<typename F, template<typename F> class T>
void FillLinearTrajectory(Position<F> const& initial,
                          Velocity<F> const& velocity,
                          Time const& Δt,
                          int const steps,
                          not_null<T<F>*> const trajectory) {
  for (int i = 0; i < steps; ++i) {
    Time const t_i = i * Δt;
    Displacement<F> const displacement_i = velocity * t_i;
    trajectory->Append(Instant(t_i),
                       DegreesOfFreedom<F>(initial + displacement_i,
                                           velocity));
  }
}

// This code is derived from Plugin::RenderTrajectory.
std::vector<std::pair<Position<World1>,
                      Position<World1>>> ApplyDynamicFrame(
    not_null<Body const*> const body,
    not_null<DynamicFrame<World1, World2>*> const dynamic_frame,
    DiscreteTrajectory<World1>::TransformingIterator<World2> const& actual_it) {
  std::vector<std::pair<Position<World1>,
                        Position<World1>>> result;

  // First build the trajectory resulting from the first transform.
  DiscreteTrajectory<World2> intermediate_trajectory;
  for (auto it = actual_it; !it.at_end(); ++it) {
    intermediate_trajectory.Append(it.time(), it.degrees_of_freedom());
  }

  // Then build the final result using the second transform.
  std::experimental::optional<Position<World1>> last_position;
  for (auto it = transforms->second(intermediate_trajectory.last().time(),
                                    intermediate_trajectory);
       !it.at_end();
       ++it) {
    Position<World1> const& position = it.degrees_of_freedom().position();
    if (last_position) {
      result.emplace_back(*last_position, position);
    }
    last_position = position;
  }
  return result;
}

void BM_BodyCentredNonRotatingDynamicFrame(
    benchmark::State& state) {  // NOLINT(runtime/references)

  Time const Δt = 1 * Hour;
  int const steps = state.range_x();

  MassiveBody earth(EarthMass);
  Position<World1> center = World1::origin;
  Position<World1> earth_initial_position =
      World1::origin + Displacement<World1>({1 * AstronomicalUnit,
                                             0 * AstronomicalUnit,
                                             0 * AstronomicalUnit});
  AngularVelocity<World1> earth_angular_velocity =
      AngularVelocity<World1>({0 * SIUnit<AngularFrequency>(),
                               0 * SIUnit<AngularFrequency>(),
                               2 * π * Radian / JulianYear});
  ContinuousTrajectory<World1> earth_trajectory(Δt, kTolerance);
  FillCircularTrajectory<World1, ContinuousTrajectory>(center,
                                                       earth_initial_position,
                                                       earth_angular_velocity,
                                                       Δt,
                                                       steps,
                                                       &earth_trajectory);

  MasslessBody probe;
  Position<World1> probe_initial_position =
      World1::origin + Displacement<World1>({0.5 * AstronomicalUnit,
                                             -1 * AstronomicalUnit,
                                             0 * AstronomicalUnit});
  Velocity<World1> probe_velocity =
      Velocity<World1>({0 * SIUnit<Speed>(),
                        100 * Kilo(Metre) / Second,
                        0 * SIUnit<Speed>()});
  DiscreteTrajectory<World1> probe_trajectory;
  FillLinearTrajectory<World1, DiscreteTrajectory>(probe_initial_position,
                                                   probe_velocity,
                                                   Δt,
                                                   steps,
                                                   &probe_trajectory);

  BodyCentredNonRotatingDynamicFrame<World1, World2> dynamic_frame(
      earth, earth_trajectory, earth_trajectory);
  while (state.KeepRunning()) {
    auto v = ApplyDynamicFrame(&probe,
                                &dynamic_frame,
                            transforms->first(probe_trajectory));
  }
}

void BM_BarycentricRotatingDynamicFrame(
    benchmark::State& state) {  // NOLINT(runtime/references)

  Time const Δt = 1 * Hour;
  int const steps = state.range_x();

  MassiveBody earth(EarthMass);
  Position<World1> earth_center = World1::origin;
  Position<World1> earth_initial_position =
      World1::origin + Displacement<World1>({1 * AstronomicalUnit,
                                             0 * AstronomicalUnit,
                                             0 * AstronomicalUnit});
  AngularVelocity<World1> earth_angular_velocity =
      AngularVelocity<World1>({0 * SIUnit<AngularFrequency>(),
                               0 * SIUnit<AngularFrequency>(),
                               2 * π * Radian / JulianYear});
  ContinuousTrajectory<World1> earth_trajectory(Δt, kTolerance);
  FillCircularTrajectory<World1, ContinuousTrajectory>(earth_center,
                                                       earth_initial_position,
                                                       earth_angular_velocity,
                                                       Δt,
                                                       steps,
                                                       &earth_trajectory);

  MassiveBody thera(EarthMass);
  Position<World1> thera_center =
      World1::origin + Displacement<World1>({2 * AstronomicalUnit,
                                             0 * AstronomicalUnit,
                                             0 * AstronomicalUnit});
  Position<World1> thera_initial_position =
      World1::origin + Displacement<World1>({-0.5 * AstronomicalUnit,
                                             0 * AstronomicalUnit,
                                             0 * AstronomicalUnit});
  AngularVelocity<World1> thera_angular_velocity =
      AngularVelocity<World1>({0 * SIUnit<AngularFrequency>(),
                               0 * SIUnit<AngularFrequency>(),
                               6 * Radian / JulianYear});
  ContinuousTrajectory<World1> thera_trajectory(Δt, kTolerance);
  FillCircularTrajectory<World1, ContinuousTrajectory>(thera_center,
                                                       thera_initial_position,
                                                       thera_angular_velocity,
                                                       Δt,
                                                       steps,
                                                       &thera_trajectory);

  MasslessBody probe;
  Position<World1> probe_initial_position =
      World1::origin + Displacement<World1>({0.5 * AstronomicalUnit,
                                             -1 * AstronomicalUnit,
                                             0 * AstronomicalUnit});
  Velocity<World1> probe_velocity =
      Velocity<World1>({0 * SIUnit<Speed>(),
                        100 * Kilo(Metre) / Second,
                        0 * SIUnit<Speed>()});
  DiscreteTrajectory<World1> probe_trajectory;
  FillLinearTrajectory<World1, DiscreteTrajectory>(probe_initial_position,
                                                   probe_velocity,
                                                   Δt,
                                                   steps,
                                                   &probe_trajectory);

  BarycentricRotatingDynamicFrame<World1, World2> dynamic_frame(earth,
                                                                earth_trajectory,
                                                                earth_trajectory,
                                                                thera,
                                                                thera_trajectory,
                                                                thera_trajectory);
  while (state.KeepRunning()) {
    auto v = ApplyDynamicFrame(&probe,
                                 &dynamic_frame,
                              transforms->first(probe_trajectory));
  }
}

int const kIter = (1000 << 10) + 1;

BENCHMARK(BM_BodyCentredNonRotatingDynamicFrame)->Arg(kIter);
BENCHMARK(BM_BarycentricRotatingDynamicFrame)->Arg(kIter);

}  // namespace benchmarks
}  // namespace principia
