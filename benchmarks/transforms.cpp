
// .\Release\benchmarks.exe --benchmark_filter=Rotating --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/05/10-13:47:31
// Benchmark                                       Time(ns)    CPU(ns) Iterations  // NOLINT(whitespace/line_length)
// ------------------------------------------------------------------------------  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<false>/1000k         642744994  639604100          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<false>/1000k        1387956357 1388408900          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<false>/1000k        1382258868 1388408900          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<false>/1000k        1383145785 1372808800          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<false>/1000k        1396267152 1388408900          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<false>/1000k_mean   1238474631 1235527920          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<false>/1000k_stddev  297906266  298023161          0  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<true>/1000k         1058213210 1045206700          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<true>/1000k         1207992053 1201207700          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<true>/1000k         1192008948 1201207700          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<true>/1000k         1200190759 1201207700          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<true>/1000k         1192056871 1201207700          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<true>/1000k_mean    1170092368 1170007500          1  // NOLINT(whitespace/line_length)
// BM_BodyCentredNonRotating<true>/1000k_stddev    56253037   62400400          0  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<false>/1000k           1890934682 1887612100          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<false>/1000k           3090300775 3073219700          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<false>/1000k           3099311090 3073219700          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<false>/1000k           3096378779 3104419900          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<false>/1000k           3086464858 3057619600          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<false>/1000k_mean      2852678037 2839218200          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<false>/1000k_stddev     480892720  476045940          0  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<true>/1000k            2319008565 2308814800          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<true>/1000k            2002243018 1996812800          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<true>/1000k            1999951100 2012412900          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<true>/1000k            1993252492 1996812800          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<true>/1000k            1980318046 1965612600          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<true>/1000k_mean       2058954644 2056093180          1  // NOLINT(whitespace/line_length)
// BM_BarycentricRotating<true>/1000k_stddev      130250647  127272343          0  // NOLINT(whitespace/line_length)

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
#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/trajectory.hpp"
#include "physics/transforms.hpp"
#include "serialization/geometry.pb.h"

// This must come last because apparently it redefines CDECL.
#include "benchmark/benchmark.h"

namespace principia {

using astronomy::EarthMass;
using astronomy::JulianYear;
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
using physics::Body;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using physics::MasslessBody;
using physics::Trajectory;
using physics::Transforms;
using si::AstronomicalUnit;
using si::Hour;
using si::Kilo;
using si::Metre;
using si::Radian;
using si::Second;

namespace benchmarks {

using World1 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST1, true>;
using World2 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST2, false>;

class TrajectoryHolder {
 public:
  explicit TrajectoryHolder(
      not_null<std::unique_ptr<Trajectory<World1>>> trajectory);

  Trajectory<World1> const& trajectory() const;

 private:
  std::unique_ptr<Trajectory<World1>> trajectory_;
};

TrajectoryHolder::TrajectoryHolder(
    not_null<std::unique_ptr<Trajectory<World1>>> trajectory)
    // TODO(phl): Y U NO MOV!
    : trajectory_(trajectory.release()) {}

Trajectory<World1> const& TrajectoryHolder::trajectory() const {
  return *trajectory_;
}

not_null<std::unique_ptr<Trajectory<World1>>> NewCircularTrajectory(
    not_null<Body const*> const body,
    Position<World1> const& center,
    Position<World1> const& initial,
    AngularVelocity<World1> const& angular_velocity,
    Time const& Δt,
    int const steps) {
  not_null<std::unique_ptr<Trajectory<World1>>> trajectory =
      std::make_unique<Trajectory<World1>>(body);
  Displacement<World1> const radius = initial - center;
  for (int i = 0; i < steps; ++i) {
    Time const t_i = i * Δt;
    Displacement<World1> const displacement_i =
        Exp(angular_velocity * t_i)(radius);
    Velocity<World1> const velocity_i =
        angular_velocity * displacement_i / Radian;
    trajectory->Append(Instant(t_i),
                       DegreesOfFreedom<World1>(initial + displacement_i,
                                                velocity_i));
  }
  return trajectory;
}

not_null<std::unique_ptr<Trajectory<World1>>> NewLinearTrajectory(
    not_null<Body const*> const body,
    Position<World1> const& initial,
    Velocity<World1> const& velocity,
    Time const& Δt,
    int const steps) {
  not_null<std::unique_ptr<Trajectory<World1>>> trajectory =
      std::make_unique<Trajectory<World1>>(body);
  for (int i = 0; i < steps; ++i) {
    Time const t_i = i * Δt;
    Displacement<World1> const displacement_i = velocity * t_i;
    trajectory->Append(Instant(t_i),
                       DegreesOfFreedom<World1>(initial + displacement_i,
                                                velocity));
  }
  return trajectory;
}

// This code is derived from Plugin::RenderTrajectory.
std::vector<std::pair<Position<World1>,
                      Position<World1>>> ApplyTransform(
    not_null<Body const*> const body,
    not_null<Transforms<TrajectoryHolder, World1, World2, World1>*> const
        transforms,
    Trajectory<World1>::TransformingIterator<World2> const& actual_it) {
  std::vector<std::pair<Position<World1>,
                        Position<World1>>> result;

  // First build the trajectory resulting from the first transform.
  Trajectory<World2> intermediate_trajectory(body);
  for (auto it = actual_it; !it.at_end(); ++it) {
    intermediate_trajectory.Append(it.time(), it.degrees_of_freedom());
  }

  // Then build the final result using the second transform.
  auto initial_it = transforms->second(intermediate_trajectory);
  if (!initial_it.at_end()) {
    for (auto final_it = initial_it;
         ++final_it, !final_it.at_end();
         initial_it = final_it) {
      result.emplace_back(initial_it.degrees_of_freedom().position(),
                          final_it.degrees_of_freedom().position());
    }
  }
  return result;
}

template<bool cache>
void BM_BodyCentredNonRotating(
    benchmark::State& state) {  // NOLINT(runtime/references)
  state.PauseTiming();

  Time const Δt = 1 * Hour;
  int const steps = state.range_x();

  MassiveBody earth(astronomy::EarthMass);
  Position<World1> center = World1::origin;
  Position<World1> earth_initial_position =
      World1::origin + Displacement<World1>({1 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  AngularVelocity<World1> earth_angular_velocity =
      AngularVelocity<World1>({0 * SIUnit<AngularFrequency>(),
                               0 * SIUnit<AngularFrequency>(),
                               2 * π * Radian / JulianYear});
  TrajectoryHolder earth_holder(NewCircularTrajectory(&earth,
                                                      center,
                                                      earth_initial_position,
                                                      earth_angular_velocity,
                                                      Δt,
                                                      steps));

  MasslessBody probe;
  Position<World1> probe_initial_position =
      World1::origin + Displacement<World1>({0.5 * si::AstronomicalUnit,
                                             -1 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  Velocity<World1> probe_velocity =
      Velocity<World1>({0 * SIUnit<Speed>(),
                        100 * Kilo(Metre) / Second,
                        0 * SIUnit<Speed>()});
  TrajectoryHolder probe_holder(NewLinearTrajectory(&probe,
                                                    probe_initial_position,
                                                    probe_velocity,
                                                    Δt,
                                                    steps));

  auto transforms = Transforms<TrajectoryHolder, World1, World2, World1>::
      BodyCentredNonRotating(earth_holder,
                             &TrajectoryHolder::trajectory);
  if (cache) {
    transforms->set_cacheable(&TrajectoryHolder::trajectory);
  }

  state.ResumeTiming();
  while (state.KeepRunning()) {
    auto v = ApplyTransform(&probe,
                   transforms.get(),
                   transforms->first(probe_holder,
                                     &TrajectoryHolder::trajectory));
  }
}

template<bool cache>
void BM_BarycentricRotating(
    benchmark::State& state) {  // NOLINT(runtime/references)
  state.PauseTiming();

  Time const Δt = 1 * Hour;
  int const steps = state.range_x();

  MassiveBody earth(astronomy::EarthMass);
  Position<World1> earth_center = World1::origin;
  Position<World1> earth_initial_position =
      World1::origin + Displacement<World1>({1 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  AngularVelocity<World1> earth_angular_velocity =
      AngularVelocity<World1>({0 * SIUnit<AngularFrequency>(),
                               0 * SIUnit<AngularFrequency>(),
                               2 * π * Radian / JulianYear});
  TrajectoryHolder earth_holder(NewCircularTrajectory(&earth,
                                                      earth_center,
                                                      earth_initial_position,
                                                      earth_angular_velocity,
                                                      Δt,
                                                      steps));

  MassiveBody thera(astronomy::EarthMass);
  Position<World1> thera_center =
      World1::origin + Displacement<World1>({2 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  Position<World1> thera_initial_position =
      World1::origin + Displacement<World1>({-0.5 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  AngularVelocity<World1> thera_angular_velocity =
      AngularVelocity<World1>({0 * SIUnit<AngularFrequency>(),
                               0 * SIUnit<AngularFrequency>(),
                               6 * Radian / JulianYear});
  TrajectoryHolder thera_holder(NewCircularTrajectory(&thera,
                                                      thera_center,
                                                      thera_initial_position,
                                                      thera_angular_velocity,
                                                      Δt,
                                                      steps));

  MasslessBody probe;
  Position<World1> probe_initial_position =
      World1::origin + Displacement<World1>({0.5 * si::AstronomicalUnit,
                                             -1 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  Velocity<World1> probe_velocity =
      Velocity<World1>({0 * SIUnit<Speed>(),
                        100 * Kilo(Metre) / Second,
                        0 * SIUnit<Speed>()});
  TrajectoryHolder probe_holder(NewLinearTrajectory(&probe,
                                                    probe_initial_position,
                                                    probe_velocity,
                                                    Δt,
                                                    steps));

  auto transforms = Transforms<TrajectoryHolder, World1, World2, World1>::
      BarycentricRotating(earth_holder,
                          thera_holder,
                          &TrajectoryHolder::trajectory);
  if (cache) {
    transforms->set_cacheable(&TrajectoryHolder::trajectory);
  }

  state.ResumeTiming();
  while (state.KeepRunning()) {
    auto v = ApplyTransform(&probe,
                   transforms.get(),
                   transforms->first(probe_holder,
                                     &TrajectoryHolder::trajectory));
  }
}

int const kIter = 1000 << 10;

BENCHMARK_TEMPLATE(BM_BodyCentredNonRotating, false)->Arg(kIter);
BENCHMARK_TEMPLATE(BM_BodyCentredNonRotating, true)->Arg(kIter);
BENCHMARK_TEMPLATE(BM_BarycentricRotating, false)->Arg(kIter);
BENCHMARK_TEMPLATE(BM_BarycentricRotating, true)->Arg(kIter);

}  // namespace benchmarks
}  // namespace principia
