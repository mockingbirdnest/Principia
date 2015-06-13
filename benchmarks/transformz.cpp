
// .\release\benchmarks.exe --benchmark_filter=Transformz --benchmark_repetitions=5          // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/06/13-16:59:51
// Benchmark                                                 Time(ns)    CPU(ns) Iterations  // NOLINT(whitespace/line_length)
// ----------------------------------------------------------------------------------------  // NOLINT(whitespace/line_length)
// BM_TransformzBodyCentredNonRotating<false>/1000k         786032176  780005000          1  // NOLINT(whitespace/line_length)
// BM_TransformzBodyCentredNonRotating<false>/1000k         786073899  795605100          1  // NOLINT(whitespace/line_length)
// BM_TransformzBodyCentredNonRotating<false>/1000k         786113954  780005000          1  // NOLINT(whitespace/line_length)
// BM_TransformzBodyCentredNonRotating<false>/1000k         781030154  780005000          1  // NOLINT(whitespace/line_length)
// BM_TransformzBodyCentredNonRotating<false>/1000k         787456489  780005000          1  // NOLINT(whitespace/line_length)
// BM_TransformzBodyCentredNonRotating<false>/1000k_mean    785341334  783125020          1  // NOLINT(whitespace/line_length)
// BM_TransformzBodyCentredNonRotating<false>/1000k_stddev    2221306    6240040          0  // NOLINT(whitespace/line_length)
// BM_TransformzBarycentricRotating<false>/1000k           2289228892 2293214700          1  // NOLINT(whitespace/line_length)
// BM_TransformzBarycentricRotating<false>/1000k           2297225213 2277614600          1  // NOLINT(whitespace/line_length)
// BM_TransformzBarycentricRotating<false>/1000k           2280222869 2277614600          1  // NOLINT(whitespace/line_length)
// BM_TransformzBarycentricRotating<false>/1000k           2289222217 2277614600          1  // NOLINT(whitespace/line_length)
// BM_TransformzBarycentricRotating<false>/1000k           2278222776 2277614600          1  // NOLINT(whitespace/line_length)
// BM_TransformzBarycentricRotating<false>/1000k_mean      2286824393 2280734620          1  // NOLINT(whitespace/line_length)
// BM_TransformzBarycentricRotating<false>/1000k_stddev       6888776    6240040          0  // NOLINT(whitespace/line_length)

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
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/trajectory.hpp"
#include "physics/transformz.hpp"
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
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using physics::MasslessBody;
using physics::Trajectory;
using physics::Transformz;
using si::AstronomicalUnit;
using si::Hour;
using si::Kilo;
using si::Metre;
using si::Radian;
using si::Second;

namespace benchmarks {

namespace {
const Length kLowTolerance = 0.001 * Metre;
const Length kHighTolerance = 0.01 * Metre;
}  // namespace

using World1 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST1, true>;
using World2 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST2, false>;

class TrajectoryHolderz {
 public:
  explicit TrajectoryHolderz(not_null<Trajectory<World1>*> trajectory);

  Trajectory<World1> const& trajectory() const;

 private:
  not_null<Trajectory<World1>*> trajectory_;
};

TrajectoryHolderz::TrajectoryHolderz(not_null<Trajectory<World1>*> trajectory)
    : trajectory_(trajectory) {}

Trajectory<World1> const& TrajectoryHolderz::trajectory() const {
  return *trajectory_;
}

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
                      Position<World1>>> ApplyTransform(
    not_null<Body const*> const body,
    not_null<Transformz<TrajectoryHolderz, World1, World2, World1>*> const
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
  std::unique_ptr<Position<World1>> last_position;  // std::optional.
  for (auto it = transforms->second(intermediate_trajectory.last().time(),
                                    intermediate_trajectory);
       !it.at_end();
       ++it) {
    Position<World1> const& position = it.degrees_of_freedom().position();
    if (last_position == nullptr) {
      last_position = std::make_unique<Position<World1>>(position);
    } else {
      result.emplace_back(*last_position, position);
      *last_position = position;
    }
  }
  return result;
}

template<bool cache>
void BM_TransformzBodyCentredNonRotating(
    benchmark::State& state) {  // NOLINT(runtime/references)

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
  ContinuousTrajectory<World1> earth_trajectory(
      Δt, kLowTolerance, kHighTolerance);
  FillCircularTrajectory<World1, ContinuousTrajectory>(center,
                                                       earth_initial_position,
                                                       earth_angular_velocity,
                                                       Δt,
                                                       steps,
                                                       &earth_trajectory);

  MasslessBody probe;
  Position<World1> probe_initial_position =
      World1::origin + Displacement<World1>({0.5 * si::AstronomicalUnit,
                                             -1 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  Velocity<World1> probe_velocity =
      Velocity<World1>({0 * SIUnit<Speed>(),
                        100 * Kilo(Metre) / Second,
                        0 * SIUnit<Speed>()});
  Trajectory<World1> probe_trajectory(&probe);
  TrajectoryHolderz probe_holder(&probe_trajectory);
  FillLinearTrajectory<World1, Trajectory>(probe_initial_position,
                                           probe_velocity,
                                           Δt,
                                           steps,
                                           &probe_trajectory);

  auto transforms = Transformz<TrajectoryHolderz, World1, World2, World1>::
      BodyCentredNonRotating(earth, earth_trajectory, earth_trajectory);
  if (cache) {
    transforms->set_cacheable(&TrajectoryHolderz::trajectory);
  }

  while (state.KeepRunning()) {
    auto v = ApplyTransform(&probe,
                   transforms.get(),
                   transforms->first(probe_holder,
                                     &TrajectoryHolderz::trajectory));
  }
}

template<bool cache>
void BM_TransformzBarycentricRotating(
    benchmark::State& state) {  // NOLINT(runtime/references)

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
  ContinuousTrajectory<World1> earth_trajectory(
      Δt, kLowTolerance, kHighTolerance);
  FillCircularTrajectory<World1, ContinuousTrajectory>(earth_center,
                                                       earth_initial_position,
                                                       earth_angular_velocity,
                                                       Δt,
                                                       steps,
                                                       &earth_trajectory);

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
  ContinuousTrajectory<World1> thera_trajectory(
      Δt, kLowTolerance, kHighTolerance);
  FillCircularTrajectory<World1, ContinuousTrajectory>(thera_center,
                                                       thera_initial_position,
                                                       thera_angular_velocity,
                                                       Δt,
                                                       steps,
                                                       &thera_trajectory);

  MasslessBody probe;
  Position<World1> probe_initial_position =
      World1::origin + Displacement<World1>({0.5 * si::AstronomicalUnit,
                                             -1 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  Velocity<World1> probe_velocity =
      Velocity<World1>({0 * SIUnit<Speed>(),
                        100 * Kilo(Metre) / Second,
                        0 * SIUnit<Speed>()});
  Trajectory<World1> probe_trajectory(&probe);
  TrajectoryHolderz probe_holder(&probe_trajectory);
  FillLinearTrajectory<World1, Trajectory>(probe_initial_position,
                                           probe_velocity,
                                           Δt,
                                           steps,
                                           &probe_trajectory);

  auto transforms = Transformz<TrajectoryHolderz, World1, World2, World1>::
      BarycentricRotating(earth,
                          earth_trajectory,
                          earth_trajectory,
                          thera,
                          thera_trajectory,
                          thera_trajectory);
  if (cache) {
    transforms->set_cacheable(&TrajectoryHolderz::trajectory);
  }

  while (state.KeepRunning()) {
    auto v = ApplyTransform(&probe,
                   transforms.get(),
                   transforms->first(probe_holder,
                                     &TrajectoryHolderz::trajectory));
  }
}

int const kIter = (1000 << 10) + 1;

BENCHMARK_TEMPLATE(BM_TransformzBodyCentredNonRotating, false)->Arg(kIter);
// BENCHMARK_TEMPLATE(BM_TransformzBodyCentredNonRotating, true)->Arg(kIter);
BENCHMARK_TEMPLATE(BM_TransformzBarycentricRotating, false)->Arg(kIter);
// BENCHMARK_TEMPLATE(BM_TransformzBarycentricRotating, true)->Arg(kIter);

}  // namespace benchmarks
}  // namespace principia
