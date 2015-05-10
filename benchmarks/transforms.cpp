
#include <memory>

#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
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

struct TrajectoryHolder {
  Trajectory<World1> const& primary();
  Trajectory<World1> const& secondary();

  not_null<std::unique_ptr<Trajectory<World1>>> primary;
  not_null<std::unique_ptr<Trajectory<World1>>> secondary;
};

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

void BM_BodyCentredNonRotating(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  TrajectoryHolder holder;

  Time const Δt = 1 * Hour;
  int const steps = 100000;


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
  holder.primary = NewCircularTrajectory(&earth,
                                         center,
                                         earth_initial_position,
                                         earth_angular_velocity,
                                         Δt,
                                         steps);

  MasslessBody probe;
  Position<World1> probe_initial_position =
      World1::origin + Displacement<World1>({0.5 * si::AstronomicalUnit,
                                             -1 * si::AstronomicalUnit,
                                             0 * si::AstronomicalUnit});
  Velocity<World1> probe_velocity =
      Velocity<World1>({0 * SIUnit<Speed>(),
                        100 * Kilo(Metre) / Second,
                        0 * SIUnit<Speed>()});
  auto const probe_trajectory = NewLinearTrajectory(&probe,
                                                    probe_initial_position,
                                                    probe_velocity,
                                                    Δt,
                                                    steps);

  auto transforms = Transforms<TrajectoryHolder, World1, World2, World1>::BodyCentredNonRotating(holder, );

  while (state.KeepRunning()) {
  }
}
BENCHMARK(BM_BodyCentredNonRotating);

}  // namespace benchmarks
}  // namespace principia
