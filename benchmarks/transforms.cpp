
#include <memory>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/trajectory.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::Exp;
using geometry::Frame;
using geometry::Position;
using geometry::Velocity;
using quantities::Time;
using physics::Body;
using physics::DegreesOfFreedom;
using physics::Trajectory;
using si::Radian;

namespace benchmarks {

using World1 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST1, true>;
using World2 = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST2, false>;

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

}  // namespace benchmarks
}  // namespace principia
