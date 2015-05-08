#include <memory>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using base::not_null;
using geometry::AngularVelocity;
using geometry::Frame;
using geometry::Position;
using quantities::Time;
using physics::Body;
using physics::Trajectory;

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
    Time const& delta_t) {
  not_null<std::unique_ptr<Trajectory<World1>>> trajectory =
      std::make_unique<Trajectory<World1>>(body);
  Displacement<World1> const radius = initial - center;
  auto delta_theta = delta_t * angular_velocity;
}

}  // namespace benchmarks
}  // namespace principia
