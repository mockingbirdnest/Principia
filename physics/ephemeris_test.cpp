#include "physics/ephemeris.hpp"

#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using integrators::BlanesMoan2002SRKN11B;
using si::Metre;
using si::Second;

namespace physics {

class EphemerisTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  EphemerisTest()
      : ephemeris_(std::vector<not_null<std::unique_ptr<MassiveBody>>>(),
                   std::vector<DegreesOfFreedom<World>>(),
                   t0_,
                   BlanesMoan2002SRKN11B<Position<World>>(),
                   1 * Second,
                   1 * Metre,
                   10 * Metre) {}

  Ephemeris<World> ephemeris_;
  Instant t0_;
};

TEST_F(EphemerisTest, Test) {
}

}  // namespace physics
}  // namespace principia
