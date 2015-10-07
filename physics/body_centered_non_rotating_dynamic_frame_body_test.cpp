#include "physics/body_centered_non_rotating_dynamic_frame.hpp"

#include <memory>

#include "astronomy/frames.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::Instant;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;

namespace physics {

class BodyCentredNonRotatingDynamicFrameTest : public ::testing::Test {
 protected:
  // The non-rotating frame centred on the big body.
  using Big = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, false /*inertial*/>;

  BodyCentredNonRotatingDynamicFrameTest()
      : period_(10 * ? * sqrt(5.0 / 7.0) * Second) {
    solar_system_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model_two_bodies_test.proto.txt",
        SOLUTION_DIR / "astronomy" / "initial_state_two_bodies_test.proto.txt");
    ephemeris_ = solar_system_.MakeEphemeris(
                     integrators::McLachlanAtela1992Order4Optimal<
                          Position<ICRFJ2000Equator>>(),
                     1 * Second,
                     1 * Metre);
    ephemeris_->Prolong(t0_ + 2 * period_);
    big_frame_ = std::make_unique<
                     BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Big>>(
                         ephemeris_.get(),
                         solar_system_.massive_body(*ephemeris_, "Big"));
  }

  Time const period_;
  Instant const t0_;
  std::unique_ptr<
      BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Big>> big_frame_;
  std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  SolarSystem<ICRFJ2000Equator> solar_system_;
};

TEST_F(BodyCentredNonRotatingDynamicFrameTest, ToBigFrame) {
  auto const to_big_frame = big_frame_->ToThisFrameAtTime(t0_ + 0 * period_);
}

}  // namespace physics
}  // namespace principia
