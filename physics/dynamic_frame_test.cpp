#include "physics/dynamic_frame.hpp"

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {

using geometry::Frame;
using geometry::InnerProduct;
using quantities::GravitationalParameter;
using quantities::Sqrt;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;

namespace physics {

namespace {

using Circular =
    Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;
using Helical =
    Frame<serialization::Frame::TestTag, serialization::Frame::TEST1, true>;

Vector<Acceleration, Circular> Gravity(Instant const& t,
                                       Position<Circular> const& q) {
  Displacement<Circular> const r = q - Circular::origin;
  auto const r2 = InnerProduct(r, r);
  return -SIUnit<GravitationalParameter>() * r / (Sqrt(r2) * r2);
}

}  // namespace

class DynamicFrameTest : public testing::Test {
 protected:
  DegreesOfFreedom<Circular> circular_degrees_of_freedom_ = {
      Circular::origin +
          Displacement<Circular>({1 * Metre, 0 * Metre, 0 * Metre}),
      Velocity<Circular>(
          {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second})};

  InertialFrame<Circular, Helical> helix_frame_ =
      InertialFrame<Circular, Helical>(
          {Circular::origin,
           Velocity<Circular>(
               {0 * Metre / Second, 0 * Metre / Second, 1 * Metre / Second})},
          Instant() /*epoch*/,
          OrthogonalMap<Circular, Helical>::Identity(),
          &Gravity);
};

TEST_F(DynamicFrameTest, Helix) {
  auto helix_frenet_frame = helix_frame_.FrenetFrame(
      Instant(),
      helix_frame_.ToThisFrameAtTime(Instant())(circular_degrees_of_freedom_));
  Vector<double, Helical> tangent =
      helix_frenet_frame(Vector<double, Frenet<Helical>>({1, 0, 0}));
  Vector<double, Helical> normal =
      helix_frenet_frame(Vector<double, Frenet<Helical>>({0, 1, 0}));
  EXPECT_THAT(normal, Componentwise(-1, 0, 0));
  EXPECT_THAT(tangent,
              Componentwise(0,
                            AlmostEquals(Sqrt(0.5), 1),
                            AlmostEquals(-Sqrt(0.5), 1)));
}

}  // namespace physics
}  // namespace principia
