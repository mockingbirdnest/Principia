#include "physics/dynamic_frame.hpp"

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"

namespace principia {

using geometry::Frame;
using geometry::InnerProduct;
using quantities::GravitationalParameter;
using quantities::Sqrt;

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
          Velocity<Circular>(
              {0 * Metre / Second, 0 * Metre / Second, 1 * Metre / Second}),
          Circular::origin /*origin_at_epoch*/,
          Instant() /*epoch*/,
          OrthogonalMap<Circular, Helical>::Identity(),
          &Gravity);
};

TEST_F(DynamicFrameTest, Helix) {
  auto helix_frenet_frame = helix_frame_.FrenetFrame(
      Instant(),
      helix_frame_.ToThisFrameAtTime(Instant())(circular_degrees_of_freedom_));
  Vector<double, Circular> tangent =
      helix_frenet_frame(Vector<double, Circular>({1, 0, 0}));
  Vector<double, Circular> normal =
      helix_frenet_frame(Vector<double, Circular>({0, 1, 0}));
  LOG(ERROR)<<tangent;
  LOG(ERROR)<<normal;
}

}  // namespace physics
}  // namespace principia
