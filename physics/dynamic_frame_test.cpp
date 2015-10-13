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

// An inertial frame.
template<typename OtherFrame, typename ThisFrame>
class InertialFrame : public DynamicFrame<OtherFrame, ThisFrame> {
 public:
  InertialFrame(
      DegreesOfFreedom<OtherFrame> const& origin_degrees_of_freedom_at_epoch,
      Instant const& epoch,
      OrthogonalMap<OtherFrame, ThisFrame> const& orthogonal_map,
      std::function<Vector<Acceleration, OtherFrame>(
          Instant const& t,
          Position<OtherFrame> const& q)> gravity);

  RigidMotion<OtherFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;
  RigidMotion<ThisFrame, OtherFrame> FromThisFrameAtTime(
      Instant const& t) const override;

  // The acceleration due to gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;

 private:
  DegreesOfFreedom<OtherFrame> const origin_degrees_of_freedom_at_epoch_;
  Instant const epoch_;
  OrthogonalMap<OtherFrame, ThisFrame> const orthogonal_map_;
  std::function<Vector<Acceleration, OtherFrame>(
      Instant const& t,
      Position<OtherFrame> const& q)> gravity_;
};

template<typename OtherFrame, typename ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::InertialFrame(
    DegreesOfFreedom<OtherFrame> const& origin_degrees_of_freedom_at_epoch,
    Instant const& epoch,
    OrthogonalMap<OtherFrame, ThisFrame> const& orthogonal_map,
    std::function<Vector<Acceleration, OtherFrame>(
        Instant const& t,
        Position<OtherFrame> const& q)> gravity)
    : origin_degrees_of_freedom_at_epoch_(origin_degrees_of_freedom_at_epoch),
      epoch_(epoch),
      orthogonal_map_(orthogonal_map),
      gravity_(std::move(gravity)) {}

template<typename OtherFrame, typename ThisFrame>
RigidMotion<OtherFrame, ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  return RigidMotion<OtherFrame, ThisFrame>(
             RigidTransformation<OtherFrame, ThisFrame>(
                 origin_degrees_of_freedom_at_epoch_.position() +
                     (t - epoch_) *
                         origin_degrees_of_freedom_at_epoch_.velocity(),
                 ThisFrame::origin,
                 orthogonal_map_),
             AngularVelocity<ThisFrame>(),
             -orthogonal_map_(origin_degrees_of_freedom_at_epoch_.velocity()));
}

template<typename OtherFrame, typename ThisFrame>
RigidMotion<ThisFrame, OtherFrame>
InertialFrame<OtherFrame, ThisFrame>::FromThisFrameAtTime(
    Instant const& t) const {
  return ToThisFrameAtTime(t).Inverse();
}

template<typename OtherFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  return orthogonal_map_(
      gravity_(t, FromThisFrameAtTime(t)(degrees_of_freedom).position()));
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
