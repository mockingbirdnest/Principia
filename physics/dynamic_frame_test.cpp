
#include "physics/dynamic_frame.hpp"

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/mock_continuous_trajectory.hpp"
#include "physics/mock_ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace physics {
namespace internal_dynamic_frame {

using astronomy::ICRS;
using base::check_not_null;
using geometry::AngularVelocity;
using geometry::Barycentre;
using geometry::Displacement;
using geometry::Frame;
using geometry::Inertial;
using geometry::InfiniteFuture;
using geometry::InfinitePast;
using geometry::InnerProduct;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Velocity;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::McLachlanAtela1992Order4Optimal;
using quantities::GravitationalParameter;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using ::testing::A;
using ::testing::Eq;
using ::testing::InSequence;
using ::testing::Return;
using ::testing::StrictMock;
namespace si = quantities::si;

namespace {

using Circular = Frame<enum class CircularTag, Inertial>;
using Helical = Frame<enum class HelicalTag, Inertial>;

char constexpr big[] = "Big";
char constexpr small[] = "Small";

Vector<Acceleration, Circular> Gravity(Instant const& t,
                                       Position<Circular> const& q) {
  Displacement<Circular> const r = q - Circular::origin;
  auto const r² = r.Norm²();
  return -si::Unit<GravitationalParameter> * r / (Sqrt(r²) * r²);
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

  Instant t_min() const override {
    return InfinitePast;
  }

  Instant t_max() const override {
    return InfiniteFuture;
  }

  RigidMotion<OtherFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::DynamicFrame*> message) const override;

 private:
  Vector<Acceleration, OtherFrame> GravitationalAcceleration(
      Instant const& t,
      Position<OtherFrame> const& q) const override;
  AcceleratedRigidMotion<OtherFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const override;

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
             OtherFrame::nonrotating,
             origin_degrees_of_freedom_at_epoch_.velocity());
}

template<typename OtherFrame, typename ThisFrame>
void InertialFrame<OtherFrame, ThisFrame>::WriteToMessage(
    not_null<serialization::DynamicFrame*> message) const {}

template<typename OtherFrame, typename ThisFrame>
Vector<Acceleration, OtherFrame>
InertialFrame<OtherFrame, ThisFrame>::GravitationalAcceleration(
    Instant const& t,
    Position<OtherFrame> const& q) const {
  return gravity_(t, q);
}

template<typename OtherFrame, typename ThisFrame>
AcceleratedRigidMotion<OtherFrame, ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::MotionOfThisFrame(
    Instant const& t) const {
  return AcceleratedRigidMotion<OtherFrame, ThisFrame>(
      ToThisFrameAtTime(t),
      /*angular_acceleration_of_to_frame=*/{},
      /*acceleration_of_to_frame_origin=*/{});
}

}  // namespace

class DynamicFrameTest : public testing::Test {
 protected:
  using MockFrame = Frame<serialization::Frame::TestTag,
                          Arbitrary,
                          Handedness::Right,
                          serialization::Frame::TEST1>;

  DynamicFrameTest()
      : solar_system_(SOLUTION_DIR / "astronomy" /
                          "test_gravity_model_two_bodies.proto.txt",
                      SOLUTION_DIR / "astronomy" /
                          "test_initial_state_two_bodies_circular.proto.txt"),
        t0_(solar_system_.epoch()),
        ephemeris_(solar_system_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymplecticRungeKuttaNyströmIntegrator<
                    McLachlanAtela1992Order4Optimal,
                    Position<ICRS>>(),
                /*step=*/10 * Milli(Second)))),
        big_(solar_system_.massive_body(*ephemeris_, big)),
        big_initial_state_(solar_system_.degrees_of_freedom(big)),
        big_gravitational_parameter_(
            solar_system_.gravitational_parameter(big)),
        small_(solar_system_.massive_body(*ephemeris_, small)),
        small_initial_state_(solar_system_.degrees_of_freedom(small)),
        small_gravitational_parameter_(
            solar_system_.gravitational_parameter(small)) {
    EXPECT_CALL(mock_ephemeris_,
                trajectory(solar_system_.massive_body(*ephemeris_, big)))
        .WillOnce(Return(&mock_big_trajectory_));
    EXPECT_CALL(mock_ephemeris_,
                trajectory(solar_system_.massive_body(*ephemeris_, small)))
        .WillOnce(Return(&mock_small_trajectory_));
    mock_frame_ =
        std::make_unique<BarycentricRotatingDynamicFrame<ICRS, MockFrame>>(
            &mock_ephemeris_, big_, small_);
  }

  SolarSystem<ICRS> solar_system_;
  Instant const t0_;
  std::unique_ptr<Ephemeris<ICRS>> const ephemeris_;
  MassiveBody const* const big_;
  DegreesOfFreedom<ICRS> const big_initial_state_;
  GravitationalParameter const big_gravitational_parameter_;
  MassiveBody const* const small_;
  DegreesOfFreedom<ICRS> const small_initial_state_;
  GravitationalParameter const small_gravitational_parameter_;
  StrictMock<MockEphemeris<ICRS>> mock_ephemeris_;
  StrictMock<MockContinuousTrajectory<ICRS>> mock_big_trajectory_;
  StrictMock<MockContinuousTrajectory<ICRS>> mock_small_trajectory_;
  std::unique_ptr<BarycentricRotatingDynamicFrame<ICRS, MockFrame>> mock_frame_;

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
          /*epoch=*/Instant() ,
          OrthogonalMap<Circular, Helical>::Identity(),
          &Gravity);
};

// Two bodies in rotation with their barycentre at rest.  The test point is at
// the origin and in motion.  The acceleration is purely due to Coriolis.
TEST_F(DynamicFrameTest, CoriolisAcceleration) {
  Instant const t = t0_ + 0 * Second;
  // The velocity is opposed to the motion and away from the centre.
  DegreesOfFreedom<MockFrame> const point_dof =
      {Displacement<MockFrame>({0 * Metre, 0 * Metre, 0 * Metre}) +
           MockFrame::origin,
       Velocity<MockFrame>({(80 - 30) * Metre / Second,
                            (-60 - 40) * Metre / Second,
                            0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const big_dof =
      {Displacement<ICRS>({0.8 * Metre, -0.6 * Metre, 0 * Metre}) +
       ICRS::origin,
       Velocity<ICRS>(
           {-16 * Metre / Second, 12 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const small_dof =
      {Displacement<ICRS>({5 * Metre, 5 * Metre, 0 * Metre}) + ICRS::origin,
       Velocity<ICRS>(
           {40 * Metre / Second, -30 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const barycentre_dof =
      Barycentre<DegreesOfFreedom<ICRS>, GravitationalParameter>(
          {big_dof, small_dof},
          {big_gravitational_parameter_, small_gravitational_parameter_});
  EXPECT_THAT(barycentre_dof.position() - ICRS::origin,
              Eq(Displacement<ICRS>({2 * Metre, 1 * Metre, 0 * Metre})));
  EXPECT_THAT(barycentre_dof.velocity(), Eq(ICRS::unmoving));

  EXPECT_CALL(mock_big_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(big_dof));
  EXPECT_CALL(mock_small_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(small_dof));
  {
    InSequence s;
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(big_), t))
        .WillOnce(Return(Vector<Acceleration, ICRS>({
                             120 * Metre / Pow<2>(Second),
                             160 * Metre / Pow<2>(Second),
                             0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(small_), t))
        .WillOnce(Return(Vector<Acceleration, ICRS>({
                             -300 * Metre / Pow<2>(Second),
                             -400 * Metre / Pow<2>(Second),
                             0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMasslessBody(
                    A<Position<ICRS> const&>(), t))
        .WillOnce(Return(Vector<Acceleration, ICRS>()));
  }

  // The Coriolis acceleration is towards the centre and opposed to the motion.
  EXPECT_THAT(mock_frame_->GeometricAcceleration(t, point_dof),
              AlmostEquals(Vector<Acceleration, MockFrame>({
                               (-1200 - 800) * Metre / Pow<2>(Second),
                               (-1600 + 600) * Metre / Pow<2>(Second),
                               0 * Metre / Pow<2>(Second)}), 0));
}

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

}  // namespace internal_dynamic_frame
}  // namespace physics
}  // namespace principia
