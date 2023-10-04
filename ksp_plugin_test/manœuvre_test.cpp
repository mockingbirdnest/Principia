#include "ksp_plugin/manœuvre.hpp"

#include <cmath>
#include <memory>
#include <utility>

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/mock_ephemeris.hpp"  // 🧙 For MockEphemeris.
#include "physics/mock_rigid_reference_frame.hpp"  // 🧙 For MockRigidReferenceFrame.  // NOLINT
#include "physics/reference_frame.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/constants.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/ksp_plugin.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/make_not_null.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Return;
using ::testing::SetArgPointee;
using ::testing::StrictMock;
using ::testing::_;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_manœuvre;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::quantities::_uk;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_make_not_null;
using namespace principia::testing_utilities::_numerics;

class ManœuvreTest : public ::testing::Test {
 protected:
  using World = Barycentric;
  using Rendering = Frame<serialization::Frame::TestTag,
                          Arbitrary,
                          Handedness::Right,
                          serialization::Frame::TEST>;

  not_null<
      std::unique_ptr<StrictMock<MockRigidReferenceFrame<World, Rendering>>>>
  MakeMockReferenceFrame() {
    auto owned_mock_reference_frame = make_not_null_unique<
        StrictMock<MockRigidReferenceFrame<World, Rendering>>>();
    mock_reference_frame_ = owned_mock_reference_frame.get();
    return owned_mock_reference_frame;
  }

  Instant const t0_;
  StrictMock<MockRigidReferenceFrame<World, Rendering>> const*
      mock_reference_frame_;
  DiscreteTrajectory<World> discrete_trajectory_;
  DegreesOfFreedom<World> const dof_ = {
      World::origin + Displacement<World>({1 * Metre, 9 * Metre, 5 * Metre}),
      Velocity<World>(
          {8 * Metre / Second, 10 * Metre / Second, 4 * Metre / Second})};
  DegreesOfFreedom<Rendering> const rendering_dof_ = {
      Rendering::origin +
          Displacement<Rendering>({1 * Metre, 9 * Metre, 5 * Metre}),
      Velocity<Rendering>(
          {8 * Metre / Second, 10 * Metre / Second, 4 * Metre / Second})};
  Vector<Acceleration, World> const gravitational_acceleration_{
      {0 * Metre / Second / Second,
       0 * Metre / Second / Second,
       -StandardGravity}};
  RigidMotion<World, Rendering> const rigid_motion_{
      RigidTransformation<World, Rendering>(
          World::origin,
          Rendering::origin,
          OrthogonalMap<World, Rendering>::Identity()),
      World::nonrotating,
      World::unmoving};
  AcceleratedRigidMotion<World, Rendering> const accelerated_rigid_motion_{
      rigid_motion_,
      Bivector<AngularAcceleration, World>(),
      Vector<Acceleration, World>()};
};

TEST_F(ManœuvreTest, TimedBurn) {
  Vector<double, Frenet<Rendering>> e_y({0, 1, 0});

  Manœuvre<World, Rendering>::Intensity intensity;
  intensity.direction = e_y;
  intensity.duration = 1 * Second;
  Manœuvre<World, Rendering>::Timing timing;
  timing.initial_time = t0_;
  Manœuvre<World, Rendering>::Burn const burn{
      intensity,
      timing,
      /*thrust=*/1 * Newton,
      /*specific_impulse=*/1 * Newton * Second / Kilogram,
      MakeMockReferenceFrame(),
      /*is_inertially_fixed=*/true};
  Manœuvre<World, Rendering> manœuvre(/*initial_mass=*/2 * Kilogram, burn);
  EXPECT_EQ(1 * Newton, manœuvre.thrust());
  EXPECT_EQ(2 * Kilogram, manœuvre.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre.specific_impulse());
// Fixed in 192'027'508.
#if PRINCIPIA_COMPILER_MSVC && (_MSC_FULL_VER == 191'627'024 || \
                                _MSC_FULL_VER == 191'627'025 || \
                                _MSC_FULL_VER == 191'627'027)
  EXPECT_TRUE(e_y == manœuvre.direction());
#else
  EXPECT_EQ(e_y, manœuvre.direction());
#endif
  EXPECT_EQ(1 * Kilogram / Second, manœuvre.mass_flow());

  EXPECT_EQ(1 * Second, manœuvre.duration());
  EXPECT_EQ(std::log(2) * Metre / Second, manœuvre.Δv().Norm());
  EXPECT_EQ((2 - Sqrt(2)) * Second, manœuvre.time_to_half_Δv());
  EXPECT_EQ(1 * Kilogram, manœuvre.final_mass());

  EXPECT_EQ(t0_, manœuvre.initial_time());
  EXPECT_EQ(t0_ + 1 * Second, manœuvre.final_time());
  EXPECT_EQ(t0_ + (2 - Sqrt(2)) * Second, manœuvre.time_of_half_Δv());

  EXPECT_OK(discrete_trajectory_.Append(manœuvre.initial_time(), dof_));
  EXPECT_CALL(*mock_reference_frame_,
              ToThisFrameAtTime(manœuvre.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              MotionOfThisFrame(manœuvre.initial_time()))
      .WillOnce(Return(accelerated_rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              GravitationalAcceleration(manœuvre.initial_time(), _))
      .WillOnce(Return(gravitational_acceleration_));
  manœuvre.set_coasting_trajectory(discrete_trajectory_.segments().begin());
  auto const acceleration = manœuvre.InertialIntrinsicAcceleration();
  EXPECT_EQ(
      0 * Metre / Pow<2>(Second),
      acceleration(manœuvre.initial_time() - 1 * Second).Norm());
  EXPECT_THAT(acceleration(manœuvre.initial_time()).Norm(),
              AlmostEquals(0.5 * Metre / Pow<2>(Second), 1));
  EXPECT_THAT(acceleration(manœuvre.time_of_half_Δv()).Norm(),
              AlmostEquals(Sqrt(0.5) * Metre / Pow<2>(Second), 1));
  EXPECT_THAT(acceleration(manœuvre.final_time()).Norm(),
              AlmostEquals(1 * Metre / Pow<2>(Second), 1));
  EXPECT_EQ(0 * Metre / Pow<2>(Second),
            acceleration(manœuvre.final_time() + 1 * Second).Norm());
}

TEST_F(ManœuvreTest, TargetΔv) {
  Vector<double, Frenet<Rendering>> e_y({0, 1, 0});

  Manœuvre<World, Rendering>::Intensity intensity;
  intensity.Δv = e_y * Metre / Second;
  Manœuvre<World, Rendering>::Timing timing;
  timing.time_of_half_Δv = t0_;
  Manœuvre<World, Rendering>::Burn const burn{
      intensity,
      timing,
      /*thrust=*/1 * Newton,
      /*specific_impulse=*/1 * Newton * Second / Kilogram,
      MakeMockReferenceFrame(),
      /*is_inertially_fixed=*/true};
  Manœuvre<World, Rendering> manœuvre(/*initial_mass=*/2 * Kilogram, burn);
  EXPECT_EQ(1 * Newton, manœuvre.thrust());
  EXPECT_EQ(2 * Kilogram, manœuvre.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre.specific_impulse());
// Fixed in 192'027'508.
#if PRINCIPIA_COMPILER_MSVC && (_MSC_FULL_VER == 191'627'024 || \
                                _MSC_FULL_VER == 191'627'025 || \
                                _MSC_FULL_VER == 191'627'027)
  EXPECT_TRUE(e_y == manœuvre.direction());
#else
  EXPECT_EQ(e_y, manœuvre.direction());
#endif
  EXPECT_EQ(1 * Kilogram / Second, manœuvre.mass_flow());

  EXPECT_EQ(1 * Metre / Second, manœuvre.Δv().Norm());
  EXPECT_EQ((2 - 2 / e) * Second, manœuvre.duration());
  EXPECT_EQ((2 - 2 / Sqrt(e)) * Second, manœuvre.time_to_half_Δv());
  EXPECT_EQ((2 / e) * Kilogram, manœuvre.final_mass());

  EXPECT_EQ(t0_ - (2 - 2 / Sqrt(e)) * Second, manœuvre.initial_time());
  EXPECT_EQ(t0_ + (2 / Sqrt(e) - 2 / e) * Second, manœuvre.final_time());
  EXPECT_EQ(t0_, manœuvre.time_of_half_Δv());

  EXPECT_OK(discrete_trajectory_.Append(manœuvre.initial_time(), dof_));
  EXPECT_CALL(*mock_reference_frame_,
              ToThisFrameAtTime(manœuvre.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              MotionOfThisFrame(manœuvre.initial_time()))
      .WillOnce(Return(accelerated_rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              GravitationalAcceleration(manœuvre.initial_time(), _))
      .WillOnce(Return(gravitational_acceleration_));
  manœuvre.set_coasting_trajectory(discrete_trajectory_.segments().begin());
  auto const acceleration = manœuvre.InertialIntrinsicAcceleration();
  EXPECT_EQ(
      0 * Metre / Pow<2>(Second),
      acceleration(manœuvre.initial_time() - 1 * Second).Norm());
  EXPECT_THAT(acceleration(manœuvre.initial_time()).Norm(),
              AlmostEquals(0.5 * Metre / Pow<2>(Second), 1));
  EXPECT_EQ(Sqrt(e) / 2 * Metre / Pow<2>(Second),
            acceleration(manœuvre.time_of_half_Δv()).Norm());
  EXPECT_EQ((e / 2) * Metre / Pow<2>(Second),
            acceleration(manœuvre.final_time()).Norm());
  EXPECT_EQ(0 * Metre / Pow<2>(Second),
            acceleration(manœuvre.final_time() + 1 * Second).Norm());
}

TEST_F(ManœuvreTest, Apollo8SIVB) {
  // Data from [Sat69].  We use the reconstructed or actual values.

  // Table 2-2. Significant Event Times Summary.
  Instant const range_zero;
  Instant const s_ivb_1st_90_percent_thrust = range_zero +    530.53 * Second;
  Instant const s_ivb_1st_eco               = range_zero +    684.98 * Second;
  // Initiate S-IVB Restart Sequence and Start of Time Base 6 (T6).
  Instant const t6                          = range_zero +   9659.54 * Second;
  Instant const s_ivb_2nd_90_percent_thrust = range_zero + 10'240.02 * Second;
  Instant const s_ivb_2nd_eco               = range_zero + 10'555.51 * Second;

  // From Table 7-2. S-IVB Steady State Performance - First Burn.
  Force thrust_1st                  = 901'557 * Newton;
  Speed specific_impulse_1st        = 4'204.1 * Newton * Second / Kilogram;
  Variation<Mass> lox_flowrate_1st  = 178.16 * Kilogram / Second;
  Variation<Mass> fuel_flowrate_1st = 36.30 * Kilogram / Second;

  // From Table 7-7. S-IVB Steady State Performance - Second Burn.
  Force thrust_2nd                  = 897'548 * Newton;
  Speed specific_impulse_2nd        = 4199.2 * Newton * Second / Kilogram;
  Variation<Mass> lox_flowrate_2nd  = 177.70 * Kilogram / Second;
  Variation<Mass> fuel_flowrate_2nd = 36.01 * Kilogram / Second;

  // Table 21-5. Total Vehicle Mass, S-IVB First Burn Phase, Kilograms.
  Mass total_vehicle_at_s_ivb_1st_90_percent_thrust = 161143 * Kilogram;
  Mass total_vehicle_at_s_ivb_1st_eco               = 128095 * Kilogram;

  // Table 21-7. Total Vehicle Mass, S-IVB Second Burn Phase, Kilograms.
  Mass total_vehicle_at_s_ivb_2nd_90_percent_thrust = 126780 * Kilogram;
  Mass total_vehicle_at_s_ivb_2nd_eco               =  59285 * Kilogram;

  // An arbitrary direction, we're not testing this.
  Vector<double, Frenet<Rendering>> e_y({0, 1, 0});

  Manœuvre<World, Rendering>::Intensity first_burn_intensity;
  first_burn_intensity.direction = e_y;
  first_burn_intensity.duration = s_ivb_1st_eco - s_ivb_1st_90_percent_thrust;
  Manœuvre<World, Rendering>::Timing first_burn_timing;
  first_burn_timing.initial_time = s_ivb_1st_90_percent_thrust;
  Manœuvre<World, Rendering>::Burn const first_burn{
      first_burn_intensity,
      first_burn_timing,
      thrust_1st,
      specific_impulse_1st,
      MakeMockReferenceFrame(),
      /*is_inertially_fixed=*/true};
  Manœuvre<World, Rendering> first_manœuvre(
      total_vehicle_at_s_ivb_1st_90_percent_thrust,
      first_burn);

  EXPECT_THAT(RelativeError(lox_flowrate_1st + fuel_flowrate_1st,
                            first_manœuvre.mass_flow()),
              Lt(1e-4));
  EXPECT_THAT(RelativeError(total_vehicle_at_s_ivb_1st_eco,
                            first_manœuvre.final_mass()),
              Lt(1e-3));
  EXPECT_EQ(s_ivb_1st_eco, first_manœuvre.final_time());

  // Accelerations from Figure 4-4. Ascent Trajectory Acceleration Comparison.
  // Final acceleration from Table 4-2. Comparison of Significant Trajectory
  // Events.
  EXPECT_OK(discrete_trajectory_.Append(first_manœuvre.initial_time(), dof_));
  EXPECT_CALL(*mock_reference_frame_,
              ToThisFrameAtTime(first_manœuvre.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              MotionOfThisFrame(first_manœuvre.initial_time()))
      .WillOnce(Return(accelerated_rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              GravitationalAcceleration(first_manœuvre.initial_time(), _))
      .WillOnce(Return(gravitational_acceleration_));
  first_manœuvre.set_coasting_trajectory(
      discrete_trajectory_.segments().begin());
  auto const first_acceleration =
      first_manœuvre.InertialIntrinsicAcceleration();
  EXPECT_THAT(
      first_acceleration(first_manœuvre.initial_time()).Norm(),
      IsNear(5.6_(1) * Metre / Pow<2>(Second)));
  EXPECT_THAT(first_acceleration(range_zero + 600 * Second).Norm(),
              IsNear(6.16_(1) * Metre / Pow<2>(Second)));
  EXPECT_THAT(first_acceleration(first_manœuvre.final_time()).Norm(),
              IsNear(7.04_(1) * Metre / Pow<2>(Second)));

  Manœuvre<World, Rendering>::Intensity second_burn_intensity;
  second_burn_intensity.direction = e_y;
  second_burn_intensity.duration = s_ivb_2nd_eco - s_ivb_2nd_90_percent_thrust;
  Manœuvre<World, Rendering>::Timing second_burn_timing;
  second_burn_timing.initial_time = s_ivb_2nd_90_percent_thrust;
  Manœuvre<World, Rendering>::Burn const second_burn{
      second_burn_intensity,
      second_burn_timing,
      thrust_2nd,
      specific_impulse_2nd,
      MakeMockReferenceFrame(),
      /*is_inertially_fixed=*/true};
  Manœuvre<World, Rendering> second_manœuvre(
      total_vehicle_at_s_ivb_2nd_90_percent_thrust,
      second_burn);

  EXPECT_THAT(RelativeError(lox_flowrate_2nd + fuel_flowrate_2nd,
                            second_manœuvre.mass_flow()),
              Lt(2e-4));
  EXPECT_THAT(RelativeError(total_vehicle_at_s_ivb_2nd_eco,
                            second_manœuvre.final_mass()),
              Lt(2e-3));
  EXPECT_EQ(s_ivb_2nd_eco, second_manœuvre.final_time());

  // Accelerations from Figure 4-9. Injection Phase Acceleration Comparison.
  // Final acceleration from Table 4-2. Comparison of Significant Trajectory
  // Events.
  EXPECT_OK(discrete_trajectory_.Append(second_manœuvre.initial_time(), dof_));
  EXPECT_CALL(*mock_reference_frame_,
              ToThisFrameAtTime(second_manœuvre.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              MotionOfThisFrame(second_manœuvre.initial_time()))
      .WillOnce(Return(accelerated_rigid_motion_));
  EXPECT_CALL(*mock_reference_frame_,
              GravitationalAcceleration(second_manœuvre.initial_time(), _))
      .WillOnce(Return(gravitational_acceleration_));
  second_manœuvre.set_coasting_trajectory(
      discrete_trajectory_.segments().begin());
  auto const second_acceleration =
      second_manœuvre.InertialIntrinsicAcceleration();
  EXPECT_THAT(second_acceleration(second_manœuvre.initial_time()).Norm(),
              IsNear(7.08_(1) * Metre / Pow<2>(Second)));
  EXPECT_THAT(second_acceleration(t6 + 650 * Second).Norm(),
              IsNear(8.01_(1) * Metre / Pow<2>(Second)));
  EXPECT_THAT(second_acceleration(t6 + 700 * Second).Norm(),
              IsNear(8.9_(1) * Metre / Pow<2>(Second)));
  EXPECT_THAT(second_acceleration(t6 + 750 * Second).Norm(),
              IsNear(9.9_(1) * Metre / Pow<2>(Second)));
  EXPECT_THAT(second_acceleration(t6 + 850 * Second).Norm(),
              IsNear(12.97_(1) * Metre / Pow<2>(Second)));
  EXPECT_THAT(second_acceleration(second_manœuvre.final_time()).Norm(),
              IsNear(15.12_(1) * Metre / Pow<2>(Second)));

  EXPECT_THAT(second_manœuvre.Δv().Norm(),
              IsNear(3.2_(1) * Kilo(Metre) / Second));

  // From the Apollo 8 flight journal.
  EXPECT_THAT(AbsoluteError(10'519.6 * Foot / Second,
                            second_manœuvre.Δv().Norm()),
              Lt(20 * Metre / Second));
}

TEST_F(ManœuvreTest, Serialization) {
  auto mock_reference_frame = MakeMockReferenceFrame();
  auto const unowned_reference_frame = mock_reference_frame.get();
  Vector<double, Frenet<Rendering>> const e_y({0, 1, 0});

  Manœuvre<World, Rendering>::Intensity intensity;
  intensity.Δv = e_y * Metre / Second;
  Manœuvre<World, Rendering>::Timing timing;
  timing.time_of_half_Δv = t0_;
  Manœuvre<World, Rendering>::Burn const burn{
      intensity,
      timing,
      /*thrust=*/1 * Newton,
      /*specific_impulse=*/1 * Newton * Second / Kilogram,
      std::move(mock_reference_frame),
      /*is_inertially_fixed=*/true};
  Manœuvre<World, Rendering> manœuvre(/*initial_mass=*/2 * Kilogram, burn);

  serialization::ReferenceFrame serialized_mock_reference_frame;
  serialized_mock_reference_frame.MutableExtension(
      serialization::BodyCentredNonRotatingReferenceFrame::
          extension)->set_centre(666);
  EXPECT_CALL(*unowned_reference_frame, WriteToMessage(_))
      .WillOnce(SetArgPointee<0>(serialized_mock_reference_frame));
  serialization::Manoeuvre message;
  manœuvre.WriteToMessage(&message);

  EXPECT_TRUE(message.has_thrust());
  EXPECT_TRUE(message.has_initial_mass());
  EXPECT_TRUE(message.has_specific_impulse());
  EXPECT_TRUE(message.has_direction());
  EXPECT_TRUE(message.has_duration());
  EXPECT_TRUE(message.has_initial_time());
  EXPECT_TRUE(message.has_frame());

  MockEphemeris<World> ephemeris;
  MassiveBody const body(si::Unit<GravitationalParameter>);
  EXPECT_CALL(ephemeris, body_for_serialization_index(666))
      .WillOnce(Return(&body));
  EXPECT_CALL(ephemeris, trajectory(_))
      .WillOnce(Return(make_not_null<ContinuousTrajectory<World> const*>()));
  Manœuvre<World, Rendering> const manœuvre_read =
      Manœuvre<World, Rendering>::ReadFromMessage(message, &ephemeris);

  EXPECT_EQ(1 * Newton, manœuvre_read.thrust());
  EXPECT_EQ(2 * Kilogram, manœuvre_read.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre_read.specific_impulse());
// Fixed in 192'027'508.
#if PRINCIPIA_COMPILER_MSVC && (_MSC_FULL_VER == 191'627'024 || \
                                _MSC_FULL_VER == 191'627'025 || \
                                _MSC_FULL_VER == 191'627'027)
  EXPECT_TRUE(e_y == manœuvre.direction());
#else
  EXPECT_EQ(e_y, manœuvre.direction());
#endif
  EXPECT_EQ(1 * Kilogram / Second, manœuvre_read.mass_flow());

  EXPECT_EQ(1 * Metre / Second, manœuvre.Δv().Norm());
  EXPECT_EQ((2 - 2 / e) * Second, manœuvre_read.duration());
  EXPECT_EQ((2 - 2 / Sqrt(e)) * Second, manœuvre_read.time_to_half_Δv());
  EXPECT_EQ((2 / e) * Kilogram, manœuvre_read.final_mass());

  EXPECT_EQ(t0_ - (2 - 2 / Sqrt(e)) * Second, manœuvre_read.initial_time());
  EXPECT_EQ(t0_ + (2 / Sqrt(e) - 2 / e) * Second, manœuvre_read.final_time());
  EXPECT_EQ(t0_, manœuvre_read.time_of_half_Δv());
}

}  // namespace ksp_plugin
}  // namespace principia
