
#include "physics/body.hpp"

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace physics {
namespace internal_body {

using astronomy::ICRS;
using astronomy::J2000;
using geometry::AngleBetween;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Displacement;
using geometry::Frame;
using geometry::Instant;
using geometry::Normalize;
using geometry::OrientedAngleBetween;
using geometry::Position;
using geometry::RadiusLatitudeLongitude;
using geometry::Vector;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::QuinlanTremaine1990Order12;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Degree2SphericalHarmonicCoefficient;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::SIUnit;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using ::testing::IsNull;
using ::testing::NotNull;

class BodyTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  // We need that so the comma doesn't get caught in macros.
  using Direction = Vector<double, World>;

  template<typename Tag, Tag tag>
  void TestRotatingBody() {
    using F = Frame<Tag, tag, true>;

    auto const rotating_body =
        RotatingBody<F>(17 * SIUnit<GravitationalParameter>(),
                        typename RotatingBody<F>::Parameters(
                            2 * Metre,
                            3 * Radian,
                            Instant() + 4 * Second,
                            angular_frequency_,
                            right_ascension_of_pole_,
                            declination_of_pole_));

    serialization::Body message;
    RotatingBody<F> const* cast_rotating_body;
    rotating_body.WriteToMessage(&message);
    EXPECT_TRUE(message.has_massive_body());
    EXPECT_FALSE(message.has_massless_body());
    EXPECT_TRUE(message.massive_body().HasExtension(
                    serialization::RotatingBody::extension));

    not_null<std::unique_ptr<MassiveBody const>> const massive_body =
        MassiveBody::ReadFromMessage(message);
    EXPECT_EQ(rotating_body.gravitational_parameter(),
              massive_body->gravitational_parameter());
    cast_rotating_body =
        dynamic_cast_not_null<RotatingBody<F> const*>(massive_body.get());
    EXPECT_THAT(cast_rotating_body, NotNull());
  }

  AngularFrequency const angular_frequency_ = -1.5 * Radian / Second;
  Angle const right_ascension_of_pole_ = 37 * Degree;
  Angle const declination_of_pole_ = 123 * Degree;

  MasslessBody massless_body_;
  MassiveBody massive_body_ =
      MassiveBody(42 * SIUnit<GravitationalParameter>());
  RotatingBody<World> rotating_body_ =
      RotatingBody<World>(17 * SIUnit<GravitationalParameter>(),
                          RotatingBody<World>::Parameters(
                              1 * Metre,
                              3 * Radian,
                              Instant() + 4 * Second,
                              angular_frequency_,
                              right_ascension_of_pole_,
                              declination_of_pole_));
  OblateBody<World> oblate_body_ =
      OblateBody<World>(17 * SIUnit<GravitationalParameter>(),
                        RotatingBody<World>::Parameters(
                            1 * Metre,
                            3 * Radian,
                            Instant() + 4 * Second,
                            angular_frequency_,
                            right_ascension_of_pole_,
                            declination_of_pole_),
                        OblateBody<World>::Parameters(
                            6,
                            1 * Metre));
};

using BodyDeathTest = BodyTest;

TEST_F(BodyTest, MasslessSerializationSuccess) {
  EXPECT_TRUE(massless_body_.is_massless());
  EXPECT_FALSE(massless_body_.is_oblate());

  serialization::Body message;
  MasslessBody const* cast_massless_body;
  massless_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massless_body());
  EXPECT_FALSE(message.has_massive_body());

  // Direct deserialization.
  // No members to test in this class, we just check that it doesn't crash.
  massless_body_ = *MasslessBody::ReadFromMessage(message);

  // Dispatching from |Body|.  Need two steps to add const and remove
  // |not_null|.
  not_null<std::unique_ptr<Body const>> body = Body::ReadFromMessage(message);
  cast_massless_body = dynamic_cast_not_null<MasslessBody const*>(body.get());
  EXPECT_THAT(cast_massless_body, NotNull());
}

// The best serialization revenge.
TEST_F(BodyTest, MassiveSerializationSuccess) {
  EXPECT_FALSE(massive_body_.is_massless());
  EXPECT_FALSE(massive_body_.is_oblate());

  serialization::Body message;
  MassiveBody const* cast_massive_body;
  massive_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_EQ(42, message.massive_body().gravitational_parameter().magnitude());

  // Direct deserialization.
  MassiveBody const massive_body = *MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(massive_body_.gravitational_parameter(),
            massive_body.gravitational_parameter());

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body>> body = Body::ReadFromMessage(message);
  cast_massive_body = dynamic_cast_not_null<MassiveBody*>(body.get());
  EXPECT_THAT(cast_massive_body, NotNull());
  EXPECT_EQ(massive_body_.gravitational_parameter(),
            cast_massive_body->gravitational_parameter());
}

TEST_F(BodyTest, RotatingSerializationSuccess) {
  EXPECT_FALSE(rotating_body_.is_massless());
  EXPECT_FALSE(rotating_body_.is_oblate());

  serialization::Body message;
  RotatingBody<World> const* cast_rotating_body;
  rotating_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_TRUE(message.massive_body().HasExtension(
                  serialization::RotatingBody::extension));
  EXPECT_EQ(17, message.massive_body().gravitational_parameter().magnitude());
  serialization::RotatingBody const rotating_body_extension =
      message.massive_body().GetExtension(
          serialization::RotatingBody::extension);
  EXPECT_EQ(3, rotating_body_extension.reference_angle().magnitude());
  EXPECT_EQ(4,
            rotating_body_extension.reference_instant().scalar().magnitude());
  EXPECT_EQ(angular_frequency_,
            AngularFrequency::ReadFromMessage(
                rotating_body_extension.angular_frequency()));
  EXPECT_EQ(right_ascension_of_pole_,
            Angle::ReadFromMessage(
                rotating_body_extension.right_ascension_of_pole()));
  EXPECT_EQ(declination_of_pole_,
            Angle::ReadFromMessage(
                rotating_body_extension.declination_of_pole()));

  // Dispatching from |MassiveBody|.
  not_null<std::unique_ptr<MassiveBody const>> const massive_body =
      MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(rotating_body_.gravitational_parameter(),
            massive_body->gravitational_parameter());
  cast_rotating_body =
      dynamic_cast_not_null<RotatingBody<World> const*>(massive_body.get());
  EXPECT_THAT(cast_rotating_body, NotNull());
  EXPECT_EQ(rotating_body_.gravitational_parameter(),
            cast_rotating_body->gravitational_parameter());
  EXPECT_EQ(rotating_body_.angular_velocity(),
            cast_rotating_body->angular_velocity());
  EXPECT_EQ(rotating_body_.AngleAt(Instant()),
            cast_rotating_body->AngleAt(Instant()));

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_rotating_body =
      dynamic_cast_not_null<RotatingBody<World> const*>(body.get());
  EXPECT_THAT(cast_rotating_body, NotNull());
  EXPECT_EQ(rotating_body_.gravitational_parameter(),
            cast_rotating_body->gravitational_parameter());
  EXPECT_EQ(rotating_body_.angular_velocity(),
            cast_rotating_body->angular_velocity());
  EXPECT_EQ(rotating_body_.AngleAt(Instant()),
            cast_rotating_body->AngleAt(Instant()));
}

TEST_F(BodyTest, OblateSerializationSuccess) {
  EXPECT_FALSE(oblate_body_.is_massless());
  EXPECT_TRUE(oblate_body_.is_oblate());

  serialization::Body message;
  OblateBody<World> const* cast_oblate_body;
  oblate_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_TRUE(message.massive_body().GetExtension(
                  serialization::RotatingBody::extension).
                      HasExtension(serialization::OblateBody::extension));
  EXPECT_EQ(17, message.massive_body().gravitational_parameter().magnitude());
  serialization::OblateBody const oblate_body_extension =
      message.massive_body().GetExtension(
                  serialization::RotatingBody::extension).
                      GetExtension(serialization::OblateBody::extension);
  EXPECT_EQ(6, oblate_body_extension.j2());

  // Dispatching from |MassiveBody|.
  not_null<std::unique_ptr<MassiveBody const>> const massive_body =
      MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            massive_body->gravitational_parameter());
  cast_oblate_body =
      dynamic_cast_not_null<OblateBody<World> const*>(massive_body.get());
  EXPECT_THAT(cast_oblate_body, NotNull());
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.polar_axis(), cast_oblate_body->polar_axis());

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_oblate_body =
      dynamic_cast_not_null<OblateBody<World> const*>(body.get());
  EXPECT_THAT(cast_oblate_body, NotNull());
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.polar_axis(), cast_oblate_body->polar_axis());
}

TEST_F(BodyTest, OblateSerializationCompatibility) {
  EXPECT_FALSE(oblate_body_.is_massless());
  EXPECT_TRUE(oblate_body_.is_oblate());

  // Construct a pre-Διόφαντος message.
  serialization::Body message;
  OblateBody<World> const* cast_oblate_body;
  oblate_body_.WriteToMessage(&message);
  serialization::OblateBody* const oblate_body_extension =
      message.mutable_massive_body()->MutableExtension(
                  serialization::RotatingBody::extension)->
                      MutableExtension(serialization::OblateBody::extension);
  oblate_body_extension->clear_reference_radius();
  Degree2SphericalHarmonicCoefficient pre_διόφαντος_j2 =
      7 * SIUnit<Degree2SphericalHarmonicCoefficient>();
  pre_διόφαντος_j2.WriteToMessage(
      oblate_body_extension->mutable_pre_diophantos_j2());

  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_oblate_body =
      dynamic_cast_not_null<OblateBody<World> const*>(body.get());
  Length const reference_radius = 1 * Metre;
  EXPECT_EQ(7 / (cast_oblate_body->gravitational_parameter() /
                 SIUnit<GravitationalParameter>()),
            cast_oblate_body->j2());
  EXPECT_EQ(reference_radius, cast_oblate_body->reference_radius());
  EXPECT_EQ(7 * SIUnit<Degree2SphericalHarmonicCoefficient>() /
                cast_oblate_body->gravitational_parameter(),
            cast_oblate_body->j2_over_μ());
}

TEST_F(BodyTest, AllFrames) {
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::ALICE_SUN>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::ALICE_WORLD>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::BARYCENTRIC>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::NAVIGATION>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::WORLD>();
  TestRotatingBody<serialization::Frame::PluginTag,
                   serialization::Frame::WORLD_SUN>();

  TestRotatingBody<serialization::Frame::SolarSystemTag,
                   serialization::Frame::GCRS>();
  TestRotatingBody<serialization::Frame::SolarSystemTag,
                   serialization::Frame::ICRS>();
  TestRotatingBody<serialization::Frame::SolarSystemTag,
                   serialization::Frame::ITRS>();

  TestRotatingBody<serialization::Frame::TestTag, serialization::Frame::TEST>();
  TestRotatingBody<serialization::Frame::TestTag,
                   serialization::Frame::TEST1>();
  TestRotatingBody<serialization::Frame::TestTag,
                   serialization::Frame::TEST2>();
  TestRotatingBody<serialization::Frame::TestTag, serialization::Frame::FROM>();
  TestRotatingBody<serialization::Frame::TestTag,
                   serialization::Frame::THROUGH>();
  TestRotatingBody<serialization::Frame::TestTag, serialization::Frame::TO>();
}

TEST_F(BodyTest, SurfaceFrame) {
  struct SurfaceFrame;
  const Vector<double, SurfaceFrame> x({1, 0, 0});
  const Bivector<double, SurfaceFrame> z({0, 0, 1});

  SolarSystem<ICRS> solar_system_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const ephemeris = solar_system_j2000.MakeEphemeris(
      /*fitting_tolerance=*/5 * Milli(Metre),
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
          /*step=*/10 * Minute));
  ephemeris->Prolong(J2000 + 1 * Day);

  auto const earth = solar_system_j2000.rotating_body(*ephemeris, "Earth");
  auto const sun = solar_system_j2000.rotating_body(*ephemeris, "Sun");
  auto const earth_trajectory = ephemeris->trajectory(earth);
  auto const sun_trajectory = ephemeris->trajectory(sun);

  auto const from_surface_frame_j2000 =
      earth->FromSurfaceFrame<SurfaceFrame>(J2000);
  auto const earth_centre_j2000 = earth_trajectory->EvaluatePosition(J2000);
  auto const sun_centre_j2000 = sun_trajectory->EvaluatePosition(J2000);
  auto const x_j2000 = from_surface_frame_j2000(x);
  auto const earth_sun_j2000 = sun_centre_j2000 - earth_centre_j2000;
  Displacement<ICRS> const projected_earth_sun_j2000(
      {earth_sun_j2000.coordinates().x,
       earth_sun_j2000.coordinates().y,
       0 * Metre});
  EXPECT_THAT(AngleBetween(x_j2000, projected_earth_sun_j2000),
              IsNear(1.1 * Degree));

  Instant const j2000_3h = J2000 + 3 * Hour;
  auto const from_surface_frame_j2000_3h =
      earth->FromSurfaceFrame<SurfaceFrame>(j2000_3h);
  auto const earth_centre_j2000_3h =
      earth_trajectory->EvaluatePosition(j2000_3h);
  auto const sun_centre_j2000_3h = sun_trajectory->EvaluatePosition(j2000_3h);
  auto const x_j2000_3h = from_surface_frame_j2000_3h(x);
  auto const earth_sun_j2000_3h = sun_centre_j2000_3h - earth_centre_j2000_3h;
  EXPECT_THAT(
      OrientedAngleBetween(
          earth_sun_j2000_3h, x_j2000_3h, from_surface_frame_j2000_3h(z)),
      IsNear(48.3 * Degree));
}

}  // namespace internal_body
}  // namespace physics
}  // namespace principia
