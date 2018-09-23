
#include "physics/geopotential.hpp"

#include "astronomy/fortran_astrodynamics_toolkit.hpp"
#include "astronomy/frames.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/legendre.hpp"
#include "physics/solar_system.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using astronomy::ICRS;
using geometry::Frame;
using numerics::LegendreNormalizationFactor;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Degree2SphericalHarmonicCoefficient;
using quantities::Degree3SphericalHarmonicCoefficient;
using quantities::ParseQuantity;
using quantities::Pow;
using quantities::SIUnit;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::IsNear;
using testing_utilities::VanishesBefore;
using ::testing::An;
using ::testing::Gt;
using ::testing::Lt;

class GeopotentialTest : public ::testing::Test {
 protected:
  using Surface = Frame<serialization::Frame::TestTag,
                        serialization::Frame::TEST1, false>;
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST2, true>;

  GeopotentialTest()
      : massive_body_parameters_(17 * SIUnit<GravitationalParameter>()),
        rotating_body_parameters_(1 * Metre,
                                  3 * Radian,
                                  Instant() + 4 * Second,
                                  angular_frequency_,
                                  right_ascension_of_pole_,
                                  declination_of_pole_) {}

  Vector<Quotient<Acceleration, GravitationalParameter>, World>
  SphericalHarmonicsAcceleration(Geopotential<World> const& geopotential,
                                 Instant const& t,
                                 Displacement<World> const& r) {
    auto const r² = r.Norm²();
    auto const one_over_r³ = 1.0 / (r² * r.Norm());
    return geopotential.SphericalHarmonicsAcceleration(t, r, r², one_over_r³);
  }

  Vector<Quotient<Acceleration, GravitationalParameter>, World>
  GeneralSphericalHarmonicsAcceleration(Geopotential<World> const& geopotential,
                                     Instant const& t,
                                     Displacement<World> const& r) {
    auto const r² = r.Norm²();
    auto const one_over_r³ = 1.0 / (r² * r.Norm());
    return geopotential.GeneralSphericalHarmonicsAcceleration(
        t, r, r², one_over_r³);
  }

  // The axis of rotation is along the z axis for ease of testing.
  AngularFrequency const angular_frequency_ = -1.5 * Radian / Second;
  Angle const right_ascension_of_pole_ = 0 * Degree;
  Angle const declination_of_pole_ = 90 * Degree;

  MassiveBody::Parameters const massive_body_parameters_;
  RotatingBody<World>::Parameters const rotating_body_parameters_;
};

TEST_F(GeopotentialTest, J2) {
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(/*j2=*/6, 1 * Metre));
  Geopotential<World> const geopotential(&body);

  // The acceleration at a point located on the axis is along the axis.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}));
    EXPECT_THAT(acceleration,
                Componentwise(VanishesBefore(1 * Pow<-2>(Metre), 0),
                              VanishesBefore(1 * Pow<-2>(Metre), 0),
                              An<Exponentiation<Length, -2>>()));
  }

  // The acceleration at a point located in the equatorial plane is directed to
  // the centre.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({30 * Metre, 40 * Metre, 0 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x / acceleration.coordinates().y,
                AlmostEquals(0.75, 0));
    EXPECT_THAT(acceleration.coordinates().z,
                VanishesBefore(1 * Pow<-2>(Metre), 0));
  }

  // The acceleration at a random point nudges the overall force away from the
  // centre and towards the equatorial plane.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({1e2 * Metre, 0 * Metre, 1e2 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x, Gt(0 * Pow<-2>(Metre)));
    EXPECT_THAT(acceleration.coordinates().z, Lt(0 * Pow<-2>(Metre)));
  }
}

TEST_F(GeopotentialTest, C22S22) {
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(
                            /*j2=*/6, /*c22=*/10, /*s22=*/-13, 1 * Metre));
  Geopotential<World> const geopotential(&body);

  // The acceleration at a point located on the axis is along the axis for the
  // (2, 2) harmonics.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}));
    EXPECT_THAT(acceleration,
                Componentwise(VanishesBefore(1 * Pow<-2>(Metre), 0),
                              VanishesBefore(1 * Pow<-2>(Metre), 0),
                              An<Exponentiation<Length, -2>>()));
  }

  // The acceleration at a point located in the equatorial plane is in the plane
  // but not directed to the centre.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({30 * Metre, 40 * Metre, 0 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x / acceleration.coordinates().y,
                Not(IsNear(0.75)));
    EXPECT_THAT(acceleration.coordinates().z,
                VanishesBefore(1 * Pow<-2>(Metre), 0));
  }
}

TEST_F(GeopotentialTest, J3) {
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(
                            /*j2=*/6, /*c22=*/1e-20, /*s22=*/1e-20,
                            /*j3=*/-5, 1 * Metre));
  Geopotential<World> const geopotential(&body);

  // The acceleration at a point located on the axis is along the axis.
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}));
    EXPECT_THAT(acceleration,
                Componentwise(VanishesBefore(1 * Pow<-2>(Metre), 0),
                              VanishesBefore(1 * Pow<-2>(Metre), 0),
                              An<Exponentiation<Length, -2>>()));
  }

  // The acceleration at a point located in the equatorial plane points towards
  // the north, as it does on Earth (I think).
  // TODO(phl): I don't know what I think anymore.  Oh the humanity!
  {
    auto const acceleration = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({30 * Metre, 40 * Metre, 0 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x / acceleration.coordinates().y,
                AlmostEquals(0.75, 0));
    EXPECT_THAT(acceleration.coordinates().z,
                Not(VanishesBefore(1 * Pow<-2>(Metre), 0)));
    EXPECT_THAT(acceleration.coordinates().z, Lt(0 * Pow<-2>(Metre)));
  }
}

TEST_F(GeopotentialTest, VerifyJ2) {
  OblateBody<World> const body1 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(/*j2=*/6, 1 * Metre));
  Geopotential<World> const geopotential1(&body1);
  serialization::OblateBody::Geopotential message;
  {
    auto* const degree2 = message.add_row();
    degree2->set_degree(2);
    auto* const order0 = degree2->add_column();
    order0->set_order(0);
    order0->set_cos(-6 / LegendreNormalizationFactor(2, 0));
    order0->set_sin(0);
  }
  OblateBody<World> const body2 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters::ReadFromMessage(
                            message, 1 * Metre));
  Geopotential<World> const geopotential2(&body2);

  // Check that the accelerations computed according to both methods are
  // consistent.
  {
    Displacement<World> const displacement(
        {0 * Metre, 0 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2,
                Componentwise(AlmostEquals(0 / Metre / Metre, 0),
                              AlmostEquals(0 / Metre / Metre, 0),
                              An<Exponentiation<Length, -2>>()));
  }
  {
    Displacement<World> const displacement(
        {1e-5 * Metre, 1e-5 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 0, 182019));
  }
  {
    Displacement<World> const displacement(
        {5 * Metre, 7 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 2, 54));
  }
}

TEST_F(GeopotentialTest, VerifyC22) {
  OblateBody<World> const body1 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(/*j2=*/1e-20,
                                                      /*c22=*/6,
                                                      /*s22=*/1e-20,
                                                      1 * Metre));
  Geopotential<World> const geopotential1(&body1);
  serialization::OblateBody::Geopotential message;
  {
    auto* const degree2 = message.add_row();
    degree2->set_degree(2);
    auto* const order0 = degree2->add_column();
    order0->set_order(0);
    order0->set_cos(-1e-20 / LegendreNormalizationFactor(2, 0));
    order0->set_sin(0);
    auto* const order2 = degree2->add_column();
    order2->set_order(2);
    order2->set_cos(6 / LegendreNormalizationFactor(2, 2));
    order2->set_sin(1e-20 / LegendreNormalizationFactor(2, 2));
  }
  OblateBody<World> const body2 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters::ReadFromMessage(
                            message, 1 * Metre));
  Geopotential<World> const geopotential2(&body2);

  // Check that the accelerations computed according to both methods are
  // consistent.
  {
    Displacement<World> const displacement(
        {1e-5 * Metre, 1e-5 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 1, 34));
  }
  {
    Displacement<World> const displacement(
        {5 * Metre, 7 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 2, 54));
  }
}

TEST_F(GeopotentialTest, VerifyS22) {
  OblateBody<World> const body1 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(/*j2=*/1e-20,
                                                      /*c22=*/1e-20,
                                                      /*s22=*/6,
                                                      1 * Metre));
  Geopotential<World> const geopotential1(&body1);
  serialization::OblateBody::Geopotential message;
  {
    auto* const degree2 = message.add_row();
    degree2->set_degree(2);
    auto* const order0 = degree2->add_column();
    order0->set_order(0);
    order0->set_cos(-1e-20 / LegendreNormalizationFactor(2, 0));
    order0->set_sin(0);
    auto* const order2 = degree2->add_column();
    order2->set_order(2);
    order2->set_cos(1e-20 / LegendreNormalizationFactor(2, 2));
    order2->set_sin(6 / LegendreNormalizationFactor(2, 2));
  }
  OblateBody<World> const body2 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters::ReadFromMessage(
                            message, 1 * Metre));
  Geopotential<World> const geopotential2(&body2);

  // Check that the accelerations computed according to both methods are
  // consistent.
  {
    Displacement<World> const displacement(
        {1e-5 * Metre, 1e-5 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 0, 14));
  }
  {
    Displacement<World> const displacement(
        {5 * Metre, 7 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 5, 6));
  }
}

TEST_F(GeopotentialTest, VerifyJ3) {
  OblateBody<World> const body1 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters(/*j2=*/1e-20,
                                                      /*c22=*/1e-20,
                                                      /*s22=*/1e-20,
                                                      /*j3=*/6,
                                                      1 * Metre));
  Geopotential<World> const geopotential1(&body1);
  serialization::OblateBody::Geopotential message;
  {
    auto* const degree2 = message.add_row();
    degree2->set_degree(2);
    auto* const order0 = degree2->add_column();
    order0->set_order(0);
    order0->set_cos(-1e-20 / LegendreNormalizationFactor(2, 0));
    order0->set_sin(0);
    auto* const order2 = degree2->add_column();
    order2->set_order(2);
    order2->set_cos(1e-20 / LegendreNormalizationFactor(2, 2));
    order2->set_sin(1e-20 / LegendreNormalizationFactor(2, 2));
  }
  {
    auto* const degree3 = message.add_row();
    degree3->set_degree(3);
    auto* const order0 = degree3->add_column();
    order0->set_order(0);
    order0->set_cos(-6 / LegendreNormalizationFactor(3, 0));
  }
  OblateBody<World> const body2 =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters::ReadFromMessage(
                            message, 1 * Metre));
  Geopotential<World> const geopotential2(&body2);

  // Check that the accelerations computed according to both methods are
  // consistent.
  {
    Displacement<World> const displacement(
        {1e-5 * Metre, 1e-5 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 0, 264755));
  }
  {
    Displacement<World> const displacement(
        {5 * Metre, 7 * Metre, 11 * Metre});
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential1, Instant(), displacement);
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential2, Instant(), displacement);
    EXPECT_THAT(acceleration2, AlmostEquals(acceleration1, 3, 6));
  }
}

TEST_F(GeopotentialTest, VerifyFortran) {
  MassiveBody::Parameters const massive_body_parameters(
      1 * SIUnit<GravitationalParameter>());
  RotatingBody<World>::Parameters rotating_body_parameters(
      /*mean_radius=*/1 * Metre,
      /*reference_angle=*/0 * Radian,
      /*reference_instant=*/Instant(),
      /*angular_frequency=*/1e-20 * Radian / Second,
      right_ascension_of_pole_,
      declination_of_pole_);
  serialization::OblateBody::Geopotential message;
  {
    auto* const degree2 = message.add_row();
    degree2->set_degree(2);
    auto* const order0 = degree2->add_column();
    order0->set_order(0);
    order0->set_cos(0 / LegendreNormalizationFactor(2, 0));
    order0->set_sin(0);
    auto* const order2 = degree2->add_column();
    order2->set_order(2);
    order2->set_cos(0 / LegendreNormalizationFactor(2, 2));
    order2->set_sin(0 / LegendreNormalizationFactor(2, 2));
  }
  {
    auto* const degree3 = message.add_row();
    degree3->set_degree(3);
    auto* const order0 = degree3->add_column();
    order0->set_order(0);
    order0->set_cos(0 / LegendreNormalizationFactor(3, 0));
    auto* const order1 = degree3->add_column();
    order1->set_order(1);
    order1->set_cos(6 / LegendreNormalizationFactor(3, 1));
    order1->set_sin(0 / LegendreNormalizationFactor(3, 1));
    auto* const order2 = degree3->add_column();
    order2->set_order(2);
    order2->set_cos(0 / LegendreNormalizationFactor(3, 2));
    order2->set_sin(0 / LegendreNormalizationFactor(3, 2));
    auto* const order3 = degree3->add_column();
    order3->set_order(3);
    order3->set_cos(0 / LegendreNormalizationFactor(3, 3));
    order3->set_sin(0 / LegendreNormalizationFactor(3, 3));
  }
  {
    auto* const degree4 = message.add_row();
    degree4->set_degree(4);
    auto* const order0 = degree4->add_column();
    order0->set_order(0);
    order0->set_cos(0 / LegendreNormalizationFactor(4, 0));
  }
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters,
                        rotating_body_parameters,
                        OblateBody<World>::Parameters::ReadFromMessage(
                            message, 1 * Metre));
  Geopotential<World> const geopotential(&body);
  {
    geometry::R3Element<double> rgr{5, 0, 0};
    double mu = 1;
    double rbar = 1;
    numerics::FixedMatrix<double, 5, 5> cnm;
    numerics::FixedMatrix<double, 5, 5> snm;
    cnm[2][0] = 0;
    cnm[2][2] = 0;
    snm[2][2] = 0;
    cnm[3][0] = 0;
    cnm[3][1] = 6;
    snm[3][1] = 0;
    cnm[3][2] = 0;
    snm[3][2] = 0;
    cnm[3][3] = 0;
    snm[3][3] = 0;
    cnm[4][0] = 0;
    LOG(ERROR)<<cnm[2][0]<<" "<<cnm[2][1]<<" "<<cnm[2][2];
    LOG(ERROR)<<snm[2][0]<<" "<<snm[2][1]<<" "<<snm[2][2];
    LOG(ERROR)<<cnm[3][0]<<" "<<cnm[3][1]<<" "<<cnm[3][2]<<" "<<cnm[3][3];
    LOG(ERROR)<<snm[3][0]<<" "<<snm[3][1]<<" "<<snm[3][2]<<" "<<snm[3][3];
    auto const acceleration2 = Vector<Acceleration, Surface>(
        1 * SIUnit<Acceleration>() *
        astronomy::fortran_astrodynamics_toolkit::Grav<4, 4>(
            rgr, mu, rbar, cnm, snm));

    Displacement<World> const displacement =
        body.FromSurfaceFrame<Surface>(Instant())(
            Displacement<Surface>(rgr * Metre));
    auto const acceleration1 = 1 * SIUnit<GravitationalParameter>() *
                              (GeneralSphericalHarmonicsAcceleration(
                                   geopotential, Instant(), displacement) -
                               displacement / Pow<3>(displacement.Norm()));
    auto const acceleration3 =
        body.ToSurfaceFrame<Surface>(Instant())(acceleration1);

    auto a = 1 * SIUnit<GravitationalParameter>() *
             GeneralSphericalHarmonicsAcceleration(
                 geopotential, Instant(), displacement);
    auto b = 1 * SIUnit<GravitationalParameter>() *
             (-displacement / Pow<3>(displacement.Norm()));

    EXPECT_THAT(acceleration3, AlmostEquals(acceleration2, 0));

    auto cx =
        (acceleration2.coordinates().x - b.coordinates().x) / a.coordinates().x;
    auto cy =
        (acceleration2.coordinates().y - b.coordinates().y) / a.coordinates().y;
    auto cz =
        (acceleration2.coordinates().z - b.coordinates().z) / a.coordinates().z;
    LOG(ERROR)<<cx<<" "<<cy<<" "<<cz;
  }
}

TEST_F(GeopotentialTest, TestVector) {
  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const earth_message = solar_system_2000.gravity_model_message("Earth");

  auto const earth_μ = solar_system_2000.gravitational_parameter("Earth");
  auto const earth_reference_radius =
      ParseQuantity<Length>(earth_message.reference_radius());
  MassiveBody::Parameters const massive_body_parameters(earth_μ);
  RotatingBody<World>::Parameters rotating_body_parameters(
      /*mean_radius=*/solar_system_2000.mean_radius("Earth"),
      /*reference_angle=*/0 * Radian,
      /*reference_instant=*/Instant(),
      /*angular_frequency=*/1e-20 * Radian / Second,
      right_ascension_of_pole_,
      declination_of_pole_);
  OblateBody<World> const body = OblateBody<World>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<World>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  Geopotential<World> const geopotential(&body);

  geometry::R3Element<double> rgr{6000000, -4000000, 5000000};
  {
    Displacement<World> const displacement =
        body.FromSurfaceFrame<Surface>(Instant())(
            Displacement<Surface>(rgr * Metre));
    auto const acceleration =
        earth_μ * (GeneralSphericalHarmonicsAcceleration(
                       geopotential, Instant(), displacement) -
                   displacement / Pow<3>(displacement.Norm()));
    auto const acceleration2 =
        body.ToSurfaceFrame<Surface>(Instant())(acceleration);
    // Result shoud read: 9  -3.5377058876337  2.3585194144421  -2.9531441870790
    LOG(ERROR)<<acceleration2;
  }
  {
    double mu = earth_μ / SIUnit<GravitationalParameter>();
    double rbar = earth_reference_radius / Metre;
    numerics::FixedMatrix<double, 10, 10> cnm;
    numerics::FixedMatrix<double, 10, 10> snm;
    for (int n = 0; n <= 9; ++n) {
      for (int m = 0; m <= n; ++m) {
        cnm[n][m] = body.cos()[n][m] * LegendreNormalizationFactor(n, m);
        snm[n][m] = body.sin()[n][m] * LegendreNormalizationFactor(n, m);
      }
    }
    auto const acceleration = Vector<Acceleration, World>(
        1 * SIUnit<Acceleration>() *
        astronomy::fortran_astrodynamics_toolkit::Grav<9, 9>(
            rgr, mu, rbar, cnm, snm));
    LOG(ERROR)<<acceleration;
  }
}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
