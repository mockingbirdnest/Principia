
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
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using astronomy::ICRS;
using astronomy::ITRS;
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
using testing_utilities::RelativeError;
using testing_utilities::VanishesBefore;
using ::testing::An;
using ::testing::Gt;
using ::testing::Lt;

class GeopotentialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  GeopotentialTest()
      : massive_body_parameters_(17 * SIUnit<GravitationalParameter>()),
        rotating_body_parameters_(1 * Metre,
                                  3 * Radian,
                                  Instant() + 4 * Second,
                                  angular_frequency_,
                                  right_ascension_of_pole_,
                                  declination_of_pole_) {}

  template<typename Frame>
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  SphericalHarmonicsAcceleration(Geopotential<Frame> const& geopotential,
                                 Instant const& t,
                                 Displacement<Frame> const& r) {
    auto const r² = r.Norm²();
    auto const one_over_r³ = 1.0 / (r² * r.Norm());
    return geopotential.SphericalHarmonicsAcceleration(t, r, r², one_over_r³);
  }

  template<typename Frame>
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  GeneralSphericalHarmonicsAcceleration(Geopotential<Frame> const& geopotential,
                                        Instant const& t,
                                        Displacement<Frame> const& r) {
    auto const r² = r.Norm²();
    auto const r_norm = Sqrt(r²);
    auto const one_over_r³ = r_norm / (r² * r²);
    return geopotential.GeneralSphericalHarmonicsAcceleration(
        t, r, r_norm, r², one_over_r³);
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

  // Consistency between the general implementation in the zonal case and the
  // J2-specific one.
  {
    auto const acceleration1 = SphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({6 * Metre, -4 * Metre, 5 * Metre}));
    auto const acceleration2 = GeneralSphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({6 * Metre, -4 * Metre, 5 * Metre}));
    EXPECT_THAT(acceleration1, AlmostEquals(acceleration2, 1));
  }
}

TEST_F(GeopotentialTest, C22S22) {
  serialization::OblateBody::Geopotential message;
  {
    auto* const degree2 = message.add_row();
    degree2->set_degree(2);
    auto* const order0 = degree2->add_column();
    order0->set_order(0);
    order0->set_cos(-6 / LegendreNormalizationFactor(2, 0));
    order0->set_sin(0);
    auto* const order2 = degree2->add_column();
    order2->set_order(2);
    order2->set_cos(10 / LegendreNormalizationFactor(2, 2));
    order2->set_sin(-13 / LegendreNormalizationFactor(2, 2));
  }
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters::ReadFromMessage(
                            message, 1 * Metre));
  Geopotential<World> const geopotential(&body);

  // The acceleration at a point located on the axis is along the axis for the
  // (2, 2) harmonics.
  {
    auto const acceleration = GeneralSphericalHarmonicsAcceleration(
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
    auto const acceleration = GeneralSphericalHarmonicsAcceleration(
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
  serialization::OblateBody::Geopotential message;
  {
    auto* const degree2 = message.add_row();
    degree2->set_degree(2);
    auto* const order0 = degree2->add_column();
    order0->set_order(0);
    order0->set_cos(-6 / LegendreNormalizationFactor(2, 0));
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
    order0->set_cos(5 / LegendreNormalizationFactor(3, 0));
  }
  OblateBody<World> const body =
      OblateBody<World>(massive_body_parameters_,
                        rotating_body_parameters_,
                        OblateBody<World>::Parameters::ReadFromMessage(
                            message, 1 * Metre));
  Geopotential<World> const geopotential(&body);

  // The acceleration at a point located on the axis is along the axis.
  {
    auto const acceleration = GeneralSphericalHarmonicsAcceleration(
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
    auto const acceleration = GeneralSphericalHarmonicsAcceleration(
        geopotential,
        Instant(),
        Displacement<World>({30 * Metre, 40 * Metre, 0 * Metre}));
    EXPECT_THAT(acceleration.coordinates().x / acceleration.coordinates().y,
                AlmostEquals(0.75, 1));
    EXPECT_THAT(acceleration.coordinates().z,
                Not(VanishesBefore(1 * Pow<-2>(Metre), 0)));
    EXPECT_THAT(acceleration.coordinates().z, Lt(0 * Pow<-2>(Metre)));
  }
}

TEST_F(GeopotentialTest, TestVector) {
  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto earth_message = solar_system_2000.gravity_model_message("Earth");
  earth_message.mutable_geopotential()->set_max_degree(9);
  earth_message.mutable_geopotential()->clear_zonal();

  auto const earth_μ = solar_system_2000.gravitational_parameter("Earth");
  auto const earth_reference_radius =
      ParseQuantity<Length>(earth_message.reference_radius());
  MassiveBody::Parameters const massive_body_parameters(earth_μ);
  RotatingBody<ICRS>::Parameters rotating_body_parameters(
      /*mean_radius=*/solar_system_2000.mean_radius("Earth"),
      /*reference_angle=*/0 * Radian,
      /*reference_instant=*/Instant(),
      /*angular_frequency=*/1e-20 * Radian / Second,
      right_ascension_of_pole_,
      declination_of_pole_);
  OblateBody<ICRS> const earth = OblateBody<ICRS>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  Geopotential<ICRS> const geopotential(&earth);

  // This test vector is from Kuga & Carrara, "Fortran- and C-codes for higher
  // order and degree geopotential computation",
  // www.dem.inpe.br/~hkk/software/high_geopot.html.  They use EGM2008, so their
  // actual results are expected to differ a bit from ours.
  Displacement<ITRS> const displacement(
      {6000000 * Metre, -4000000 * Metre, 5000000 * Metre});
  Vector<Acceleration, ITRS> const expected_acceleration(
      {-3.5377058876337 * Metre / Second / Second,
       2.3585194144421 * Metre / Second / Second,
       -2.9531441870790 * Metre / Second / Second});

  Vector<Acceleration, ITRS> actual_acceleration_cpp;
  Vector<Acceleration, ITRS> actual_acceleration_f90;

  {
    Displacement<ICRS> const icrs_displacement =
        earth.FromSurfaceFrame<ITRS>(Instant())(
            Displacement<ITRS>(displacement));
    auto const icrs_acceleration =
        earth_μ * (GeneralSphericalHarmonicsAcceleration(
                       geopotential, Instant(), icrs_displacement) -
                   icrs_displacement / Pow<3>(icrs_displacement.Norm()));
    actual_acceleration_cpp =
        earth.ToSurfaceFrame<ITRS>(Instant())(icrs_acceleration);
    EXPECT_THAT(RelativeError(actual_acceleration_cpp, expected_acceleration),
                IsNear(1.3e-8));
  }
  {
    double mu = earth_μ / SIUnit<GravitationalParameter>();
    double rbar = earth_reference_radius / Metre;
    numerics::FixedMatrix<double, 10, 10> cnm;
    numerics::FixedMatrix<double, 10, 10> snm;
    for (int n = 0; n <= 9; ++n) {
      for (int m = 0; m <= n; ++m) {
        cnm[n][m] = earth.cos()[n][m] * LegendreNormalizationFactor(n, m);
        snm[n][m] = earth.sin()[n][m] * LegendreNormalizationFactor(n, m);
      }
    }
    actual_acceleration_f90 = Vector<Acceleration, ITRS>(
        SIUnit<Acceleration>() *
        astronomy::fortran_astrodynamics_toolkit::
            ComputeGravityAccelerationLear<9, 9>(
                displacement.coordinates() / Metre, mu, rbar, cnm, snm));
    EXPECT_THAT(RelativeError(actual_acceleration_f90, expected_acceleration),
                IsNear(1.3e-8));
  }
  EXPECT_THAT(actual_acceleration_cpp,
              AlmostEquals(actual_acceleration_f90, 4, 8));
}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
