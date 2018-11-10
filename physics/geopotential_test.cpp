
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
using base::make_not_null_unique;
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
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::VanishesBefore;
using ::testing::AllOf;
using ::testing::An;
using ::testing::Each;
using ::testing::Eq;
using ::testing::ElementsAre;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Property;

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
  Geopotential<World> const geopotential(&body, /*tolerance=*/0);

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
    EXPECT_THAT(acceleration1, AlmostEquals(acceleration2, 3));
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
  Geopotential<World> const geopotential(&body, /*tolerance=*/0);

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
  Geopotential<World> const geopotential(&body, /*tolerance=*/0);

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
                AlmostEquals(0.75, 0));
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
  Geopotential<ICRS> const geopotential(&earth, /*tolerance=*/0);

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

TEST_F(GeopotentialTest, HarmonicDamping) {
  HarmonicDamping σ(1 * Metre);
  EXPECT_THAT(σ.inner_threshold(), Eq(1 * Metre));
  EXPECT_THAT(σ.outer_threshold(), Eq(3 * Metre));
  Vector<double, World> x({1, 0, 0});
  Inverse<Square<Length>> const ℜ_over_r = 5 / Pow<2>(Metre);
  Inverse<Square<Length>> const ℜʹ = 17 / Pow<2>(Metre);
  Inverse<Square<Length>> σℜ_over_r;
  Vector<Inverse<Square<Length>>, World> grad_σℜ;

  {
    Length const r = 3 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    EXPECT_THAT(σℜ_over_r, Eq(0 / Pow<2>(Metre)));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(0 / Pow<2>(Metre)));
  }
  {
    Length const r = 2 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    auto const ℜ = ℜ_over_r * r;
    auto const σ = 0.5;
    auto const σʹ = -3 / (4 * Metre);
    EXPECT_THAT(σℜ_over_r, Eq(ℜ_over_r / 2));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(σʹ * ℜ + ℜʹ * σ));
  }
  {
    Length const r = 1 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    EXPECT_THAT(σℜ_over_r, Eq(ℜ_over_r));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(ℜʹ));
  }
  {
    Length const r = 0.5 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    EXPECT_THAT(σℜ_over_r, Eq(ℜ_over_r));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(ℜʹ));
  }
}

TEST_F(GeopotentialTest, ThresholdComputation) {
  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto earth_message = solar_system_2000.gravity_model_message("Earth");
  earth_message.mutable_geopotential()->set_max_degree(5);
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
  auto earth = make_not_null_unique<OblateBody<ICRS>>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  Geopotential<ICRS> geopotential(earth.get(), /*tolerance=*/0x1p-24);

  EXPECT_THAT(
      geopotential.degree_damping(),
      ElementsAre(
          /*0=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*1=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*2=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(1'500'000 * Kilo(Metre))),
          /*3=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(43'000 * Kilo(Metre))),
          /*4=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(23'000 * Kilo(Metre))),
          /*5=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(18'000 * Kilo(Metre)))));
  EXPECT_THAT(geopotential.tesseral_damping().inner_threshold(),
              IsNear(110'000 * Kilo(Metre)));
  EXPECT_THAT(geopotential.first_tesseral_degree(), Eq(3));

  geopotential = Geopotential<ICRS>(earth.get(), /*tolerance=*/0);

  EXPECT_THAT(geopotential.degree_damping(),
              Each(Property(&HarmonicDamping::inner_threshold,
                            Eq(Infinity<Length>()))));
  EXPECT_THAT(geopotential.tesseral_damping().inner_threshold(),
              Eq(Infinity<Length>()));
  EXPECT_THAT(geopotential.first_tesseral_degree(), Eq(0));

  // TODO(eggrobin): This is brittle; we should have |SolarSystem| utilities for
  // that.
  double const earth_c20 = earth_message.geopotential().row(0).column(0).cos();
  earth_message.mutable_geopotential()
      ->mutable_row(0)
      ->mutable_column(0)
      ->clear_cos();
  earth = make_not_null_unique<OblateBody<ICRS>>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  geopotential = Geopotential<ICRS>(earth.get(), /*tolerance=*/0x1p-24);

  EXPECT_THAT(
      geopotential.degree_damping(),
      ElementsAre(
          /*0=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*1=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*2=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(110'000 * Kilo(Metre))),
          /*3=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(43'000 * Kilo(Metre))),
          /*4=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(23'000 * Kilo(Metre))),
          /*5=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(18'000 * Kilo(Metre)))));
  EXPECT_THAT(geopotential.tesseral_damping().inner_threshold(),
              IsNear(110'000 * Kilo(Metre)));
  EXPECT_THAT(geopotential.first_tesseral_degree(), Eq(2));

  earth_message.mutable_geopotential()
      ->mutable_row(0)
      ->mutable_column(0)
      ->set_cos(earth_c20);
  double const earth_c30 = earth_message.geopotential().row(1).column(0).cos();
  earth_message.mutable_geopotential()
      ->mutable_row(1)
      ->mutable_column(0)
      ->set_cos(earth_c20);
  earth = make_not_null_unique<OblateBody<ICRS>>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  geopotential = Geopotential<ICRS>(earth.get(), /*tolerance=*/0x1p-24);

  EXPECT_THAT(
      geopotential.degree_damping(),
      ElementsAre(
          /*0=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*1=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*2=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(1'500'000 * Kilo(Metre))),
          /*3=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(280'000 * Kilo(Metre))),
          /*4=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(23'000 * Kilo(Metre))),
          /*5=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(18'000 * Kilo(Metre)))));
  EXPECT_THAT(geopotential.tesseral_damping().inner_threshold(),
              IsNear(110'000 * Kilo(Metre)));
  EXPECT_THAT(geopotential.first_tesseral_degree(), Eq(4));

  earth_message.mutable_geopotential()
      ->mutable_row(1)
      ->mutable_column(0)
      ->set_cos(earth_c30);
  for (auto& row : *earth_message.mutable_geopotential()->mutable_row()) {
    for (auto& column : *row.mutable_column()) {
      if (column.order() != 0) {
        column.clear_cos();
        column.clear_sin();
      }
    }
  }
  earth = make_not_null_unique<OblateBody<ICRS>>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  geopotential = Geopotential<ICRS>(earth.get(), /*tolerance=*/0x1p-24);

  EXPECT_THAT(
      geopotential.degree_damping(),
      ElementsAre(
          /*0=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*1=*/Property(&HarmonicDamping::inner_threshold,
                         Eq(Infinity<Length>())),
          /*2=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(1'500'000 * Kilo(Metre))),
          /*3=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(35'000 * Kilo(Metre))),
          /*4=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(22'000 * Kilo(Metre))),
          /*5=*/Property(&HarmonicDamping::inner_threshold,
                         IsNear(12'000 * Kilo(Metre)))));
  EXPECT_THAT(geopotential.tesseral_damping().inner_threshold(), Eq(0 * Metre));
  EXPECT_THAT(geopotential.first_tesseral_degree(), Eq(6));

  geopotential = Geopotential<ICRS>(earth.get(), /*tolerance=*/0);

  EXPECT_THAT(geopotential.degree_damping(),
              Each(Property(&HarmonicDamping::inner_threshold,
                            Eq(Infinity<Length>()))));
  EXPECT_THAT(geopotential.tesseral_damping().inner_threshold(), Eq(0 * Metre));
  EXPECT_THAT(geopotential.first_tesseral_degree(), Eq(6));
}

TEST_F(GeopotentialTest, DampedForces) {
  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const earth_message = [&solar_system_2000]() {
    auto message = solar_system_2000.gravity_model_message("Earth");
    message.mutable_geopotential()->set_max_degree(5);
    message.mutable_geopotential()->clear_zonal();
    return message;
  }();

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
  OblateBody<ICRS> const earth(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  Geopotential<ICRS> const earth_geopotential(&earth, /*tolerance=*/0x1p-24);

  // The bodies underlying the geopotentials below.
  std::vector<not_null<std::unique_ptr<OblateBody<ICRS>>>> earths;

  // |geopotential_degree[n]| has the degree n terms of the geopotential, with
  // tolerance 0.
  std::array<std::optional<Geopotential<ICRS>>, 6> geopotential_degree;
  // |geopotential_j2| has the degree 2 zonal term term with tolerance 0.
  std::optional<Geopotential<ICRS>> geopotential_j2;
  // |geopotential_c22_s22| has degree 2 sectoral terms with tolerance 0.
  std::optional<Geopotential<ICRS>> geopotential_c22_s22;

  for (int n = 2; n <= 5; ++n) {
    auto earth_degree_n_message = earth_message;
    for (auto& row :
         *earth_degree_n_message.mutable_geopotential()->mutable_row()) {
      if (row.degree() != n) {
        row.clear_column();
      }
    }
    earths.push_back(make_not_null_unique<OblateBody<ICRS>>(
        massive_body_parameters,
        rotating_body_parameters,
        OblateBody<ICRS>::Parameters::ReadFromMessage(
            earth_degree_n_message.geopotential(), earth_reference_radius)));
    geopotential_degree[n] =
        Geopotential<ICRS>(earths.back().get(), /*tolerance=*/0);
  }
  {
    auto earth_c22_s22_message = earth_message;
    for (auto& row :
         *earth_c22_s22_message.mutable_geopotential()->mutable_row()) {
      for (auto& column : *row.mutable_column()) {
        if (column.order() == 0 || row.degree() != 2) {
          column.clear_cos();
          column.clear_sin();
        }
      }
    }
    earths.push_back(make_not_null_unique<OblateBody<ICRS>>(
        massive_body_parameters,
        rotating_body_parameters,
        OblateBody<ICRS>::Parameters::ReadFromMessage(
            earth_c22_s22_message.geopotential(), earth_reference_radius)));
    geopotential_c22_s22 =
        Geopotential<ICRS>(earths.back().get(), /*tolerance=*/0);
  }
  {
    auto earth_j2_message = earth_message;
    for (auto& row : *earth_j2_message.mutable_geopotential()->mutable_row()) {
      for (auto& column : *row.mutable_column()) {
        if (column.order() != 0 || row.degree() != 2) {
          column.clear_cos();
          column.clear_sin();
        }
      }
    }
    earths.push_back(make_not_null_unique<OblateBody<ICRS>>(
        massive_body_parameters,
        rotating_body_parameters,
        OblateBody<ICRS>::Parameters::ReadFromMessage(
            earth_j2_message.geopotential(), earth_reference_radius)));
    geopotential_j2 = Geopotential<ICRS>(earths.back().get(), /*tolerance=*/0);
}

  Vector<double, ICRS> direction({1, 2, 5});
  direction = direction / direction.Norm();

  auto const get_acceleration =
      [&direction](Geopotential<ICRS> const& geopotential, Length const& r) {
        return geopotential.GeneralSphericalHarmonicsAcceleration(
            Instant(), r * direction, r, r * r, 1 / Pow<3>(r));
      };
  auto const get_radial_acceleration =
      [&direction, &get_acceleration](Geopotential<ICRS> const& geopotential, Length const& r) {
        Vector<double, ICRS> const down = -direction;
        return InnerProduct(get_acceleration(geopotential, r), down);
      };
  auto const get_latitudinal_acceleration =
      [&direction, &get_acceleration](Geopotential<ICRS> const& geopotential,
                                      Length const& r) {
        Vector<double, ICRS> const celestial_north({0, 0, 1});
        Vector<double, ICRS> local_north =
            celestial_north.OrthogonalizationAgainst(direction);
        return InnerProduct(get_acceleration(geopotential, r), local_north);
      };
  auto const get_longitudinal_acceleration =
      [&direction, &get_acceleration](Geopotential<ICRS> const& geopotential,
                                      Length const& r) {
        Vector<double, ICRS> const down = -direction;
        Bivector<double, ICRS> const north({0, 0, 1});
        Vector<double, ICRS> local_east = down * north;
        return InnerProduct(get_acceleration(geopotential, r), local_east);
      };

  // Above the  outer threshold for J2.
  EXPECT_THAT(
      get_acceleration(earth_geopotential, 5'000'000 * Kilo(Metre)).Norm(),
      Eq(0 / Pow<2>(Metre)));

  {
    // Inspect the J2 sigmoid.
    Length const s0 = earth_geopotential.degree_damping()[2].inner_threshold();
    Length const s1 = earth_geopotential.degree_damping()[2].outer_threshold();
    EXPECT_THAT(s0, IsNear(1'500'000 * Kilo(Metre)));

    // The radial component grows beyond the undamped one.  We check the ratio
    // at the arithmetic mean of the thresholds, and at its maximum.
    EXPECT_THAT(
        get_radial_acceleration(earth_geopotential, (s0 + s1) / 2) /
            get_radial_acceleration(*geopotential_j2, (s0 + s1) / 2),
        AlmostEquals(1, 0));
    EXPECT_THAT(
        get_radial_acceleration(earth_geopotential, 2 * s0 * s1 / (s0 + s1)) /
            get_radial_acceleration(*geopotential_j2, 2 * s0 * s1 / (s0 + s1)),
        AlmostEquals(1.125, 2));

    // The latitudinal component remains below the undamped one (it is simply
    // multiplied by σ, instead of involving σ′).
    EXPECT_THAT(
        get_latitudinal_acceleration(earth_geopotential, (s0 + s1) / 2) /
            get_latitudinal_acceleration(*geopotential_j2, (s0 + s1) / 2),
        AlmostEquals(0.5, 2));
  }

  // Below the inner threshold for J2, but still above all other outer thresholds.
  EXPECT_THAT(earth_geopotential.degree_damping()[2].inner_threshold(), Gt(1'000'000 * Kilo(Metre)));
  EXPECT_THAT(earth_geopotential.tesseral_damping().outer_threshold(), Lt(1'000'000 * Kilo(Metre)));
  EXPECT_THAT(earth_geopotential.degree_damping()[3].outer_threshold(), Lt(1'000'000 * Kilo(Metre)));
  EXPECT_THAT(get_acceleration(earth_geopotential, 1'000'000 * Kilo(Metre)),
              Eq(get_acceleration(*geopotential_j2, 1'000'000 * Kilo(Metre))));

  {
    // Inspect the C22 and S22 sigmoid.
    Length const s0 = earth_geopotential.tesseral_damping().inner_threshold();
    Length const s1 = earth_geopotential.tesseral_damping().outer_threshold();
    EXPECT_THAT(s0, IsNear(101'000 * Kilo(Metre)));

    // Although this sigmoid overlaps with the degree 3 one, the midpoint still
    // lies above the outer threshold for degree 3.
    EXPECT_THAT(s0,
                Lt(earth_geopotential.degree_damping()[3].outer_threshold()));
    EXPECT_THAT((s0 + s1) / 2,
                Gt(earth_geopotential.degree_damping()[3].outer_threshold()));

    EXPECT_THAT(
        (get_radial_acceleration(earth_geopotential, (s0 + s1) / 2) -
         get_radial_acceleration(*geopotential_j2, (s0 + s1) / 2)) /
            get_radial_acceleration(*geopotential_c22_s22, (s0 + s1) / 2),
        AlmostEquals(1, 212));
    EXPECT_THAT(
        (get_latitudinal_acceleration(earth_geopotential, (s0 + s1) / 2) -
         get_latitudinal_acceleration(*geopotential_j2, (s0 + s1) / 2)) /
            get_latitudinal_acceleration(*geopotential_c22_s22, (s0 + s1) / 2),
        AlmostEquals(0.5, 1773));
    EXPECT_THAT(
        get_longitudinal_acceleration(earth_geopotential, (s0 + s1) / 2) /
            get_longitudinal_acceleration(*geopotential_c22_s22, (s0 + s1) / 2),
        AlmostEquals(0.5, 103));
  }

  {
    // Inspect the degree 3 sigmoid.
    Length const s0 = earth_geopotential.degree_damping()[3].inner_threshold();
    Length const s1 = earth_geopotential.degree_damping()[3].outer_threshold();
    EXPECT_THAT(s0, IsNear(43'000 * Kilo(Metre)));

    // Although this sigmoid overlaps with the degree 3 and tesseral ones, the
    // midpoint still lies above the outer threshold for degree 3, and below the
    // inner tesseral threshold.
    EXPECT_THAT(s0,
                Lt(earth_geopotential.degree_damping()[4].outer_threshold()));
    EXPECT_THAT(s1, Gt(earth_geopotential.tesseral_damping().inner_threshold()));
    EXPECT_THAT(
        (s0 + s1) / 2,
        AllOf(Gt(earth_geopotential.degree_damping()[4].outer_threshold()),
              Lt(earth_geopotential.tesseral_damping().inner_threshold())));

    // Note that the radial acceleration ratio at the midpoint is 7/8 rather
    // than 1, as it depends on the degree: it is (4 + n) / (2n + 2).
    EXPECT_THAT(
        (get_radial_acceleration(earth_geopotential, (s0 + s1) / 2) -
         get_radial_acceleration(*geopotential_degree[2], (s0 + s1) / 2)) /
            get_radial_acceleration(*geopotential_degree[3], (s0 + s1) / 2),
        AlmostEquals(0.875, 278));

    EXPECT_THAT(
        (get_latitudinal_acceleration(earth_geopotential, (s0 + s1) / 2) -
         get_latitudinal_acceleration(*geopotential_degree[2], (s0 + s1) / 2)) /
            get_latitudinal_acceleration(*geopotential_degree[3], (s0 + s1) / 2),
        AlmostEquals(0.5, 26512));
    EXPECT_THAT(
        (get_longitudinal_acceleration(earth_geopotential, (s0 + s1) / 2) -
         get_longitudinal_acceleration(*geopotential_degree[2],
                                       (s0 + s1) / 2)) /
            get_longitudinal_acceleration(*geopotential_degree[3],
                                          (s0 + s1) / 2),
        AlmostEquals(0.5, 2130));
  }

  // The outer threshold for degree 5 is above the inner threshold for degree 3,
  // so degrees 4 and above have mixed sigmoids, which are tricky to test.
    EXPECT_THAT(
        earth_geopotential.degree_damping()[5].outer_threshold(),
        Gt(earth_geopotential.degree_damping()[3].inner_threshold()));
}

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
