#include <memory>

#include "astronomy/mercury_orbiter.hpp"
#include "astronomy/orbital_elements.hpp"
#include "astronomy/frames.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/solar_system.hpp"
#include "mathematica/logger.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace astronomy {

using base::not_null;
using geometry::Frame;
using geometry::Instant;
using geometry::Interval;
using geometry::NonRotating;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MassiveBody;
using physics::MasslessBody;
using physics::SolarSystem;
using physics::Trajectory;
using quantities::Cos;
using quantities::Length;
using quantities::Pow;
using quantities::Sin;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::operator""_;
using testing_utilities::IsNear;

// A test that showcases the eccentricity-inclination exchange mechanism
// described in [Лид61] and [Koz62].  We follow the treatment in [Лид61].
class Лидов古在Test : public ::testing::Test {
 protected:
  using MercuryCentredInertial =
      Frame<enum class MercuryCentredInertialTag, NonRotating>;

  Лидов古在Test()
      : solar_system_1950_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2433282_500000000.proto.txt"),
        ephemeris_(solar_system_1950_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<
                    QuinlanTremaine1990Order12,
                    Ephemeris<ICRS>::NewtonianMotionEquation>(),
                /*step=*/10 * Minute))),
        mercury_(*solar_system_1950_.massive_body(*ephemeris_, "Mercury")),
        mercury_frame_(ephemeris_.get(), &mercury_) {
    google::SetStderrLogging(google::INFO);
  }

  ~Лидов古在Test() override {
    google::SetStderrLogging(FLAGS_stderrthreshold);
  }

  SolarSystem<ICRS> solar_system_1950_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> ephemeris_;
  MassiveBody const& mercury_;
  BodyCentredNonRotatingDynamicFrame<ICRS, MercuryCentredInertial>
      mercury_frame_;
};

#if !_DEBUG
TEST_F(Лидов古在Test, MercuryOrbiter) {
  DiscreteTrajectory<ICRS> icrs_trajectory;
  EXPECT_OK(icrs_trajectory.Append(
      MercuryOrbiterInitialTime, MercuryOrbiterInitialDegreesOfFreedom<ICRS>));
  auto& icrs_segment = icrs_trajectory.segments().front();
  icrs_segment.SetDownsampling({.max_dense_intervals = 10'000,
                                .tolerance = 10 * Metre});
  auto const instance =
      ephemeris_->NewInstance({&icrs_trajectory},
                              Ephemeris<ICRS>::NoIntrinsicAccelerations,
                              {SymmetricLinearMultistepIntegrator<
                                   Quinlan1999Order8A,
                                   Ephemeris<ICRS>::NewtonianMotionEquation>(),
                               /*step=*/10 * Second});

  for (int year = 1967; year < 1980; ++year) {
    Instant const t = DateTimeAsTT(date_time::DateTime::BeginningOfDay(
        date_time::Date::Ordinal(year, 1)));
    LOG(INFO) << "Flowing to " << t;
    auto const status = ephemeris_->FlowWithFixedStep(t, *instance);
    if (!status.ok()) {
      LOG(INFO) << status << " at " << icrs_trajectory.back().time;
      break;
    }
  }
  mathematica::Logger logger(
      SOLUTION_DIR / "mathematica" /
          PRINCIPIA_UNICODE_PATH("лидов_古在.generated.wl"),
      /*make_unique=*/false);

  DiscreteTrajectory<MercuryCentredInertial> mercury_centred_trajectory;
  for (auto const& [t, dof] : icrs_trajectory) {
    EXPECT_OK(mercury_centred_trajectory.Append(
        t, mercury_frame_.ToThisFrameAtTime(t)(dof)));
    logger.Append(
        "q",
        mercury_centred_trajectory.back().degrees_of_freedom.position(),
        mathematica::ExpressIn(Metre));
  }

  OrbitalElements const elements = OrbitalElements::ForTrajectory(
      mercury_centred_trajectory, mercury_, MasslessBody{}).value();
  // The constants c₁ and c₂ are defined in [Лид61], equations (58) and (59)
  // respectively.
  Interval<double> c₁;
  Interval<double> c₂;
  for (auto const& elements : elements.mean_elements()) {
    double const ε = 1 - Pow<2>(elements.eccentricity);
    double const cos²_i = Pow<2>(Cos(elements.inclination));
    double const sin²_i = Pow<2>(Sin(elements.inclination));
    double const sin²_ω = Pow<2>(Sin(elements.argument_of_periapsis));
    c₁.Include(ε * cos²_i);
    c₂.Include((1 - ε) * (2.0 / 5.0 - sin²_i * sin²_ω));
    logger.Append("t", elements.time, mathematica::ExpressIn(Second));
    logger.Append("a", elements.semimajor_axis, mathematica::ExpressIn(Metre));
    logger.Append("e", elements.eccentricity);
    logger.Append("i", elements.inclination, mathematica::ExpressIn(Radian));
    logger.Append(R"(\[Omega])",
                  elements.argument_of_periapsis,
                  mathematica::ExpressIn(Radian));
  }
  // The elements e, i, and ω all vary quite a lot.
  EXPECT_THAT(elements.mean_eccentricity_interval().min, IsNear(0.40_(1)));
  EXPECT_THAT(elements.mean_eccentricity_interval().max, IsNear(0.88_(1)));
  EXPECT_THAT(elements.mean_inclination_interval().min,
              IsNear(62_(1) * Degree));
  EXPECT_THAT(elements.mean_inclination_interval().max,
              IsNear(77_(1) * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().min,
              IsNear(51_(1) * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().max,
              IsNear(129_(1) * Degree));

  // The conservation of the “тривиального интеграла a = const” [Лид61, p. 25]
  // is excellent: while the sun is nudging and deforming the orbit, it is not
  // pumping energy into nor out of it.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().min,
              IsNear(14'910.0_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().max,
              IsNear(14'910.3_(1) * Kilo(Metre)));

  // The integral c₁ is preserved quite well: we have an exchange between
  // inclination and eccentricity.
  EXPECT_THAT(c₁.min, IsNear(0.042_(1)));
  EXPECT_THAT(c₁.max, IsNear(0.050_(1)));

  // The integral c₂ is also conserved: the long-term evolution of the orbital
  // elements is as described in [Лид61].
  EXPECT_THAT(c₂.min, IsNear(-0.091_(1)));
  EXPECT_THAT(c₂.max, IsNear(-0.083_(1)));

  // TODO(egg): The above are integrals of motion only when averaging over an
  // orbit of the perturbing body (so here, over the orbit of Mercury); see what
  // things look like under a moving average.
}
#endif

}  // namespace astronomy
}  // namespace principia
