#define PRINCIPIA_LOG_TO_MATHEMATICA 0

#include <memory>

#include "astronomy/date_time.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/mercury_orbiter.hpp"
#include "astronomy/orbital_elements.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

#if PRINCIPIA_LOG_TO_MATHEMATICA
#include "mathematica/logger.hpp"
#endif

namespace principia {
namespace astronomy {

using ::testing::AnyOf;
using ::testing::Eq;
using namespace principia::astronomy::_date_time;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_mercury_orbiter;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;

// A test that showcases the eccentricity-inclination exchange mechanism
// described in [Лид61] and [Koz62].  We follow the treatment in [Лид61].
class Лидов古在Test : public ::testing::Test {
 protected:
  using MercuryCentredInertial =
      Frame<struct MercuryCentredInertialTag, NonRotating>;

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
  BodyCentredNonRotatingReferenceFrame<ICRS, MercuryCentredInertial>
      mercury_frame_;
};

#if !_DEBUG
TEST_F(Лидов古在Test, MercuryOrbiter) {
  DiscreteTrajectory<ICRS> icrs_trajectory;
  EXPECT_OK(icrs_trajectory.Append(
      MercuryOrbiterInitialTime, MercuryOrbiterInitialDegreesOfFreedom<ICRS>));
  auto& icrs_segment = icrs_trajectory.segments().front();
  // Carefully tuned.
  icrs_segment.SetDownsampling({.max_dense_intervals = 10'000,
                                .tolerance = 1 * Milli(Metre)});
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
#if PRINCIPIA_LOG_TO_MATHEMATICA
  mathematica::Logger logger(
      SOLUTION_DIR / "mathematica" /
          PRINCIPIA_UNICODE_PATH("лидов_古在.generated.wl"),
      /*make_unique=*/false);
#endif

  DiscreteTrajectory<MercuryCentredInertial> mercury_centred_trajectory;
  for (auto const& [t, dof] : icrs_trajectory) {
    EXPECT_OK(mercury_centred_trajectory.Append(
        t, mercury_frame_.ToThisFrameAtTime(t)(dof)));
#if PRINCIPIA_LOG_TO_MATHEMATICA
    logger.Append(
        "q",
        mercury_centred_trajectory.back().degrees_of_freedom.position(),
        mathematica::ExpressIn(Metre));
#endif
  }

  EXPECT_THAT(mercury_centred_trajectory.size(),
              AnyOf(Eq(1'534'438),    // Windows, Ubuntu.
                    Eq(1'534'680)));  // macOS.
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
#if PRINCIPIA_LOG_TO_MATHEMATICA
    logger.Append("t", elements.time, mathematica::ExpressIn(Second));
    logger.Append("a", elements.semimajor_axis, mathematica::ExpressIn(Metre));
    logger.Append("e", elements.eccentricity);
    logger.Append("i", elements.inclination, mathematica::ExpressIn(Radian));
    logger.Append(R"(\[Omega])",
                  elements.argument_of_periapsis,
                  mathematica::ExpressIn(Radian));
#endif
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
  // pumping energy into nor out of it.  The true values are 14'910.01 and
  // 14'910.28 km.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().min,
              AnyOf(IsNear(14'910.01_(1) * Kilo(Metre)),    // Windows, macOS.
                    IsNear(14'909.96_(1) * Kilo(Metre))));  // Ubuntu.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().max,
              AnyOf(IsNear(14'910.28_(1) * Kilo(Metre)),    // Windows, macOS.
                    IsNear(14'910.29_(1) * Kilo(Metre))));  // Ubuntu.

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
#undef PRINCIPIA_LOG_TO_MATHEMATICA
