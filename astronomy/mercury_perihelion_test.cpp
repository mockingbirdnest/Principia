
#include <numeric>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using base::not_null;
using geometry::AngleBetween;
using integrators::McLachlanAtela1992Order5Optimal;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MassiveBody;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::si::AstronomicalUnit;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using testing_utilities::AbsoluteError;

namespace astronomy {

class MercuryPerihelionTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
            McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
            /*step=*/45 * Minute));
  }

  MercuryPerihelionTest()
      : sun_(solar_system_1950_.massive_body(*ephemeris_, "Sun")),
        mercury_(solar_system_1950_.massive_body(*ephemeris_, "Mercury")) {}

  static SolarSystem<ICRFJ2000Equator> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  not_null<MassiveBody const*> sun_;
  not_null<MassiveBody const*> mercury_;
};

SolarSystem<ICRFJ2000Equator> MercuryPerihelionTest::solar_system_1950_;
std::unique_ptr<Ephemeris<ICRFJ2000Equator>> MercuryPerihelionTest::ephemeris_;

TEST_F(MercuryPerihelionTest, PrintPerihelion) {
  DiscreteTrajectory<ICRFJ2000Equator> sun_apoapsides;
  DiscreteTrajectory<ICRFJ2000Equator> sun_periapsides;
  DiscreteTrajectory<ICRFJ2000Equator> mercury_apoapsides;
  DiscreteTrajectory<ICRFJ2000Equator> mercury_periapsides;
  ephemeris_->Prolong(solar_system_1950_.epoch() + 10 * JulianYear);
  ephemeris_->ComputeApsides(sun_,
                             mercury_,
                             sun_apoapsides,
                             sun_periapsides,
                             mercury_apoapsides,
                             mercury_periapsides);

  std::experimental::optional<Instant> previous_time;
  std::experimental::optional<Displacement<ICRFJ2000Equator>>
      previous_displacement;
  std::vector<AngularFrequency> precessions;
  for (auto sun_it = sun_periapsides.Begin(),
            mercury_it = mercury_periapsides.Begin();
       sun_it != sun_periapsides.End() &&
       mercury_it != mercury_periapsides.End();
       ++sun_it, ++mercury_it) {
    Instant const time = sun_it.time();
    Displacement<ICRFJ2000Equator> const displacement =
        sun_it.degrees_of_freedom().position() -
        mercury_it.degrees_of_freedom().position();
    //LOG(ERROR)<<time<<": "<<displacement.Norm();
    //EXPECT_LT(AbsoluteError(displacement.Norm(),
    //                        3.075030670219868e-01 * AstronomicalUnit),
    //          900 * Kilo(Metre));
    if (previous_time) {
      AngularFrequency const precession =
          AngleBetween(displacement, *previous_displacement) /
          (time - *previous_time);
      precessions.push_back(precession);
      LOG(ERROR)<<precession * 100 * JulianYear / (1 * ArcSecond);
      //EXPECT_LT(
      //    AbsoluteError(time - *previous_time, 8.796888204428582e+01 * Day),
      //    75 * Second);
    }
    previous_time = time;
    previous_displacement = displacement;
  }
  AngularFrequency average;
  for (auto const precession : precessions) {
    average += precession;
  }
  average /= precessions.size();
  LOG(ERROR) << "Average: " << average * 100 * JulianYear / (1 * ArcSecond);
}

}  // namespace astronomy
}  // namespace principia
