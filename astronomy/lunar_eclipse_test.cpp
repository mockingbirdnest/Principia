#include <memory>

#include "astronomy/frames.hpp"
#include "astronomy/time_scales.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/sign.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "numerics/root_finders.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {

using ::testing::Eq;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_time_scales;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_sign;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_astronomy;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics;

namespace {

Time const bisection_interval = 10 * Minute;

Sign const U14 = Sign::Positive();
Sign const U23 = Sign::Negative();

char const arrow[] = "-------------------> ";

}  // namespace

class LunarEclipseTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<ICRS>::NewtonianMotionEquation>(),
            /*step=*/10 * Minute));
    r_sun_ = solar_system_1950_.mean_radius("Sun");
    r_moon_ = solar_system_1950_.mean_radius("Moon");

    // We follow A. Danjon, see http://eclipse.gsfc.nasa.gov/LEcat5/shadow.html
    // and http://web.archive.org/web/20170915031958/http://www.imcce.fr:80/langues/fr/ephemerides/phenomenes/eclipses/lune/.
    // In particular, this means that we must use the equatorial radius for the
    // Earth, not the mean radius.
    r_earth_ = TerrestrialEquatorialRadius;
    atmospheric_depth_ = (1.0 / 85.0 - 1.0 / 594.0) * r_earth_;
  }

  // A positive `time_error` means that the actual contact happens after
  // `current_time`.
  void CheckLunarUmbralEclipse(Instant const& current_time,
                               Sign const moon_offset_sign,
                               ApproximateQuantity<Angle> const& angular_error,
                               ApproximateQuantity<Time> const& time_error) {
    EXPECT_OK(ephemeris_->Prolong(current_time + bisection_interval));
    auto const sun = solar_system_1950_.massive_body(*ephemeris_, "Sun");
    auto const earth = solar_system_1950_.massive_body(*ephemeris_, "Earth");
    auto const moon = solar_system_1950_.massive_body(*ephemeris_, "Moon");

    // Angle formed by a right circular cone with sides defined by tangent lines
    // between Sun and Earth, and axis running through the centres of each.
    auto const umbral_half_aperture = [earth, sun](
        Instant const& t) {
      auto const q_sun = ephemeris_->trajectory(sun)->EvaluatePosition(t);
      auto const q_earth = ephemeris_->trajectory(earth)->EvaluatePosition(t);
      return ArcSin((r_sun_ - (r_earth_ + atmospheric_depth_)) /
                    (q_sun - q_earth).Norm());
    };

    auto const earth_moon_angle = [earth,
                                   moon,
                                   moon_offset_sign,
                                   sun,
                                   umbral_half_aperture](Instant const& t) {
      auto const q_sun = ephemeris_->trajectory(sun)->EvaluatePosition(t);
      auto const q_moon = ephemeris_->trajectory(moon)->EvaluatePosition(t);
      auto const q_earth = ephemeris_->trajectory(earth)->EvaluatePosition(t);
      auto const apex_of_moon_locus_at_umbral_contact =
          q_earth +
          Normalize(q_earth - q_sun) *
              (r_earth_ + atmospheric_depth_ + moon_offset_sign * r_moon_) /
              Sin(umbral_half_aperture(t));
      // Angle between Earth and Moon as seen at
      // `apex_of_moon_locus_at_umbral_contact`.
      return AngleBetween(apex_of_moon_locus_at_umbral_contact - q_earth,
                          apex_of_moon_locus_at_umbral_contact - q_moon);
    };

    // We are at the desired contact if the angle between Earth and Moon from
    // the apex of locus of the moon at that contact is the same value as the
    // half-aperture of the umbra (Earth-Sun cone).
    EXPECT_THAT(AbsoluteError(umbral_half_aperture(current_time),
                              earth_moon_angle(current_time)),
                IsNear(angular_error))
        << NAMED(umbral_half_aperture(current_time)) << ", "
        << NAMED(earth_moon_angle(current_time)) << ", " << NAMED(current_time);

    Instant const& actual_contact_time = Bisect(
        [earth_moon_angle, umbral_half_aperture](Instant const& t) {
          return umbral_half_aperture(t) - earth_moon_angle(t);
        },
        current_time - bisection_interval,
        current_time + bisection_interval);
    EXPECT_THAT(actual_contact_time - current_time, IsNear(time_error))
        << NAMED(actual_contact_time) << ", " << NAMED(current_time);
    LOG(INFO) << arrow << AbsoluteError(umbral_half_aperture(current_time),
                                        earth_moon_angle(current_time))
              << " " << actual_contact_time - current_time;
  }

  // A positive `time_error` means that the actual contact happens after
  // `current_time`.
  void CheckLunarPenumbralEclipse(
      Instant const& current_time,
      Sign const moon_offset_sign,
      ApproximateQuantity<Angle> const& angular_error,
      ApproximateQuantity<Time> const& time_error) {
    EXPECT_OK(ephemeris_->Prolong(current_time + bisection_interval));
    auto const sun = solar_system_1950_.massive_body(*ephemeris_, "Sun");
    auto const earth = solar_system_1950_.massive_body(*ephemeris_, "Earth");
    auto const moon = solar_system_1950_.massive_body(*ephemeris_, "Moon");

    auto const penumbral_half_aperture = [earth, sun](
        Instant const& t) {
      auto const q_sun = ephemeris_->trajectory(sun)->EvaluatePosition(t);
      auto const q_earth = ephemeris_->trajectory(earth)->EvaluatePosition(t);
      return ArcSin((r_sun_ + r_earth_ + atmospheric_depth_) /
                    (q_sun - q_earth).Norm());
    };

    auto const earth_moon_angle = [earth,
                                   moon,
                                   moon_offset_sign,
                                   sun,
                                   penumbral_half_aperture](Instant const& t) {
      auto const q_sun = ephemeris_->trajectory(sun)->EvaluatePosition(t);
      auto const q_moon = ephemeris_->trajectory(moon)->EvaluatePosition(t);
      auto const q_earth = ephemeris_->trajectory(earth)->EvaluatePosition(t);

      auto const apex_of_moon_locus_at_penumbral_contact =
          q_earth +
          Normalize(q_sun - q_earth) *
              (r_earth_ + atmospheric_depth_ + moon_offset_sign * r_moon_) /
              Sin(penumbral_half_aperture(t));
      // Angle between Earth and Moon as seen at
      // apex_of_moon_locus_at_penumbral_contact.
      return AngleBetween(apex_of_moon_locus_at_penumbral_contact - q_earth,
                          apex_of_moon_locus_at_penumbral_contact - q_moon);
    };

    // We are at the desired contact if the angle between Earth and Moon from
    // the apex of locus of the moon at that contact is the same value as the
    // half-aperture of the penumbra.
    EXPECT_THAT(AbsoluteError(penumbral_half_aperture(current_time),
                              earth_moon_angle(current_time)),
                IsNear(angular_error))
        << NAMED(penumbral_half_aperture(current_time)) << ", "
        << NAMED(earth_moon_angle(current_time)) << ", " << NAMED(current_time);

    Instant const& actual_contact_time = Bisect(
        [earth_moon_angle, penumbral_half_aperture](Instant const& t) {
          return penumbral_half_aperture(t) - earth_moon_angle(t);
        },
        current_time - bisection_interval,
        current_time + bisection_interval);
    EXPECT_THAT(actual_contact_time - current_time, IsNear(time_error))
        << NAMED(actual_contact_time) << ", " << NAMED(current_time);
    LOG(INFO) << arrow
              << Sign(penumbral_half_aperture(current_time) -
                      earth_moon_angle(current_time)) *
                     AbsoluteError(penumbral_half_aperture(current_time),
                                   earth_moon_angle(current_time))
              << " " << actual_contact_time - current_time;
  }

  static SolarSystem<ICRS> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
  static Length r_sun_;
  static Length r_earth_;
  static Length r_moon_;
  static Length atmospheric_depth_;
};

SolarSystem<ICRS> LunarEclipseTest::solar_system_1950_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2433282_500000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> LunarEclipseTest::ephemeris_;
Length LunarEclipseTest::r_sun_;
Length LunarEclipseTest::r_earth_;
Length LunarEclipseTest::r_moon_;
Length LunarEclipseTest::atmospheric_depth_;

#if !defined(_DEBUG)

TEST_F(LunarEclipseTest, Year1950) {
  {
    constexpr auto P1 = "1950-04-02T18:10:49"_UT1;
    constexpr auto U1 = "1950-04-02T19:09:19"_UT1;
    constexpr auto U2 = "1950-04-02T20:30:38"_UT1;
    constexpr auto U3 = "1950-04-02T20:57:33"_UT1;
    constexpr auto U4 = "1950-04-02T22:18:54"_UT1;
    constexpr auto P4 = "1950-04-02T23:17:21"_UT1;

    CheckLunarPenumbralEclipse(P1, U14, 1.4e-5_(1) * Radian, 31_(1) * Second);
    CheckLunarUmbralEclipse(U1, U14,    2.1e-5_(1) * Radian, 33_(1) * Second);
    CheckLunarUmbralEclipse(U2, U23,    1.8e-5_(1) * Radian, 42_(1) * Second);
    CheckLunarUmbralEclipse(U3, U23,    9.6e-6_(1) * Radian, 24_(1) * Second);
    CheckLunarUmbralEclipse(U4, U14,    2.0e-5_(1) * Radian, 31_(1) * Second);
    CheckLunarPenumbralEclipse(P4, U14, 1.5e-5_(1) * Radian, 33_(1) * Second);
  }

  {
    constexpr auto P1 = "1950-09-26T01:21:43"_UT1;
    constexpr auto U1 = "1950-09-26T02:31:48"_UT1;
    constexpr auto U2 = "1950-09-26T03:54:33"_UT1;
    constexpr auto U3 = "1950-09-26T04:38:49"_UT1;
    constexpr auto U4 = "1950-09-26T06:01:33"_UT1;
    constexpr auto P4 = "1950-09-26T07:11:47"_UT1;

    CheckLunarPenumbralEclipse(P1, U14, 1.5e-5_(1) * Radian, 37_(1) * Second);
    CheckLunarUmbralEclipse(U1, U14,    2.3e-5_(1) * Radian, 39_(1) * Second);
    CheckLunarUmbralEclipse(U2, U23,    2.8e-5_(1) * Radian, 45_(1) * Second);
    CheckLunarUmbralEclipse(U3, U23,    1.9e-5_(1) * Radian, 32_(1) * Second);
    CheckLunarUmbralEclipse(U4, U14,    2.2e-5_(1) * Radian, 37_(1) * Second);
    CheckLunarPenumbralEclipse(P4, U14, 1.6e-5_(1) * Radian, 38_(1) * Second);
  }
}

TEST_F(LunarEclipseTest, Year1951) {
  {
    constexpr auto P1 = "1951-03-23T08:50:00"_UT1;
    constexpr auto P4 = "1951-03-23T12:24:19"_UT1;

    CheckLunarPenumbralEclipse(P1, U14, 9.5e-6_(1) * Radian, 33_(1) * Second);
    CheckLunarPenumbralEclipse(P4, U14, 7.8e-6_(1) * Radian, 27_(1) * Second);
  }

  {
    constexpr auto P1 = "1951-09-15T10:29:16"_UT1;
    constexpr auto P4 = "1951-09-15T14:23:52"_UT1;

    CheckLunarPenumbralEclipse(P1, U14, 9.5e-6_(1) * Radian, 30_(1) * Second);
    CheckLunarPenumbralEclipse(P4, U14, 8.1e-6_(1) * Radian, 26_(1) * Second);
  }
}

TEST_F(LunarEclipseTest, Year1952) {
  {
    constexpr auto P1 = "1952-02-10T22:08:20"_UT1;
    constexpr auto U1 = "1952-02-11T00:04:17"_UT1;
    constexpr auto U4 = "1952-02-11T01:14:24"_UT1;
    constexpr auto P4 = "1952-02-11T03:10:15"_UT1;

    CheckLunarPenumbralEclipse(P1, U14, 1.1e-5_(1) * Radian, 32_(1) * Second);
    CheckLunarUmbralEclipse(U1, U14,    4.2e-6_(1) * Radian, 21_(1) * Second);
    CheckLunarUmbralEclipse(U4, U14,    1.0e-5_(1) * Radian, 52_(1) * Second);
    CheckLunarPenumbralEclipse(P4, U14, 1.4e-5_(1) * Radian, 42_(1) * Second);
  }

  {
    constexpr auto P1 = "1952-08-05T17:28:13"_UT1;
    constexpr auto U1 = "1952-08-05T18:33:49"_UT1;
    constexpr auto U4 = "1952-08-05T21:01:00"_UT1;
    constexpr auto P4 = "1952-08-05T22:06:35"_UT1;

    CheckLunarPenumbralEclipse(P1, U14, 8.4e-6_(1) * Radian, 20_(1) * Second);
    CheckLunarUmbralEclipse(U1, U14,    9.6e-6_(1) * Radian, 20_(1) * Second);
    CheckLunarUmbralEclipse(U4, U14,    1.3e-5_(1) * Radian, 27_(1) * Second);
    CheckLunarPenumbralEclipse(P4, U14, 1.2e-5_(1) * Radian, 28_(1) * Second);
  }
}

TEST_F(LunarEclipseTest, DISABLED_Year2000) {
  constexpr auto P1 = "2000-01-21T02:04:26"_UT1;
  constexpr auto U1 = "2000-01-21T03:01:50"_UT1;
  constexpr auto U2 = "2000-01-21T04:05:01"_UT1;
  constexpr auto U3 = "2000-01-21T05:22:00"_UT1;
  constexpr auto U4 = "2000-01-21T06:25:09"_UT1;
  constexpr auto P4 = "2000-01-21T07:22:38"_UT1;

  CheckLunarPenumbralEclipse(P1, U14, 7.8e-5_(1) * Radian, -162_(1) * Second);
  CheckLunarUmbralEclipse(U1, U14,    1.1e-4_(1) * Radian, -161_(1) * Second);
  CheckLunarUmbralEclipse(U2, U23,    2.0e-4_(1) * Radian, -163_(1) * Second);
  CheckLunarUmbralEclipse(U3, U23,    2.0e-4_(1) * Radian, -157_(1) * Second);
  CheckLunarUmbralEclipse(U4, U14,    1.1e-4_(1) * Radian, -157_(1) * Second);
  CheckLunarPenumbralEclipse(P4, U14, 7.5e-5_(1) * Radian, -156_(1) * Second);
}

TEST_F(LunarEclipseTest, DISABLED_Year2048) {
  // No UT1 in the future, but's that's what NASA gives us (don't ask).  Using
  // TT plus the ΔT that they use.
  constexpr Time ΔT = 91.2 * Second;
  constexpr auto P1 = "2048-01-01T03:52:39"_TT + ΔT;
  constexpr auto U1 = "2048-01-01T05:05:17"_TT + ΔT;
  constexpr auto U2 = "2048-01-01T06:24:27"_TT + ΔT;
  constexpr auto U3 = "2048-01-01T07:20:23"_TT + ΔT;
  constexpr auto U4 = "2048-01-01T08:39:33"_TT + ΔT;
  constexpr auto P4 = "2048-01-01T09:52:05"_TT + ΔT;

  CheckLunarPenumbralEclipse(P1, U14, 1.5e-4_(1) * Radian, -358_(1) * Second);
  CheckLunarUmbralEclipse(U1, U14,    2.2e-4_(1) * Radian, -356_(1) * Second);
  CheckLunarUmbralEclipse(U2, U23,    2.7e-4_(1) * Radian, -356_(1) * Second);
  CheckLunarUmbralEclipse(U3, U23,    3.1e-4_(1) * Radian, -355_(1) * Second);
  CheckLunarUmbralEclipse(U4, U14,    2.2e-4_(1) * Radian, -355_(1) * Second);
  CheckLunarPenumbralEclipse(P4, U14, 1.5e-4_(1) * Radian, -353_(1) * Second);
}

#endif

}  // namespace astronomy
}  // namespace principia
