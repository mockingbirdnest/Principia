
#include <string>
#include <vector>
#include <utility>

#include "astronomy/orbit_ground_track.hpp"
#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "numerics/polynomial.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace astronomy {

using base::make_not_null_unique;
using base::not_null;
using base::OFStream;
using geometry::Instant;
using geometry::Interval;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::QuinlanTremaine1990Order12;
using numerics::EstrinEvaluator;
using numerics::PolynomialInMonomialBasis;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::BodySurfaceDynamicFrame;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::Angle;
using quantities::Mod;
using quantities::Pow;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::astronomy::TerrestrialEquatorialRadius;
using quantities::si::ArcMinute;
using quantities::si::ArcSecond;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::DifferenceFrom;
using testing_utilities::IsNear;
using testing_utilities::IsOk;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_⑴;
using ::testing::AllOf;
using ::testing::Field;
using ::testing::Lt;
using ::testing::Property;

namespace {

struct SP3Files {
  std::vector<std::string> names;
  StandardProduct3::Dialect dialect;

  static SP3Files const& GNSS();
  static SP3Files const& SPOT5();
  static SP3Files const& Sentinel3A();
  static SP3Files const& TOPEXPoséidon();
};

struct SP3Orbit {
  StandardProduct3::SatelliteIdentifier satellite;
  SP3Files files;
};

SP3Files const& SP3Files::GNSS() {
  static const SP3Files files = {{"WUM0MGXFIN_20190970000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20190980000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20190990000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191000000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191010000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191020000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191030000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191040000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191050000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191060000_01D_15M_ORB.SP3"},
                                 StandardProduct3::Dialect::ChineseMGEX};
  return files;
}

SP3Files const& SP3Files::SPOT5() {
  static const SP3Files files = {{"ssasp501.b10170.e10181.D__.sp3"},
                                 StandardProduct3::Dialect::Standard};
  return files;
}

SP3Files const& SP3Files::Sentinel3A() {
  static const SP3Files files = {{"ssas3a20.b18358.e19003.DG_.sp3"},
                                 StandardProduct3::Dialect::Standard};
  return files;
}

SP3Files const& SP3Files::TOPEXPoséidon() {
  static const SP3Files files = {{"grgtop03.b97344.e97348.D_S.sp3"},
                                 StandardProduct3::Dialect::GRGS};
  return files;
}

}  // namespace

class OrbitAnalysisTest : public ::testing::Test {
 protected:
  OrbitAnalysisTest()
      : earth_1957_(RemoveAllButEarth(SolarSystem<ICRS>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2436116_311504629.proto.txt"))),
        // As there is only one body left in this ephemeris, integration is
        // exact up to roundoff error regardless of the step size.  Use a long
        // step so that we do not waste time computing the ephemeris (it starts
        // in 1957, and we need it until 2019).
        ephemeris_(earth_1957_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
                /*step=*/1 * JulianYear))),
        earth_(*earth_1957_.rotating_body(*ephemeris_, "Earth")) {}

  // Returns a GCRS trajectory obtained by stitching together the trajectories
  // of |sp3_orbit.satellites| in |sp3_orbit.files|.
  not_null<std::unique_ptr<DiscreteTrajectory<GCRS>>> EarthCentredTrajectory(
      SP3Orbit const& sp3_orbit) {
    BodyCentredNonRotatingDynamicFrame<ICRS, GCRS> gcrs{ephemeris_.get(),
                                                        &earth_};
    BodySurfaceDynamicFrame<ICRS, ITRS> itrs{ephemeris_.get(), &earth_};

    auto result = make_not_null_unique<DiscreteTrajectory<GCRS>>();
    for (auto const& file : sp3_orbit.files.names) {
      StandardProduct3 sp3(
          SOLUTION_DIR / "astronomy" / "standard_product_3" / file,
          sp3_orbit.files.dialect);
      std::vector<not_null<DiscreteTrajectory<ITRS> const*>> const& orbit =
          sp3.orbit(sp3_orbit.satellite);
      CHECK_EQ(orbit.size(), 1);
      auto const& arc = *orbit.front();
      for (auto it = arc.Begin(); it != arc.End(); ++it) {
        ephemeris_->Prolong(it.time());
        result->Append(
            it.time(),
            gcrs.ToThisFrameAtTime(it.time())(
                itrs.FromThisFrameAtTime(it.time())(it.degrees_of_freedom())));
      }
    }
    return result;
  }

  std::tuple<OrbitalElements, OrbitRecurrence, OrbitGroundTrack>
  ElementsAndRecurrence(SP3Orbit const& orbit) {
    auto const earth_centred_trajectory = EarthCentredTrajectory(orbit);
    auto const elements = OrbitalElements::ForTrajectory(
                              *earth_centred_trajectory,
                              earth_,
                              MasslessBody{}).ValueOrDie();

    {
      auto const identifier = (std::stringstream() << orbit.satellite).str();
      OFStream f(SOLUTION_DIR / "mathematica" /
                 (identifier + "_elements.generated.wl"));
      f << mathematica::Assign(identifier + "osculatingEquinoctialElements",
                               elements.osculating_equinoctial_elements());
      f << mathematica::Assign(identifier + "meanEquinoctialElements",
                               elements.mean_equinoctial_elements());
    }

    auto const recurrence =
        OrbitRecurrence::ClosestRecurrence(elements.nodal_period(),
                                           elements.nodal_precession(),
                                           earth_,
                                           /*max_abs_Cᴛₒ=*/100);
    // Since our ITRS-to-ICRS conversion disregards the precession of the
    // equinoxes, it is not completely clear which year we should be dealing
    // with here.  Given that the SP3 files have data in the ITRS (whose equator
    // is the equator of date), and that we map that to the ICRS equator, our
    // IRCS here is more like an equator-of-date frame (although it is not clear
    // whether its x axis resembles the equinox of the date), so we pick
    // the tropical year.
    // We use values based Newcomb's formula for the mean longitude of the sun,
    //   L = 279°41′48″.04 + 129 602 768″.13 T + 1.089″ T²,
    // where T is in Julian centuries since 1900, January 0, Greenwich Mean noon.
    // Ephemeris time (ET) was defined by the IAU (10th general assembly (1958),
    // commissions 4 and 31, recommendation 2) from Newcomb's tables, with
    // 1900, January 0 at 12 ET being the instant at which the mean longitude of
    // the sun was 279°41′48″.08, and the ET second being 1/31 556 925.9747 of
    // the tropical year at that epoch—where 31 556 925.9747 s =
    // 2π rad / (129 602 768″.13 / 100 a).  TDT, and later TT, were then defined
    // in such a way as to achieve approximate continuity with ET, see 16th
    // general assembly (1976), commission 4, recommendation 5, note 2, and 21st
    // general assembly (1991), resolution A4, recommendation IV, note 4.
    // We can therefore work with this formula in TT.
    PolynomialInMonomialBasis<Angle, Instant, 2, EstrinEvaluator> const
        newcomb_mean_longitude(
            {279 * Degree + 41 * ArcMinute + 48.04 * ArcSecond,
             129'602'768.13 * ArcSecond / (100 * JulianYear),
             1.089 * ArcSecond / Pow<2>(100 * JulianYear)},
            "1899-12-31T12:00:00"_TT);
    auto const ground_track = OrbitGroundTrack::ForTrajectory(
        *earth_centred_trajectory,
        earth_,
        recurrence,
        {{/*epoch=*/J2000,
          /*mean_longitude_at_epoch=*/newcomb_mean_longitude.Evaluate(J2000),
          /*year*/ 2 * π * Radian /
              newcomb_mean_longitude.EvaluateDerivative(J2000)}});
    return {elements, recurrence, ground_track};
  }

  SolarSystem<ICRS> earth_1957_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> ephemeris_;
  RotatingBody<ICRS> const& earth_;

 private:
  SolarSystem<ICRS> RemoveAllButEarth(SolarSystem<ICRS> solar_system) {
    std::vector<std::string> const names = solar_system.names();
    for (auto const& name : names) {
      if (name != "Earth") {
        solar_system.RemoveMassiveBody(name);
      }
    }
    return solar_system;
  }
};

#if !defined(_DEBUG)

// COSPAR ID 2010-001A, SVN C003.
// 北斗二號 GEO01.
// PRN C01, GEO, 140.0° E.
TEST_F(OrbitAnalysisTest, 北斗GEO) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::北斗, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'166_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(1.42_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.000186_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(166_⑴ * Degree));
}

// COSPAR ID 2010-036A, SVN C005.
// 北斗二號 IGSO01.
// PRN C06, IGSO, 117°E.
TEST_F(OrbitAnalysisTest, 北斗IGSO) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::北斗, 6}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'161_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(54.19_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0078_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(232_⑴ * Degree));
}

// COSPAR ID 2010-045A, SVN J001.
// Block I-Q, みちびき初号機.
// PRN J01, quasi-zenith orbit.
TEST_F(OrbitAnalysisTest, みちびきQZO) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::みちびき, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  // Expected orbital elements from the Quasi-Zenith Satellite System
  // Performance Standard (PS-QZSS-001).
  EXPECT_THAT(
      elements.mean_semimajor_axis_interval().midpoint(),
      AbsoluteErrorFrom(42'165 * Kilo(Metre), IsNear(6.3_⑴ * Kilo(Metre))));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(41_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.075_⑴));
  // The operational range is 270° ± 2.5°.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(270_⑴ * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(0.12_⑴ * Degree));
}

// COSPAR ID 2017-048A, SVN J003.
// Block II-G, みちびき3号機.
// PRN J07, GEO.
TEST_F(OrbitAnalysisTest, みちびきGEO) {
  auto j07_files = SP3Files::GNSS();
  // J07 is missing from the last two files.
  j07_files.names.resize(j07_files.names.size() - 2);
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::みちびき, 7}, j07_files});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'166_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(0.067_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.00023_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(224_⑴ * Degree));
}

// COSPAR ID 2018-078B, SVN C216.
// 北斗三號 MEO15 (Shanghai Engineering Center for Microsatellites).
// PRN C34, slot A-7.
TEST_F(OrbitAnalysisTest, 北斗MEO) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::北斗, 34}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, -1),
                    Property(&OrbitRecurrence::Cᴛₒ, 7)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(27'906_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(55.10_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.000558_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(1.003_⑴ * Degree));
}

// COSPAR ID 2016-030A.
// Galileo-Full Operational Capability Flight Model 10 (GSAT0210) “Danielė”.
// PRN E01, slot A02.
TEST_F(OrbitAnalysisTest, GalileoNominalSlot) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::Galileo, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, -3),
                    Property(&OrbitRecurrence::Cᴛₒ, 10)));

  // Reference elements from
  // https://www.gsc-europa.eu/system-status/orbital-and-technical-parameters.
  Instant const reference_epoch = "2016-11-21T00:00:00"_UTC;
  Instant const initial_time = elements.mean_elements().front().time;
  Instant const mean_time =
      initial_time + (elements.mean_elements().back().time - initial_time) / 2;

  auto const nominal_nodal_precession = -0.02764398 * Degree / Day;
  auto const nominal_anomalistic_mean_motion =
      613.72253566 * Degree / Day;

  EXPECT_THAT(
      elements.nodal_precession(),
      AllOf(AbsoluteErrorFrom(nominal_nodal_precession,
                              IsNear(0.00032_⑴ * Degree / Day)),
            RelativeErrorFrom(nominal_nodal_precession, IsNear(0.011_⑴))));
  EXPECT_THAT(2 * π * Radian / elements.anomalistic_period(),
              AllOf(AbsoluteErrorFrom(nominal_anomalistic_mean_motion,
                                      IsNear(0.53_⑴ * Degree / Day)),
                    RelativeErrorFrom(nominal_anomalistic_mean_motion,
                                      IsNear(0.00086_⑴))));

  EXPECT_THAT(
      elements.mean_semimajor_axis_interval().midpoint(),
      AbsoluteErrorFrom(29'599.8 * Kilo(Metre), IsNear(0.35_⑴ * Kilo(Metre))));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(00'000.084_⑴ * Kilo(Metre)));

  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.000'17_⑴));  // Nominal: 0.0.
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(0.000'015_⑴));

  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              AbsoluteErrorFrom(56.0 * Degree, IsNear(0.61_⑴ * Degree)));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(00.01_⑴ * Degree));

  EXPECT_THAT(
      Mod(elements.mean_longitude_of_ascending_node_interval().midpoint() -
              nominal_nodal_precession * (mean_time - reference_epoch),
          2 * π * Radian),
      AbsoluteErrorFrom(317.632 * Degree, IsNear(0.082_⑴ * Degree)));

  // Note that the reference parameters have e = 0, and therefore conventionally
  // set ω = 0, ω′ = 0.
  // However, e is never quite 0; we can compute a mean ω.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(88_⑴ * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(6.3_⑴ * Degree));

  // Since the reference parameters conventionally set ω = 0, the given mean
  // anomaly is actually the mean argument of latitude; in order to get numbers
  // consistent with theirs, we must look at ω + M.
  EXPECT_THAT(Mod(elements.mean_elements().front().argument_of_periapsis +
                      elements.mean_elements().front().mean_anomaly -
                      nominal_anomalistic_mean_motion *
                          (initial_time - reference_epoch),
                  2 * π * Radian),
              AbsoluteErrorFrom(225.153 * Degree, IsNear(0.53_⑴ * Degree)));
}

// COSPAR ID 2014-050B, SVN E202
// Galileo-Full Operational Capability Flight Model 2 (GSAT0202) “Milena”.
// PRN E14, slot Ext02.
TEST_F(OrbitAnalysisTest, GalileoExtendedSlot) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::Galileo, 14}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, -3),
                    Property(&OrbitRecurrence::Cᴛₒ, 20)));

  // Reference elements from
  // https://www.gsc-europa.eu/system-status/orbital-and-technical-parameters.
  Instant const reference_epoch = "2016-11-21T00:00:00"_UTC;
  Instant const initial_time = elements.mean_elements().front().time;
  Instant const mean_time =
      initial_time + (elements.mean_elements().back().time - initial_time) / 2;

  auto const nominal_nodal_precession = -0.03986760 * Degree / Day;
  auto const nominal_apsidal_precession = 0.03383184 * Degree / Day;
  auto const nominal_anomalistic_mean_motion = 667.86467481 * Degree / Day;

  EXPECT_THAT(
      elements.nodal_precession(),
      AllOf(AbsoluteErrorFrom(nominal_nodal_precession,
                              IsNear(0.00023_⑴ * Degree / Day)),
            RelativeErrorFrom(nominal_nodal_precession, IsNear(0.0059_⑴))));
  EXPECT_THAT(2 * π * Radian / elements.anomalistic_period(),
              AllOf(AbsoluteErrorFrom(nominal_anomalistic_mean_motion,
                                      IsNear(0.00047_⑴ * Degree / Day)),
                    RelativeErrorFrom(nominal_anomalistic_mean_motion,
                                      IsNear(7.1e-07_⑴))));

  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              AbsoluteErrorFrom(27'977.6 * Kilo(Metre),
                                IsNear(0.0997_⑴ * Kilo(Metre))));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(00'000.096_⑴ * Kilo(Metre)));

  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              AbsoluteErrorFrom(0.162, IsNear(0.0041_⑴)));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(0.000'15_⑴));

  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              AbsoluteErrorFrom(49.850 * Degree, IsNear(0.77_⑴ * Degree)));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(00.0044_⑴ * Degree));

  EXPECT_THAT(
      Mod(elements.mean_longitude_of_ascending_node_interval().midpoint() -
              nominal_nodal_precession * (mean_time - reference_epoch),
          2 * π * Radian),
      AbsoluteErrorFrom(52.521 * Degree, IsNear(0.29_⑴ * Degree)));
  EXPECT_THAT(
      Mod(elements.mean_argument_of_periapsis_interval().midpoint() -
              nominal_apsidal_precession * (mean_time - reference_epoch),
          2 * π * Radian),
      AbsoluteErrorFrom(56.198 * Degree, IsNear(0.48_⑴ * Degree)));

  EXPECT_THAT(Mod(elements.mean_elements().front().mean_anomaly -
                      nominal_anomalistic_mean_motion *
                          (initial_time - reference_epoch),
                  2 * π * Radian),
              AbsoluteErrorFrom(136.069 * Degree, IsNear(2.5_⑴ * Degree)));
}

// COSPAR ID 2009-070A, SVN R730.
// ГЛОНАСС-М Космос 2456, Ураган-М № 730.
// PRN R01, plane 1.
TEST_F(OrbitAnalysisTest, ГЛОНАСС) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::ГЛОНАСС, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, 1),
                    Property(&OrbitRecurrence::Cᴛₒ, 8)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(25'507_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(64.20_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.00040_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(330_⑴ * Degree));
}

// COSPAR ID 2011-036A, SVN G063.
// GPS block IIF satellite.
// PRN G01, plane D, slot 2.
TEST_F(OrbitAnalysisTest, GPS) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::GPS, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(26'560_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(55.86_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0086_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(39_⑴ * Degree));
}

// COSPAR ID 1992-052A, TOPEX/Poséidon.
TEST_F(OrbitAnalysisTest, TOPEXPoséidon) {
  // The references for these orbital characteristics are:
  // — Bhat et al. (1998), TOPEX/Poseidon orbit maintenance for the first five
  //   years;
  // — Benada (1997), PO.DAAC Merged GDR (TOPEX-Poseidon) Generation B User’s
  //   handbook, version 2.0;
  // — Blanc et al. (1996), AVISO/Altimetry, AVISO User Handbook for Merged
  //   TOPEX/POSEIDON products, AVI-NT-02-101, Edition 3.0.

  auto [elements, recurrence, ground_track] =  // NOLINT(whitespace/braces)
      ElementsAndRecurrence({{StandardProduct3::SatelliteGroup::General, 1},
                             SP3Files::TOPEXPoséidon()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 13),
                    Property(&OrbitRecurrence::Dᴛₒ, -3),
                    Property(&OrbitRecurrence::Cᴛₒ, 10)));

  // Reference semimajor axis from the legend of figure 7 of Bhat et al.; that
  // value is given as 7 714.42942 in table 1 of Bhat et al., 7 714.4278 km in
  // Blanc et al., and 7 714.43 km in Benada.
  // Figure 7 of Bhat et al. shows that the mean semimajor axis varies by up to
  // 6 m either side of the reference value.  Our test data is from cycle 198;
  // reading the graph around that time shows that the mean semimajor axis was a
  // bit less than 3 m above the nominal value around that time.
  EXPECT_THAT(
      elements.mean_semimajor_axis_interval().midpoint(),
              DifferenceFrom(7714.42938 * Kilo(Metre), IsNear(2.93_⑴ * Metre)));
  // Reference inclination from the legend of figure 9 of Bhat et al.; that
  // value is given as 66.040° in table 1 of Bhat et al., 66.039° in Blanc et
  // al., and 66.04° in Benada.
  // Figure 9 of Bhat et al. shows that the observed values of the mean
  // inclination vary between 66.0375° and 66.0455° over the course of the
  // mission (66.0408° ± 0.0040° according to the abstract).  We can see from
  // the graph that during the period under test, in early December 1997, the
  // inclination is below the reference value, and varies by several thousandths
  // of a degree in a matter of days; the minimum and maximum that we compute
  // below therefore seem plausible.
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              DifferenceFrom(66.0408 * Degree, IsNear(-0.0025_⑴ * Degree)));
  EXPECT_THAT(elements.mean_inclination_interval(),
              AllOf(Field(&Interval<Angle>::min, IsNear(66.0365_⑴ * Degree)),
                    Field(&Interval<Angle>::max, IsNear(66.0401_⑴ * Degree))));

  // The same nominal values are given by Blanc et al., Benada, and Bhat et al:
  // e = 0.000095, ω = 90.0°.
  // However, these values are only meaningful if we take mean elements that are
  // free of variations over one 127-revolution ground track cycle (rather than
  // simply over one revolution).  Figure 8 of Bhat et al. shows that both the
  // theoretical and observed mean e and ω vary between 40 ppm and 140 ppm, and
  // between 60° and 120°, respectively.
  EXPECT_THAT(elements.mean_eccentricity_interval(),
              AllOf(Field(&Interval<double>::min, IsNear(88e-6_⑴)),
                    Field(&Interval<double>::max, IsNear(110e-6_⑴))));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval(),
              AllOf(Field(&Interval<Angle>::min, IsNear(74.7_⑴ * Degree)),
                    Field(&Interval<Angle>::max, IsNear(98.5_⑴ * Degree))));

  // Nominal longitude of the equatorial crossing of the first ascending pass
  // East of the ITRF zero-meridian (pass 135), as given in section 2 of Benada.
  // Blanc et al. round these longitudes to a hundredth of a degree, thus 0.71°
  // for pass 1.
  // We can see from figure 6 of Bhat et al. that, during the period under test,
  // the equatorial crossing is about 600 m east of the reference.
  EXPECT_THAT(
      ground_track.reduced_longitude_of_equator_crossing()->midpoint(),
      AbsoluteErrorFrom(0.7117 * Degree,
                        AllOf(IsNear(0.0051_⑴ * Degree),
                              IsNear(573_⑴ * Metre *
                                     (Radian / TerrestrialEquatorialRadius)))));

  // Nominal longitude of the equatorial crossing of pass 1, as given in the
  // auxiliary data table in Blanc et al.  The reference grid there lists that
  // longitude as 99.92°, and the table of equator crossing longitudes in Benada
  // lists it as 99.9249°.  However, the auxiliary data table in Benada gives a
  // longitude of 99.947° for pass 1, which looks like a typo.
  EXPECT_THAT(
      Mod(ground_track.reduced_longitude_of_equator_crossing()->midpoint() -
              ((135 - 1) / 2) * recurrence.equatorial_shift(),
          2 * π * Radian),
      AbsoluteErrorFrom(99.9242 * Degree, IsNear(0.0052_⑴ * Degree)));

  // Variability over the period under test (3.5 days).
  EXPECT_THAT(ground_track.reduced_longitude_of_equator_crossing()->measure(),
              IsNear(0.0025_⑴ * Degree));
  EXPECT_THAT(ground_track.reduced_longitude_of_equator_crossing()->measure() *
                  TerrestrialEquatorialRadius / Radian,
              IsNear(273_⑴ * Metre));

  // TOPEX/Poséidon is not sun-synchronous.
  EXPECT_THAT(ground_track.mean_solar_time_of_ascending_node()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(42_⑴ * Minute));
  EXPECT_THAT(ground_track.mean_solar_time_of_descending_node()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(42_⑴ * Minute));
}

// The following two satellites are sun-synchronous.

// COSPAR ID 2002-021A, SPOT-5 (Satellite Pour l’Observation de la Terre).
TEST_F(OrbitAnalysisTest, SPOT5) {
  auto [elements, recurrence, ground_track] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::General, 94}, SP3Files::SPOT5()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 5),
                    Property(&OrbitRecurrence::Cᴛₒ, 26)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(7'200_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(98.73_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0012_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(89.38_⑴ * Degree));

  // The nominal mean solar times of the nodes are 22:30 ascending, 10:30
  // descending.
  EXPECT_THAT(ground_track.mean_solar_time_of_ascending_node()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(22.452_⑴ * Hour));
  EXPECT_THAT(ground_track.mean_solar_time_of_descending_node()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(10.452_⑴ * Hour));
  EXPECT_THAT(ground_track.mean_solar_time_of_ascending_node()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(3.91_⑴ * Second));
  EXPECT_THAT(ground_track.mean_solar_time_of_descending_node()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(3.91_⑴ * Second));
}

// COSPAR ID 2016-011A, Sentinel-3A.
TEST_F(OrbitAnalysisTest, Sentinel3A) {
  auto [elements, recurrence, ground_track] =  // NOLINT(whitespace/braces)
      ElementsAndRecurrence({{StandardProduct3::SatelliteGroup::General, 74},
                             SP3Files::Sentinel3A()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 7),
                    Property(&OrbitRecurrence::Cᴛₒ, 27)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(7'177_⑴ * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(98.63_⑴ * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0011_⑴));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(90.01_⑴ * Degree));

  // The nominal mean solar times of the nodes are 22:00 ascending, 10:00
  // descending.
  EXPECT_THAT(ground_track.mean_solar_time_of_ascending_node()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(21.987_⑴ * Hour));
  EXPECT_THAT(ground_track.mean_solar_time_of_descending_node()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(09.987_⑴ * Hour));
  EXPECT_THAT(ground_track.mean_solar_time_of_ascending_node()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(1.70_⑴ * Second));
  EXPECT_THAT(ground_track.mean_solar_time_of_descending_node()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(1.63_⑴ * Second));
}

#endif

}  // namespace astronomy
}  // namespace principia
