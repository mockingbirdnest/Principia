
#include <string>
#include <vector>
#include <utility>

#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace astronomy {

using base::make_not_null_unique;
using base::not_null;
using geometry::Instant;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::BodySurfaceDynamicFrame;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::Mod;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::ArcMinute;
using quantities::si::ArcSecond;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::IsNear;
using testing_utilities::IsOk;
using ::testing::AllOf;
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

  std::pair<OrbitalElements, OrbitRecurrence> ElementsAndRecurrence(
      SP3Orbit const& orbit) {
    auto const elements =
        OrbitalElements::ForTrajectory(*EarthCentredTrajectory(orbit),
                                       earth_,
                                       MasslessBody{}).ValueOrDie();
    auto const recurrence =
        OrbitRecurrence::ClosestRecurrence(elements.nodal_period(),
                                           elements.nodal_precession(),
                                           earth_,
                                           /*max_abs_Cᴛₒ=*/100);
    return {elements, recurrence};
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

// COSPAR ID 2010-001A, SVN C003.
// 北斗二號 GEO01.
// PRN C01, GEO, 140.0° E.
TEST_F(OrbitAnalysisTest, 北斗GEO) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::北斗, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'166 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(1.42 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.000186));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(166 * Degree));
}

// COSPAR ID 2010-036A, SVN C005.
// 北斗二號 IGSO01.
// PRN C06, IGSO, 117°E.
TEST_F(OrbitAnalysisTest, 北斗IGSO) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::北斗, 6}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'161 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(54.19 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(), IsNear(0.0078));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(232 * Degree));
}

// COSPAR ID 2010-045A, SVN J001.
// Block I-Q, みちびき初号機.
// PRN J01, quasi-zenith orbit.
TEST_F(OrbitAnalysisTest, みちびきQZO) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::みちびき, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  // Expected orbital elements from the Quasi-Zenith Satellite System
  // Performance Standard (PS-QZSS-001).
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'165 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(41 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(), IsNear(0.075));
  // The operational range is 270° ± 2.5°.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(270 * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(0.12 * Degree));
}

// COSPAR ID 2017-048A, SVN J003.
// Block II-G, みちびき3号機.
// PRN J07, GEO.
TEST_F(OrbitAnalysisTest, みちびきGEO) {
  auto j07_files = SP3Files::GNSS();
  // J07 is missing from the last two files.
  j07_files.names.resize(j07_files.names.size() - 2);
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::みちびき, 7}, j07_files});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'166 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(0.067 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.00023));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(224 * Degree));
}

// COSPAR ID 2018-078B, SVN C216.
// 北斗三號 MEO15 (Shanghai Engineering Center for Microsatellites).
// PRN C34, slot A-7.
TEST_F(OrbitAnalysisTest, 北斗MEO) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::北斗, 34}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, -1),
                    Property(&OrbitRecurrence::Cᴛₒ, 7)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(27'906 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(55.01 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.000558));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(1.003 * Degree));
}

// COSPAR ID 2016-030A.
// Galileo-Full Operational Capability Flight Model 10 (GSAT0210) “Danielė”.
// PRN E01, slot A02.
TEST_F(OrbitAnalysisTest, GalileoNominalSlot) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
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

  EXPECT_THAT(elements.nodal_precession(), IsNear(nominal_nodal_precession));
  EXPECT_THAT(2 * π * Radian / elements.anomalistic_period(),
              IsNear(nominal_anomalistic_mean_motion));

  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(29'599.8 * Kilo(Metre), 1.000'03));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(00'000.084 * Kilo(Metre)));

  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.000'17));  // Nominal: 0.0.
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(0.000'015));

  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(56.0 * Degree, 1.03));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(00.01 * Degree));

  EXPECT_THAT(
      Mod(elements.mean_longitude_of_ascending_node_interval().midpoint() -
              nominal_nodal_precession * (mean_time - reference_epoch),
          2 * π * Radian),
      IsNear(317.632 * Degree, 1.000'6));

  // Note that the reference parameters have e = 0, and therefore conventionally
  // set ω = 0, ω′ = 0.
  // However, e is never quite 0; we can compute a mean ω.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(88 * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(6.3 * Degree));

  // Since the reference parameters conventionally set ω = 0, the given mean
  // anomaly is actually the mean argument of latitude; in order to get numbers
  // consistent with theirs, we must look at ω + M.
  EXPECT_THAT(Mod(elements.mean_elements().front().argument_of_periapsis +
                      elements.mean_elements().front().mean_anomaly -
                      nominal_anomalistic_mean_motion *
                          (initial_time - reference_epoch),
                  2 * π * Radian),
              IsNear(225.153 * Degree, 1.005));
}

// COSPAR ID 2014-050B, SVN E202
// Galileo-Full Operational Capability Flight Model 2 (GSAT0202) “Milena”.
// PRN E14, slot Ext02.
TEST_F(OrbitAnalysisTest, GalileoExtendedSlot) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
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

  EXPECT_THAT(elements.nodal_precession(), IsNear(nominal_nodal_precession));
  EXPECT_THAT(2 * π * Radian / elements.anomalistic_period(),
              IsNear(nominal_anomalistic_mean_motion));

  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(27'977.6 * Kilo(Metre), 1.000'01));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(00'000.096 * Kilo(Metre)));

  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.162, 1.06));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(0.000'15));

  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(49.850 * Degree, 1.04));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(00.0044 * Degree));

  EXPECT_THAT(
      Mod(elements.mean_longitude_of_ascending_node_interval().midpoint() -
              nominal_nodal_precession * (mean_time - reference_epoch),
          2 * π * Radian),
      IsNear(52.521 * Degree, 1.02));
  EXPECT_THAT(
      Mod(elements.mean_argument_of_periapsis_interval().midpoint() -
              nominal_apsidal_precession * (mean_time - reference_epoch),
          2 * π * Radian),
      IsNear(56.198 * Degree, 1.02));

  EXPECT_THAT(Mod(elements.mean_elements().front().mean_anomaly -
                      nominal_anomalistic_mean_motion *
                          (initial_time - reference_epoch),
                  2 * π * Radian),
              IsNear(136.069 * Degree, 1.04));
}

// COSPAR ID 2009-070A, SVN R730.
// ГЛОНАСС-М Космос 2456, Ураган-М № 730.
// PRN R01, plane 1.
TEST_F(OrbitAnalysisTest, ГЛОНАСС) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::ГЛОНАСС, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, 1),
                    Property(&OrbitRecurrence::Cᴛₒ, 8)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(25'507 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(64.20 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.00040));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(330 * Degree));
}

// COSPAR ID 2011-036A, SVN G063.
// GPS block IIF satellite.
// PRN G01, plane D, slot 2.
TEST_F(OrbitAnalysisTest, GPS) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::GPS, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(26'560 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(55.86 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(), IsNear(0.0086));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(39 * Degree));
}

// COSPAR ID 1992-052A, TOPEX/Poséidon.
TEST_F(OrbitAnalysisTest, TOPEXPoséidon) {
  auto [elements, recurrence] =  // NOLINT(whitespace/braces)
      ElementsAndRecurrence({{StandardProduct3::SatelliteGroup::General, 1},
                             SP3Files::TOPEXPoséidon()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 13),
                    Property(&OrbitRecurrence::Dᴛₒ, -3),
                    Property(&OrbitRecurrence::Cᴛₒ, 10)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(7'714 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(66.03 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(), IsNear(9.9e-5));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(86.62 * Degree));
}

// COSPAR ID 2002-021A, SPOT-5 (Satellite Pour l’Observation de la Terre).
TEST_F(OrbitAnalysisTest, SPOT5) {
  auto [elements, recurrence] = ElementsAndRecurrence(  // NOLINT
      {{StandardProduct3::SatelliteGroup::General, 94}, SP3Files::SPOT5()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 5),
                    Property(&OrbitRecurrence::Cᴛₒ, 26)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(7'200 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(98.73 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(), IsNear(0.0012));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(89.38 * Degree));
}

// COSPAR ID 2016-011A, Sentinel-3A.
TEST_F(OrbitAnalysisTest, Sentinel3A) {
  auto [elements, recurrence] =  // NOLINT(whitespace/braces)
      ElementsAndRecurrence({{StandardProduct3::SatelliteGroup::General, 74},
                             SP3Files::Sentinel3A()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 7),
                    Property(&OrbitRecurrence::Cᴛₒ, 27)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(7'177 * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(98.63 * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(), IsNear(0.0011));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(90.01 * Degree));
}

}  // namespace astronomy
}  // namespace principia
