#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/orbit_ground_track.hpp"
#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "astronomy/standard_product_3.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/angle_reduction.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace astronomy {

using ::testing::AllOf;
using ::testing::AnyOf;
using ::testing::Field;
using ::testing::Lt;
using ::testing::Property;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_orbit_ground_track;
using namespace principia::astronomy::_orbit_recurrence;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::astronomy::_standard_product_3;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::numerics::_angle_reduction;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;

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
                SymmetricLinearMultistepIntegrator<
                    QuinlanTremaine1990Order12,
                    Ephemeris<ICRS>::NewtonianMotionEquation>(),
                /*step=*/1 * JulianYear))),
        earth_(*earth_1957_.rotating_body(*ephemeris_, "Earth")) {}

  // Returns a GCRS trajectory obtained by stitching together the trajectories
  // of `sp3_orbit.satellites` in `sp3_orbit.files`.
  not_null<std::unique_ptr<DiscreteTrajectory<GCRS>>> EarthCentredTrajectory(
      SP3Orbit const& sp3_orbit) {
    BodyCentredNonRotatingReferenceFrame<ICRS, GCRS> gcrs{ephemeris_.get(),
                                                          &earth_};
    BodySurfaceReferenceFrame<ICRS, ITRS> itrs{ephemeris_.get(), &earth_};

    auto result = make_not_null_unique<DiscreteTrajectory<GCRS>>();
    for (auto const& file : sp3_orbit.files.names) {
      StandardProduct3 sp3(
          SOLUTION_DIR / "astronomy" / "standard_product_3" / file,
          sp3_orbit.files.dialect);
      std::vector<not_null<DiscreteTrajectory<ITRS> const*>> const& orbit =
          sp3.orbit(sp3_orbit.satellite);
      CHECK_EQ(orbit.size(), 1);
      auto const& arc = *orbit.front();
      for (auto const& [time, degrees_of_freedom] : arc) {
        EXPECT_OK(ephemeris_->Prolong(time));
        EXPECT_OK(result->Append(
            time,
            gcrs.ToThisFrameAtTime(time)(
                itrs.FromThisFrameAtTime(time)(degrees_of_freedom))));
      }
    }
    return result;
  }

  std::tuple<OrbitalElements, OrbitRecurrence, OrbitGroundTrack>
  ElementsAndRecurrence(SP3Orbit const& orbit) {
    auto const earth_centred_trajectory = EarthCentredTrajectory(orbit);
    auto elements = OrbitalElements::ForTrajectory(
                        *earth_centred_trajectory,
                        earth_,
                        MasslessBody{},
                        /*fill_osculating_equinoctial_elements=*/true)
                        .value();

    {
      auto const identifier = (std::stringstream() << orbit.satellite).str();
      Logger logger(SOLUTION_DIR / "mathematica" /
                        (identifier + "_elements.generated.wl"),
                    /*make_unique=*/false);
      logger.Set(identifier + "osculatingEquinoctialElements",
                 elements.osculating_equinoctial_elements(),
                 ExpressIn(Metre, Second, Radian));
      logger.Set(identifier + "meanEquinoctialElements",
                 elements.mean_equinoctial_elements(),
                 ExpressIn(Metre, Second, Radian));
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
    // We use values based on Newcomb's formula for the mean longitude of the
    // sun, L = 279°41′48″.04 + 129 602 768″.13 T + 1.089″ T², where T is in
    // Julian centuries since 1900, January 0, Greenwich Mean noon.  Ephemeris
    // time (ET) was defined by the IAU (10th general assembly (1958),
    // commissions 4 and 31, recommendation 2) from Newcomb's tables, with 1900,
    // January 0 at 12 ET being the instant at which the mean longitude of the
    // sun was 279°41′48″.08, and the ET second being 1/31 556 925.9747 of the
    // tropical year at that epoch—where 31 556 925.9747 s =
    // 2π rad / (129 602 768″.13 / 100 a).  TDT, and later TT, were then defined
    // in such a way as to achieve approximate continuity with ET, see 16th
    // general assembly (1976), commission 4, recommendation 5, note 2, and 21st
    // general assembly (1991), resolution A4, recommendation IV, note 4.  We
    // can therefore work with this formula in TT.
    PolynomialInMonomialBasis<Angle, Instant, 2> const
        newcomb_mean_longitude(
            {279 * Degree + 41 * ArcMinute + 48.04 * ArcSecond,
             129'602'768.13 * ArcSecond / (100 * JulianYear),
             1.089 * ArcSecond / Pow<2>(100 * JulianYear)},
            "1899-12-31T12:00:00"_TT);
    auto ground_track = OrbitGroundTrack::ForTrajectory(
        *earth_centred_trajectory,
        earth_,
        {{.epoch = J2000,
          .mean_longitude_at_epoch = newcomb_mean_longitude(J2000),
          .year = 2 * π * Radian /
              newcomb_mean_longitude.EvaluateDerivative(J2000)}}).value();
    return {std::move(elements), recurrence, std::move(ground_track)};
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
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::北斗, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'166_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(1.42_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.000186_(1)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(166_(1) * Degree));
}

// COSPAR ID 2010-036A, SVN C005.
// 北斗二號 IGSO01.
// PRN C06, IGSO, 117°E.
TEST_F(OrbitAnalysisTest, 北斗IGSO) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::北斗, 6}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'161_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(54.19_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0078_(1)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(232_(1) * Degree));
}

// COSPAR ID 2010-045A, SVN J001.
// Block I-Q, みちびき初号機.
// PRN J01, quasi-zenith orbit.
TEST_F(OrbitAnalysisTest, みちびきQZO) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::みちびき, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  // Expected orbital elements from the Quasi-Zenith Satellite System
  // Performance Standard (PS-QZSS-001).
  EXPECT_THAT(
      elements.mean_semimajor_axis_interval().midpoint(),
      AbsoluteErrorFrom(42'165 * Kilo(Metre), IsNear(6.3_(1) * Kilo(Metre))));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(41_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.075_(1)));
  // The operational range is 270° ± 2.5°.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(270_(1) * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(0.12_(1) * Degree));
}

// COSPAR ID 2017-048A, SVN J003.
// Block II-G, みちびき3号機.
// PRN J07, GEO.
TEST_F(OrbitAnalysisTest, みちびきGEO) {
  auto j07_files = SP3Files::GNSS();
  // J07 is missing from the last two files.
  j07_files.names.resize(j07_files.names.size() - 2);
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::みちびき, 7}, j07_files});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 1),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(42'166_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(0.067_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.00023_(1)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(224_(1) * Degree));
}

// COSPAR ID 2018-078B, SVN C216.
// 北斗三號 MEO15 (Shanghai Engineering Center for Microsatellites).
// PRN C34, slot A-7.
TEST_F(OrbitAnalysisTest, 北斗MEO) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::北斗, 34}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, -1),
                    Property(&OrbitRecurrence::Cᴛₒ, 7)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(27'906_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(55.10_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              AnyOf(IsNear(0.000554_(1)),    // Windows, Ubuntu.
                    IsNear(0.000550_(1))));  // macOS.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              AnyOf(IsNear(0.7849_(1) * Degree),      // Windows.
                    IsNear(1.133_(1) * Degree),       // Ubuntu.
                    IsNear(-0.07348_(1) * Degree)));  // macOS.
}

// COSPAR ID 2016-030A.
// Galileo-Full Operational Capability Flight Model 10 (GSAT0210) “Danielė”.
// PRN E01, slot A02.
TEST_F(OrbitAnalysisTest, GalileoNominalSlot) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::Galileo, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, -3),
                    Property(&OrbitRecurrence::Cᴛₒ, 10)));

  // Reference elements from
  // https://www.gsc-europa.eu/system-service-status/orbital-and-technical-parameters.
  Instant const reference_epoch = "2016-11-21T00:00:00"_UTC;
  Instant const initial_time = elements.mean_elements().front().time;
  Instant const mean_time =
      initial_time + (elements.mean_elements().back().time - initial_time) / 2;

  auto const nominal_nodal_precession = -0.02764398 * Degree / Day;
  auto const nominal_anomalistic_mean_motion = 613.72253566 * Degree / Day;

  EXPECT_THAT(
      elements.nodal_precession(),
      AllOf(AbsoluteErrorFrom(nominal_nodal_precession,
                              IsNear(0.00032_(1) * Degree / Day)),
            RelativeErrorFrom(nominal_nodal_precession, IsNear(0.011_(1)))));
  EXPECT_THAT(2 * π * Radian / elements.anomalistic_period(),
              AllOf(AbsoluteErrorFrom(
                        nominal_anomalistic_mean_motion,
                        AnyOf(IsNear(0.66_(1) * Degree / Day),    // Windows.
                              IsNear(0.63_(1) * Degree / Day))),  // macOS.
                    RelativeErrorFrom(nominal_anomalistic_mean_motion,
                                      AnyOf(IsNear(0.00108_(1)),  // Windows.
                                            IsNear(0.00107_(1)),  // Ubuntu.
                                            IsNear(0.00103_(1))))));  // macOS.

  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              AbsoluteErrorFrom(29'599.8 * Kilo(Metre),
                                IsNear(0.33_(1) * Kilo(Metre))));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(00'000.089_(1) * Kilo(Metre)));

  // Nominal: 0.0.
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              AnyOf(IsNear(0.000'17_(1)),    // Windows, Ubuntu.
                    IsNear(0.000'18_(1))));  // macOS.
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              AnyOf(IsNear(0.000'018_(1)),    // Windows.
                    IsNear(0.000'025_(1)),    // Ubuntu.
                    IsNear(0.000'022_(1))));  // macOS.

  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              AbsoluteErrorFrom(56.0 * Degree, IsNear(0.61_(1) * Degree)));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(00.01_(1) * Degree));

  EXPECT_THAT(
      (ReduceAngle<0, 2 * π>(
          elements.mean_longitude_of_ascending_node_interval().midpoint() -
          nominal_nodal_precession * (mean_time - reference_epoch))),
      AbsoluteErrorFrom(317.632 * Degree, IsNear(0.082_(1) * Degree)));

  // Note that the reference parameters have e = 0, and therefore conventionally
  // set ω = 0, ω′ = 0.
  // However, e is never quite 0; we can compute a mean ω.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(89_(1) * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              AnyOf(IsNear(7.8_(1) * Degree),    // Windows.
                    IsNear(7.2_(1) * Degree),    // Ubuntu.
                    IsNear(7.3_(1) * Degree)));  // macOS.

  // Since the reference parameters conventionally set ω = 0, the given mean
  // anomaly is actually the mean argument of latitude; in order to get numbers
  // consistent with theirs, we must look at ω + M.
  EXPECT_THAT(
      (ReduceAngle<0, 2 * π>(
          elements.mean_elements().front().argument_of_periapsis +
          elements.mean_elements().front().mean_anomaly -
          nominal_anomalistic_mean_motion * (initial_time - reference_epoch))),
      AbsoluteErrorFrom(225.153 * Degree, IsNear(0.53_(1) * Degree)));
}

// COSPAR ID 2014-050B, SVN E202
// Galileo-Full Operational Capability Flight Model 2 (GSAT0202) “Milena”.
// PRN E14, slot Ext02.
TEST_F(OrbitAnalysisTest, GalileoExtendedSlot) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
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
                              IsNear(0.00023_(1) * Degree / Day)),
            RelativeErrorFrom(nominal_nodal_precession, IsNear(0.0059_(1)))));
  EXPECT_THAT(2 * π * Radian / elements.anomalistic_period(),
              AllOf(AbsoluteErrorFrom(nominal_anomalistic_mean_motion,
                                      IsNear(0.0011_(1) * Degree / Day)),
                    RelativeErrorFrom(nominal_anomalistic_mean_motion,
                                      IsNear(1.7e-06_(1)))));

  EXPECT_THAT(
      elements.mean_semimajor_axis_interval().midpoint(),
      AbsoluteErrorFrom(27'977.6 * Kilo(Metre),
                        AnyOf(IsNear(0.0485_(1) * Kilo(Metre)),     // Windows.
                              IsNear(0.0531_(1) * Kilo(Metre)),     // Ubuntu.
                              IsNear(0.0534_(1) * Kilo(Metre)))));  // macOS.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              AnyOf(IsNear(00'000.101_(1) * Kilo(Metre)),    // Windows.
                    IsNear(00'000.099_(1) * Kilo(Metre)),    // Ubuntu.
                    IsNear(00'000.098_(1) * Kilo(Metre))));  // macOS.

  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              AbsoluteErrorFrom(0.162, IsNear(0.0041_(1))));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(0.000'15_(1)));

  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              AbsoluteErrorFrom(49.850 * Degree, IsNear(0.77_(1) * Degree)));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(00.0044_(1) * Degree));

  EXPECT_THAT(
      (ReduceAngle<0, 2 * π>(
          elements.mean_longitude_of_ascending_node_interval().midpoint() -
          nominal_nodal_precession * (mean_time - reference_epoch))),
      AbsoluteErrorFrom(52.521 * Degree, IsNear(0.29_(1) * Degree)));
  EXPECT_THAT((ReduceAngle<0, 2 * π>(
                  elements.mean_argument_of_periapsis_interval().midpoint() -
                  nominal_apsidal_precession * (mean_time - reference_epoch))),
              AbsoluteErrorFrom(56.198 * Degree, IsNear(0.48_(1) * Degree)));

  EXPECT_THAT(
      (ReduceAngle<0, 2 * π>(
          elements.mean_elements().front().mean_anomaly -
          nominal_anomalistic_mean_motion * (initial_time - reference_epoch))),
      AbsoluteErrorFrom(136.069 * Degree, IsNear(2.5_(1) * Degree)));
}

// COSPAR ID 2009-070A, SVN R730.
// ГЛОНАСС-М Космос 2456, Ураган-М № 730.
// PRN R01, plane 1.
TEST_F(OrbitAnalysisTest, ГЛОНАСС) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::ГЛОНАСС, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, 1),
                    Property(&OrbitRecurrence::Cᴛₒ, 8)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(25'507_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(64.20_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.00040_(1)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(330_(1) * Degree));
}

// COSPAR ID 2011-036A, SVN G063.
// GPS block IIF satellite.
// PRN G01, plane D, slot 2.
TEST_F(OrbitAnalysisTest, GPS) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::GPS, 1}, SP3Files::GNSS()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 2),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(26'560_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(55.86_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0086_(1)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              IsNear(39_(1) * Degree));
}

// COSPAR ID 1992-052A, TOPEX/Poséidon.
TEST_F(OrbitAnalysisTest, TOPEXPoséidon) {
  // The references for these orbital characteristics are [BSFL98], [Ben97] and
  // [BS96].

  auto const [elements, recurrence, ground_track] =
      ElementsAndRecurrence({{StandardProduct3::SatelliteGroup::General, 1},
                             SP3Files::TOPEXPoséidon()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 13),
                    Property(&OrbitRecurrence::Dᴛₒ, -3),
                    Property(&OrbitRecurrence::Cᴛₒ, 10)));

  // Reference semimajor axis from the legend of figure 7 of [BSFL98]; that
  // value is given as 7 714.42942 in table 1 of [BSFL98], 7 714.4278 km in
  // [BS96], and 7 714.43 km in [Ben97].
  // Figure 7 of [BSFL98] shows that the mean semimajor axis varies by up to
  // 6 m either side of the reference value.  Our test data is from cycle 198;
  // reading the graph around that time shows that the mean semimajor axis was a
  // bit less than 3 m above the nominal value around that time.
  EXPECT_THAT(
      elements.mean_semimajor_axis_interval().midpoint(),
              DifferenceFrom(7714.42938 * Kilo(Metre),
                             AnyOf(IsNear(2.63_(1) * Metre),     // Windows.
                                   IsNear(2.52_(1) * Metre),     // Ubuntu.
                                   IsNear(2.34_(1) * Metre))));  // macOS.
  // Reference inclination from the legend of figure 9 of [BSFL98]; that
  // value is given as 66.040° in table 1 of [BSFL98], 66.039° in [BS96], and
  // 66.04° in [Ben97].
  // Figure 9 of [BSFL98] shows that the observed values of the mean
  // inclination vary between 66.0375° and 66.0455° over the course of the
  // mission (66.0408° ± 0.0040° according to the abstract).  We can see from
  // the graph that during the period under test, in early December 1997, the
  // inclination is below the reference value, and varies by several thousandths
  // of a degree in a matter of days; the minimum and maximum that we compute
  // below therefore seem plausible.
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              DifferenceFrom(66.0408 * Degree, IsNear(-0.0025_(1) * Degree)));
  EXPECT_THAT(
      elements.mean_inclination_interval(),
      AllOf(Field(&Interval<Angle>::min, IsNear(66.0365_(1) * Degree)),
            Field(&Interval<Angle>::max, IsNear(66.0401_(1) * Degree))));

  // The same nominal values are given by [BS96], [Ben97], and [BSFL98]:
  // e = 0.000095, ω = 90.0°.
  // However, these values are only meaningful if we take mean elements that are
  // free of variations over one 127-revolution ground track cycle (rather than
  // simply over one revolution).  Figure 8 of [BSFL98] shows that both the
  // theoretical and observed mean e and ω vary between 40 ppm and 140 ppm, and
  // between 60° and 120°, respectively.
  EXPECT_THAT(elements.mean_eccentricity_interval(),
              AllOf(Field(&Interval<double>::min,
                          AnyOf(IsNear(83e-6_(1)),    // Windows.
                                IsNear(88e-6_(1)),    // Ubuntu.
                                IsNear(88e-6_(1)))),  // macOS.
                    Field(&Interval<double>::max,
                          AnyOf(IsNear(109e-6_(1)),      // Windows, macOS.
                                IsNear(112e-6_(1))))));  // Ubuntu.
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval(),
              AllOf(Field(&Interval<Angle>::min,
                          AnyOf(IsNear(73.8_(1) * Degree),    // Windows.
                                IsNear(74.0_(1) * Degree),    // Ubuntu.
                                IsNear(74.7_(1) * Degree))),  // macOS.
                    Field(&Interval<Angle>::max,
                          AnyOf(IsNear(99.2_(1) * Degree),      // Windows.
                                IsNear(98.9_(1) * Degree),      // Ubuntu.
                                IsNear(98.8_(1) * Degree)))));  // macOS.

  // Nominal longitude of the equatorial crossing of the first ascending pass
  // East of the ITRF zero-meridian (pass 135), as given in section 2 of
  // [Ben97].  [BS96] round these longitudes to a hundredth of a degree, thus
  // 0.71° for pass 135.
  // We can see from figure 6 of [BSFL98] that, during the period under test,
  // the equatorial crossing is about 600 m east of the reference.
  EXPECT_THAT(
      ground_track
          .equator_crossing_longitudes(recurrence,
                                       /*first_ascending_pass_index=*/135)
          .longitudes_reduced_to_pass(135)
          .midpoint(),
      DifferenceFrom(0.7117 * Degree,
                     AllOf(IsNear(0.0051_(1) * Degree),
                           IsNear(573_(1) * Metre *
                                  (Radian / TerrestrialEquatorialRadius)))));
  // Nominal longitude of the equatorial crossing of the following (descending)
  // pass (pass 136), as given in section 2 of [Ben97].  [BS96] round these
  // longitudes to a hundredth of a degree, thus 166.54° for pass 136.
  EXPECT_THAT(ground_track
                  .equator_crossing_longitudes(
                      recurrence, /*first_ascending_pass_index=*/135)
                  .longitudes_reduced_to_pass(136)
                  .midpoint(),
              DifferenceFrom(166.5385 * Degree, IsNear(0.0071_(1) * Degree)));

  // Nominal longitude of the equatorial crossing of pass 1, as given in the
  // auxiliary data table in [BS96].  The reference grid there lists that
  // longitude as 99.92°, and the table of equator crossing longitudes in
  // [Ben97] lists it as 99.9249°.  However, the auxiliary data table in [Ben97]
  // gives a longitude of 99.947° for pass 1, which looks like a typo.
  EXPECT_THAT(ground_track
                  .equator_crossing_longitudes(
                      recurrence, /*first_ascending_pass_index=*/135)
                  .longitudes_reduced_to_pass(1)
                  .midpoint(),
              DifferenceFrom(99.9242 * Degree, IsNear(0.0052_(1) * Degree)));

  // Variability over the period under test (3.5 days).
  EXPECT_THAT(ground_track
                  .equator_crossing_longitudes(
                      recurrence, /*first_ascending_pass_index=*/135)
                  .longitudes_reduced_to_pass(1)
                  .measure(),
              IsNear(0.0025_(1) * Degree));
  EXPECT_THAT(ground_track
                      .equator_crossing_longitudes(
                          recurrence, /*first_ascending_pass_index=*/135)
                      .longitudes_reduced_to_pass(1)
                      .measure() *
                  TerrestrialEquatorialRadius / Radian,
              IsNear(273_(1) * Metre));

  // TOPEX/Poséidon is not sun-synchronous.
  EXPECT_THAT(ground_track.mean_solar_times_of_ascending_nodes()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(42_(1) * Minute));
  EXPECT_THAT(ground_track.mean_solar_times_of_descending_nodes()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(42_(1) * Minute));
}

// The following two satellites are sun-synchronous.

// COSPAR ID 2002-021A, SPOT-5 (Satellite Pour l’Observation de la Terre).
TEST_F(OrbitAnalysisTest, SPOT5) {
  auto const [elements, recurrence, ground_track] = ElementsAndRecurrence(
      {{StandardProduct3::SatelliteGroup::General, 94}, SP3Files::SPOT5()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 5),
                    Property(&OrbitRecurrence::Cᴛₒ, 26)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(7'200_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(98.73_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0012_(1)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              AnyOf(IsNear(89.63_(1) * Degree),    // Windows.
                    IsNear(89.24_(1) * Degree),    // Ubuntu.
                    IsNear(89.47_(1) * Degree)));  // macOS.

  // The nominal mean solar times of the nodes are 22:30 ascending, 10:30
  // descending.
  EXPECT_THAT(ground_track.mean_solar_times_of_ascending_nodes()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(22.452_(1) * Hour));
  EXPECT_THAT(ground_track.mean_solar_times_of_descending_nodes()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(10.452_(1) * Hour));
  EXPECT_THAT(ground_track.mean_solar_times_of_ascending_nodes()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(3.91_(1) * Second));
  EXPECT_THAT(ground_track.mean_solar_times_of_descending_nodes()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(3.91_(1) * Second));
}

// COSPAR ID 2016-011A, Sentinel-3A.
TEST_F(OrbitAnalysisTest, Sentinel3A) {
  auto const [elements, recurrence, ground_track] =
      ElementsAndRecurrence({{StandardProduct3::SatelliteGroup::General, 74},
                             SP3Files::Sentinel3A()});

  EXPECT_THAT(recurrence,
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 7),
                    Property(&OrbitRecurrence::Cᴛₒ, 27)));
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              IsNear(7'177_(1) * Kilo(Metre)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              IsNear(98.63_(1) * Degree));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0011_(1)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              AnyOf(IsNear(90.00_(1) * Degree),    // Windows, Ubuntu.
                    IsNear(90.02_(1) * Degree)));  // macOS.

  // The nominal mean solar times of the nodes are 22:00 ascending, 10:00
  // descending.
  EXPECT_THAT(ground_track.mean_solar_times_of_ascending_nodes()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(21.987_(1) * Hour));
  EXPECT_THAT(ground_track.mean_solar_times_of_descending_nodes()->midpoint() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(09.987_(1) * Hour));
  EXPECT_THAT(ground_track.mean_solar_times_of_ascending_nodes()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(1.70_(1) * Second));
  EXPECT_THAT(ground_track.mean_solar_times_of_descending_nodes()->measure() *
                  (1 * Day / (2 * π * Radian)),
              IsNear(1.63_(1) * Second));
}

#endif

}  // namespace astronomy
}  // namespace principia
