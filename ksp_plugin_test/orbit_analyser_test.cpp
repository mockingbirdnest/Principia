
#include "ksp_plugin/orbit_analyser.hpp"

#include <string>
#include <vector>

#include "astronomy/standard_product_3.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/integrators.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace ksp_plugin {

using astronomy::ITRS;
using astronomy::OrbitRecurrence;
using astronomy::StandardProduct3;
using base::not_null;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodySurfaceDynamicFrame;
using physics::Ephemeris;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::astronomy::TerrestrialEquatorialRadius;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using testing_utilities::IsNear;
using testing_utilities::operator""_⑴;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::IsNull;
using ::testing::Optional;
using ::testing::Property;

class OrbitAnalyserTest : public testing::Test {
 protected:
  OrbitAnalyserTest()
      : earth_1957_(RemoveAllButEarth(SolarSystem<Barycentric>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2436116_311504629.proto.txt",
            /*ignore_frame=*/true))),
        // As there is only one body left in this ephemeris, integration is
        // exact up to roundoff error regardless of the step size.  Use a long
        // step so that we do not waste time computing the ephemeris.
        ephemeris_(earth_1957_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<Barycentric>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<Barycentric>>(),
                /*step=*/1 * JulianYear))),
        earth_(*earth_1957_.rotating_body(*ephemeris_, "Earth")),
        itrs_(ephemeris_.get(), &earth_),
        topex_poséidon_(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                            "grgtop03.b97344.e97348.D_S.sp3",
                        StandardProduct3::Dialect::GRGS) {}

  SolarSystem<Barycentric> earth_1957_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> ephemeris_;
  RotatingBody<Barycentric> const& earth_;
  BodySurfaceDynamicFrame<Barycentric, ITRS> itrs_;
  StandardProduct3 topex_poséidon_;

 private:
  SolarSystem<Barycentric> RemoveAllButEarth(
      SolarSystem<Barycentric> solar_system) {
    std::vector<std::string> const names = solar_system.names();
    for (auto const& name : names) {
      if (name != "Earth") {
        solar_system.RemoveMassiveBody(name);
      }
    }
    return solar_system;
  }
};

TEST_F(OrbitAnalyserTest, TOPEXPoséidon) {
  OrbitAnalyser analyser(ephemeris_.get(), DefaultHistoryParameters());
  EXPECT_THAT(analyser.analysis(), IsNull());
  EXPECT_THAT(analyser.progress_of_next_analysis(), Eq(0));
  auto const& arc =
      *topex_poséidon_.orbit(
          {StandardProduct3::SatelliteGroup::General, 1}).front();
  ephemeris_->Prolong(arc.begin()->time);
  analyser.RequestAnalysis(arc.begin()->time,
                           itrs_.FromThisFrameAtTime(arc.begin()->time)(
                               arc.begin()->degrees_of_freedom),
                           3 * Hour);
  while (analyser.progress_of_next_analysis() != 1) {
    absl::SleepFor(absl::Milliseconds(10));
  }
  // Since |progress_of_next_analysis| only tracks the integration, not the
  // analysis, we have no guarantee that an analysis is available immediately.
  do {
    absl::SleepFor(absl::Milliseconds(10));
    analyser.RefreshAnalysis();
  } while (analyser.analysis() == nullptr);
  EXPECT_THAT(analyser.analysis()
                  ->elements()
                  ->mean_semimajor_axis_interval()
                  .midpoint(),
              IsNear(7714_⑴ * Kilo(Metre)));
  EXPECT_THAT(analyser.analysis()->recurrence(),
              Optional(AllOf(Property(&OrbitRecurrence::νₒ, 13),
                             Property(&OrbitRecurrence::Dᴛₒ, -3),
                             Property(&OrbitRecurrence::Cᴛₒ, 10))));
  EXPECT_THAT(analyser.analysis()
                      ->equatorial_crossings()
                      ->longitudes_reduced_to_pass(1)
                      .measure() *
                  TerrestrialEquatorialRadius / Radian,
              IsNear(93_⑴ * Metre));
  // [13; -1; 3] is the subcycle of [13; -3; 10].
  analyser.analysis()->SetRecurrence({13, -1, 3});
  EXPECT_THAT(analyser.analysis()
                      ->equatorial_crossings()
                      ->longitudes_reduced_to_pass(1)
                      .measure() *
                  TerrestrialEquatorialRadius / Radian,
              IsNear(8211_⑴ * Metre));
  // Back to the auto-detected recurrence.
  analyser.analysis()->ResetRecurrence();
  EXPECT_THAT(analyser.analysis()->recurrence(),
              Optional(AllOf(Property(&OrbitRecurrence::νₒ, 13),
                             Property(&OrbitRecurrence::Dᴛₒ, -3),
                             Property(&OrbitRecurrence::Cᴛₒ, 10))));
}

}  // namespace ksp_plugin
}  // namespace principia
