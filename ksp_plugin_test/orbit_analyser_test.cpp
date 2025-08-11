#include "ksp_plugin/orbit_analyser.hpp"

#include <memory>
#include <string>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/integrators.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

namespace principia {
namespace ksp_plugin {

using ::testing::AllOf;
using ::testing::Eq;
using ::testing::IsNull;
using ::testing::Optional;
using ::testing::Property;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_orbit_recurrence;
using namespace principia::astronomy::_standard_product_3;
using namespace principia::base::_not_null;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_integrators;
using namespace principia::ksp_plugin::_orbit_analyser;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_astronomy;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;

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
                SymmetricLinearMultistepIntegrator<
                    QuinlanTremaine1990Order12,
                    Ephemeris<Barycentric>::NewtonianMotionEquation>(),
                /*step=*/1 * JulianYear))),
        earth_(*earth_1957_.rotating_body(*ephemeris_, "Earth")),
        itrs_(ephemeris_.get(), &earth_),
        topex_poséidon_(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                            "grgtop03.b97344.e97348.D_S.sp3",
                        StandardProduct3::Dialect::GRGS) {}

  SolarSystem<Barycentric> earth_1957_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> ephemeris_;
  RotatingBody<Barycentric> const& earth_;
  BodySurfaceReferenceFrame<Barycentric, ITRS> itrs_;
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

TEST_F(OrbitAnalyserTest, TinyAnalysis) {
  OrbitAnalyser analyser(ephemeris_.get(), DefaultHistoryParameters());
  EXPECT_THAT(analyser.analysis(), IsNull());
  EXPECT_THAT(analyser.progress_of_next_analysis(), Eq(0));
  auto const& arc =
      *topex_poséidon_.orbit(
          {StandardProduct3::SatelliteGroup::General, 1}).front();
  EXPECT_OK(ephemeris_->Prolong(arc.begin()->time));
  analyser.RequestAnalysis(
      {.first_time = arc.begin()->time,
       .first_degrees_of_freedom = itrs_.FromThisFrameAtTime(arc.begin()->time)(
           arc.begin()->degrees_of_freedom),
       .mission_duration = Abs(J2000 - arc.begin()->time) * 0x1p-45});
  do {
    absl::SleepFor(absl::Milliseconds(10));
    analyser.RefreshAnalysis();
  } while (analyser.analysis() == nullptr);
}


TEST_F(OrbitAnalyserTest, TOPEXPoséidon) {
  OrbitAnalyser analyser(ephemeris_.get(), DefaultHistoryParameters());
  EXPECT_THAT(analyser.analysis(), IsNull());
  EXPECT_THAT(analyser.progress_of_next_analysis(), Eq(0));
  auto const& arc =
      *topex_poséidon_.orbit(
          {StandardProduct3::SatelliteGroup::General, 1}).front();
  EXPECT_OK(ephemeris_->Prolong(arc.begin()->time));
  analyser.RequestAnalysis(
      {.first_time = arc.begin()->time,
       .first_degrees_of_freedom = itrs_.FromThisFrameAtTime(arc.begin()->time)(
           arc.begin()->degrees_of_freedom),
       .mission_duration = 3 * Hour});
  while (analyser.progress_of_next_analysis() != 1) {
    absl::SleepFor(absl::Milliseconds(10));
  }
  // Since `progress_of_next_analysis` only tracks the integration, not the
  // analysis, we have no guarantee that an analysis is available immediately.
  do {
    absl::SleepFor(absl::Milliseconds(10));
    analyser.RefreshAnalysis();
  } while (analyser.analysis() == nullptr);
  EXPECT_THAT(analyser.analysis()
                  ->elements()
                  ->mean_semimajor_axis_interval()
                  .midpoint(),
              IsNear(7714_(1) * Kilo(Metre)));
  EXPECT_THAT(analyser.analysis()->recurrence(),
              Optional(AllOf(Property(&OrbitRecurrence::νₒ, 13),
                             Property(&OrbitRecurrence::Dᴛₒ, -3),
                             Property(&OrbitRecurrence::Cᴛₒ, 10))));
  EXPECT_THAT(analyser.analysis()
                      ->equatorial_crossings()
                      ->longitudes_reduced_to_pass(1)
                      .measure() *
                  TerrestrialEquatorialRadius / Radian,
              IsNear(93_(1) * Metre));
  // [13; -1; 3] is the subcycle of [13; -3; 10].
  analyser.analysis()->SetRecurrence({13, -1, 3});
  EXPECT_THAT(analyser.analysis()
                      ->equatorial_crossings()
                      ->longitudes_reduced_to_pass(1)
                      .measure() *
                  TerrestrialEquatorialRadius / Radian,
              IsNear(8211_(1) * Metre));
  // Back to the auto-detected recurrence.
  analyser.analysis()->ResetRecurrence();
  EXPECT_THAT(analyser.analysis()->recurrence(),
              Optional(AllOf(Property(&OrbitRecurrence::νₒ, 13),
                             Property(&OrbitRecurrence::Dᴛₒ, -3),
                             Property(&OrbitRecurrence::Cᴛₒ, 10))));
}

}  // namespace ksp_plugin
}  // namespace principia
