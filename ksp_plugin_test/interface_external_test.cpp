
#include "ksp_plugin/interface.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin_test/test_plugin.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace interface {

using astronomy::ICRFJ2000Equator;
using base::make_not_null_unique;
using ksp_plugin::GUID;
using ksp_plugin::Navigation;
using ksp_plugin::PartId;
using ksp_plugin::TestPlugin;
using ksp_plugin::Vessel;
using physics::SolarSystem;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Tonne;
using testing_utilities::SolarSystemFactory;
using ::testing::Eq;

namespace {

constexpr PartId part_id = 789;
constexpr char const* vessel_guid = "123-456";
constexpr char const* part_name = "Picard's desk";
constexpr char const* vessel_name = "NCC-1701-D";

}  // namespace

class InterfaceExternalTest : public ::testing::Test {
 protected:
  InterfaceExternalTest()
      : plugin_(SolarSystem<ICRFJ2000Equator>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt")) {
    physics::KeplerianElements<Barycentric> low_earth_orbit;
    low_earth_orbit.eccentricity = 0;
    low_earth_orbit.semimajor_axis = 6783 * Kilo(Metre);
    low_earth_orbit.inclination = 53.6 * Degree;
    low_earth_orbit.longitude_of_ascending_node = 0 * Radian;
    low_earth_orbit.argument_of_periapsis = 0 * Radian;
    low_earth_orbit.mean_anomaly = 0 * Radian;
    vessel_ = &plugin_.AddVesselInEarthOrbit(
        vessel_guid, vessel_name, part_id, part_name, low_earth_orbit);
  }

  TestPlugin plugin_;
  Vessel* vessel_;
};

TEST_F(InterfaceExternalTest, GetNearestPlannedCoastDegreesOfFreedom) {
  plugin_.CreateFlightPlan(
      vessel_guid, plugin_.CurrentTime() + 24 * Hour, 1 * Tonne);
  vessel_->flight_plan().Append(ksp_plugin::Burn{
      1 * Kilo(Newton),
      1 * Kilo(Newton) * Second / Tonne,
      plugin_.NewBodyCentredNonRotatingNavigationFrame(
          SolarSystemFactory::Earth),
      plugin_.CurrentTime() + 30 * Second,
      Velocity<Frenet<Navigation>>(
          {100 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second}),
      /*is_inertially_fixed=*/false});
  QP result;
  auto const status = principia__ExternalGetNearestPlannedCoastDegreesOfFreedom(
      &plugin_,
      SolarSystemFactory::Earth,
      vessel_guid,
      /*manoeuvre_index=*/0,
      /*reference_position=*/{0, 7'000'000, 0},
      &result);
  EXPECT_THAT(status.error, Eq(static_cast<int>(base::Error::OK)));
}

}  // namespace interface
}  // namespace principia
