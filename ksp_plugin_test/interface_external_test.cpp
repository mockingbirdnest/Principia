#include "ksp_plugin/interface.hpp"

#include <string>
#include <utility>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/vessel.hpp"
#include "ksp_plugin_test/fake_plugin.hpp"
#include "physics/solar_system.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace interface {

using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_identification;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::ksp_plugin_test::_fake_plugin;
using namespace principia::physics::_solar_system;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_solar_system_factory;

namespace {

constexpr PartId part_id = 1729;
constexpr char const* vessel_guid = "NCC 1701-D";
constexpr char const* part_name = "Picard's desk";
constexpr char const* vessel_name = "Enterprise";

MATCHER(IsOk,
        std::string(negation ? "is not" : "is") + " ok") {
  return arg.error == 0;
}

}  // namespace

class InterfaceExternalTest : public ::testing::Test {
 protected:
  InterfaceExternalTest()
      : plugin_(SolarSystem<ICRS>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt")) {
    physics::_kepler_orbit::KeplerianElements<Barycentric> low_earth_orbit;
    low_earth_orbit.eccentricity = 0;
    low_earth_orbit.semimajor_axis = 6783 * Kilo(Metre);
    low_earth_orbit.inclination = 0 * Degree;
    low_earth_orbit.longitude_of_ascending_node = 0 * Radian;
    low_earth_orbit.argument_of_periapsis = 0 * Radian;
    low_earth_orbit.mean_anomaly = 0 * Radian;
    vessel_ = &plugin_.AddVesselInEarthOrbit(
        vessel_guid, vessel_name, part_id, part_name, low_earth_orbit);
  }

  FakePlugin plugin_;
  Vessel* vessel_;
};

TEST_F(InterfaceExternalTest, GetNearestPlannedCoastDegreesOfFreedom) {
  plugin_.CreateFlightPlan(
      vessel_guid, plugin_.CurrentTime() + 6 * Hour, 1 * Tonne);
  NavigationManœuvre::Intensity intensity;
  intensity.Δv = Velocity<Frenet<Navigation>>({1000 * Metre / Second,
                                               0 * Metre / Second,
                                               0 * Metre / Second});
  NavigationManœuvre::Timing timing;
  timing.initial_time =  plugin_.CurrentTime() + 30 * Second;
  NavigationManœuvre::Burn burn{
      intensity,
      timing,
      180 * Kilo(Newton),
      4.56 * Kilo(Newton) * Second / Kilogram,
      plugin_.NewBodyCentredNonRotatingNavigationFrame(
          SolarSystemFactory::Earth),
      /*is_inertially_fixed=*/false};
  EXPECT_OK(vessel_->flight_plan().Insert(std::move(burn), 0));
  QP result;
  auto const to_world =
      plugin_.renderer().BarycentricToWorld(plugin_.PlanetariumRotation());
  auto const* const status =
      principia__ExternalGetNearestPlannedCoastDegreesOfFreedom(
          &plugin_,
          SolarSystemFactory::Earth,
          vessel_guid,
          /*manoeuvre_index=*/0,
          /*reference_position=*/
          ToXYZ(to_world(Displacement<Barycentric>(
                             {-100'000 * Kilo(Metre), 0 * Metre, 0 * Metre}))
                    .coordinates() /
                Metre),
          &result);
  EXPECT_THAT(*status, IsOk());
  auto const barycentric_result =
      to_world.Inverse()(FromQP<RelativeDegreesOfFreedom<World>>(result));
  // The reference position is far above the apoapsis, so the result is roughly
  // the apoapsis.
  EXPECT_THAT(
      barycentric_result,
      Componentwise(Componentwise(IsNear(-11'972_(1) * Kilo(Metre)),
                                  IsNear(-120_(1) * Kilo(Metre)),
                                  AllOf(Gt(-50 * Metre), Lt(50 * Metre))),
                    Componentwise(IsNear(-6.7_(1) * Metre / Second),
                                  IsNear(-4.9_(1) * Kilo(Metre) / Second),
                                  AllOf(Gt(-1 * Centi(Metre) / Second),
                                        Lt(1 * Centi(Metre) / Second)))));
}

TEST_F(InterfaceExternalTest, Geopotential) {
  XY coefficient;
  double radius;
  auto const* status = principia__ExternalGeopotentialGetCoefficient(
      &plugin_,
      SolarSystemFactory::Earth,
      /*degree=*/2,
      /*order=*/0,
      &coefficient);
  EXPECT_THAT(*status, IsOk());
  EXPECT_THAT(-coefficient.x * Sqrt(5), IsNear(1.08e-3_(1)));
  EXPECT_THAT(coefficient.y, Eq(0));

  status = principia__ExternalGeopotentialGetCoefficient(
      &plugin_,
      SolarSystemFactory::Earth,
      /*degree=*/3,
      /*order=*/1,
      &coefficient);
  EXPECT_THAT(*status, IsOk());
  EXPECT_THAT(coefficient.x, IsNear(2.03e-6_(1)));
  EXPECT_THAT(coefficient.y, IsNear(0.248e-6_(1)));

  status = principia__ExternalGeopotentialGetCoefficient(
      &plugin_,
      SolarSystemFactory::Earth,
      /*degree=*/1729,
      /*order=*/163,
      &coefficient);
  EXPECT_THAT(*status, IsOk());
  EXPECT_THAT(coefficient.x, Eq(0));
  EXPECT_THAT(coefficient.y, Eq(0));

  status = principia__ExternalGeopotentialGetReferenceRadius(
      &plugin_,
      SolarSystemFactory::Earth,
      &radius);
  EXPECT_THAT(*status, IsOk());
  EXPECT_THAT(radius, Eq(6'378'136.3));

  status = principia__ExternalGeopotentialGetCoefficient(
      &plugin_,
      SolarSystemFactory::Ariel,
      /*degree=*/2,
      /*order=*/2,
      &coefficient);
  EXPECT_THAT(*status, IsOk());
  EXPECT_THAT(coefficient.x, Eq(0));
  EXPECT_THAT(coefficient.y, Eq(0));

  status = principia__ExternalGeopotentialGetCoefficient(
      &plugin_,
      SolarSystemFactory::Ariel,
      /*degree=*/0,
      /*order=*/0,
      &coefficient);
  EXPECT_THAT(*status, IsOk());
  EXPECT_THAT(coefficient.x, Eq(1));
  EXPECT_THAT(coefficient.y, Eq(0));

  status = principia__ExternalGeopotentialGetReferenceRadius(
      &plugin_,
      SolarSystemFactory::Ariel,
      &radius);
  EXPECT_THAT(*status, IsOk());
  EXPECT_THAT(radius, Eq(578'900));
}

}  // namespace interface
}  // namespace principia
