#include "ksp_plugin_test/fake_plugin.hpp"

#include <optional>
#include <string>
#include <type_traits>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "physics/massless_body.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace ksp_plugin_test {
namespace _fake_plugin {
namespace internal {

using namespace principia::physics::_massless_body;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_solar_system_factory;

FakePlugin::FakePlugin(SolarSystem<ICRS> const& solar_system)
    : Plugin(/*game_epoch=*/solar_system.epoch_literal(),
             /*solar_system_epoch=*/solar_system.epoch_literal(),
             /*planetarium_rotation=*/0 * Radian) {
  using SSId = std::underlying_type_t<SolarSystemFactory::Id>;
  absl::flat_hash_map<SSId, Index> ssid_to_index;

  for (SSId ssid = SolarSystemFactory::Sun;
       ssid <= SolarSystemFactory::LastBody;
       ++ssid) {
    std::optional<Index> const parent_index =
        ssid == SolarSystemFactory::Sun
            ? std::nullopt
            : std::make_optional(
                  ssid_to_index[SolarSystemFactory::parent(ssid)]);
    ssid_to_index[ssid] = InsertCelestialAbsoluteCartesian(
        parent_index,
        solar_system.gravity_model_message(SolarSystemFactory::name(ssid)),
        solar_system.cartesian_initial_state_message(
            SolarSystemFactory::name(ssid)));
  }
  EndInitialization();
}

Vessel& FakePlugin::AddVesselInEarthOrbit(
    GUID const& vessel_id,
    std::string const& vessel_name,
    PartId const part_id,
    std::string const& part_name,
    KeplerianElements<Barycentric> const& elements) {
  KeplerOrbit<Barycentric> const earth_orbit(
      *GetCelestial(SolarSystemFactory::Earth).body(),
      MasslessBody{},
      elements,
      CurrentTime());
  auto const barycentric_dof = earth_orbit.StateVectors(CurrentTime());
  auto const alice_dof = PlanetariumRotation()(barycentric_dof);
  bool inserted;
  InsertOrKeepVessel(vessel_id,
                     vessel_name,
                     SolarSystemFactory::Earth,
                     /*loaded=*/false,
                     inserted);
  CHECK(inserted) << "Failed to insert vessel " << vessel_name << " ("
                  << vessel_id << ")";
  InsertUnloadedPart(part_id, part_name, vessel_id, alice_dof);
  PrepareToReportCollisions();
  FreeVesselsAndPartsAndCollectPileUps(20 * Milli(Second));
  return *GetVessel(vessel_id);
}

}  // namespace internal
}  // namespace _fake_plugin
}  // namespace ksp_plugin_test
}  // namespace principia
