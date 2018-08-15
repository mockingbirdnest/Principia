
#include "ksp_plugin_test/fake_plugin.hpp"

#include <optional>
#include <string>

#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_fake_plugin {

using physics::MasslessBody;
using physics::KeplerOrbit;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::SolarSystemFactory;

FakePlugin::FakePlugin(SolarSystem<ICRS> const& solar_system)
    : Plugin(/*game_epoch=*/solar_system.epoch_literal(),
             /*solar_system_epoch=*/solar_system.epoch_literal(),
             /*planetarium_rotation=*/0 * Radian) {
  for (int index = SolarSystemFactory::Sun;
       index <= SolarSystemFactory::LastBody;
       ++index) {
    std::optional<Index> parent_index =
        index == SolarSystemFactory::Sun
            ? std::nullopt
            : std::make_optional(SolarSystemFactory::parent(index));
    InsertCelestialAbsoluteCartesian(
        index,
        parent_index,
        solar_system.gravity_model_message(SolarSystemFactory::name(index)),
        solar_system.cartesian_initial_state_message(
            SolarSystemFactory::name(index)));
  }
  EndInitialization();
}

Vessel& FakePlugin::AddVesselInEarthOrbit(
    GUID const& vessel_id,
    std::string const& vessel_name,
    PartId const part_id,
    std::string const& part_name,
    KeplerianElements<Barycentric> const& elements) {
  KeplerOrbit<Barycentric> earth_orbit(
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

}  // namespace internal_fake_plugin
}  // namespace ksp_plugin
}  // namespace principia
