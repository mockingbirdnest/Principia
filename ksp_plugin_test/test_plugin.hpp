
#pragma once

#include "astronomy/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_test_plugin {

using astronomy::ICRFJ2000Equator;
using physics::KeplerianElements;
using physics::SolarSystem;

class TestPlugin : public Plugin {
 public:
  TestPlugin(SolarSystem<ICRFJ2000Equator> const& solar_system);

  Vessel& AddVesselInEarthOrbit(GUID const& vessel_id,
                                std::string const& vessel_name,
                                PartId part_id,
                                std::string const& part_name,
                                KeplerianElements<Barycentric> const& elements);
};

}  // namespace internal_test_plugin

using internal_test_plugin::TestPlugin;

}  // namespace ksp_plugin
}  // namespace principia
