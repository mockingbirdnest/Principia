
#pragma once

#include <string>

#include "astronomy/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_fake_plugin {

using astronomy::ICRFJ2000Equator;
using physics::KeplerianElements;
using physics::SolarSystem;

class FakePlugin : public Plugin {
 public:
  // Creates a test plugin with the bodies of the given |solar_system|.  The
  // system must be the Sol system.
  explicit FakePlugin(SolarSystem<ICRFJ2000Equator> const& solar_system);

  // Adds an unloaded vessel with a single part with the given osculating
  // elements around the Earth at |CurrentTime()|.
  Vessel& AddVesselInEarthOrbit(GUID const& vessel_id,
                                std::string const& vessel_name,
                                PartId part_id,
                                std::string const& part_name,
                                KeplerianElements<Barycentric> const& elements);
};

}  // namespace internal_fake_plugin

using internal_fake_plugin::FakePlugin;

}  // namespace ksp_plugin
}  // namespace principia
