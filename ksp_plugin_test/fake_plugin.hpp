#pragma once

#include <string>

#include "astronomy/frames.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace ksp_plugin_test {
namespace _fake_plugin {
namespace internal {

using namespace principia::astronomy::_frames;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_identification;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_solar_system;

class FakePlugin : public Plugin {
 public:
  // Creates a test plugin with the bodies of the given `solar_system`.  The
  // system must be the Sol system.
  explicit FakePlugin(SolarSystem<ICRS> const& solar_system);

  // Adds an unloaded vessel with a single part with the given osculating
  // elements around the Earth at `CurrentTime()`.
  Vessel& AddVesselInEarthOrbit(GUID const& vessel_id,
                                std::string const& vessel_name,
                                PartId part_id,
                                std::string const& part_name,
                                KeplerianElements<Barycentric> const& elements);
};

}  // namespace internal

using internal::FakePlugin;

}  // namespace _fake_plugin
}  // namespace ksp_plugin_test
}  // namespace principia
