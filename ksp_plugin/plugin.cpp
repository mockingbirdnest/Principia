#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace ksp_plugin {

using physics::Body;

// Represents a KSP |CelestialBody|.
struct Plugin::Celestial {
  // Takes ownership of |body|.
  explicit Celestial(physics::Body const* const body) : body(body) {}
  std::unique_ptr<physics::Body const> const body;
  // The parent body for the 2-body approximation. Not owning, should only
  // be null for the sun.
  Celestial const* parent = nullptr;
  // The past and present trajectory of the body.
  std::unique_ptr<physics::Trajectory<World>> history;
};

// Represents a KSP |Vessel|.
struct Plugin::Vessel {
  // Constructs a vessel whose parent is initially |*parent|. |parent| should
  // not be null.
  explicit Vessel(Celestial const* parent) : parent(parent) {
    CHECK(parent != nullptr) << "null parent";
  };
  // A massless |Body|.
  std::unique_ptr<physics::Body const> const body =
      new physics::Body(GravitationalParameter());
  // The parent body for the 2-body approximation. Not owning, should not be
  // null.
  Celestial const* parent;
  // The past and present trajectory of the body.
  std::unique_ptr<physics::Trajectory<World>> history;
  // Whether to keep the |Vessel| during the next cleanup.
  bool keep = true;
};


void Plugin::InsertCelestial(int index,
                             GravitationalParameter gravitational_parameter) {
  auto const inserted = celestials_.insert(
      {index, std::make_unique<Celestial>(new Body(gravitational_parameter))});
  CHECK(inserted.second) << "Multiple bodies bearing the same index";
}

void Plugin::UpdateCelestialHierarchy(int index, int parent) {
  CHECK(celestials_.find(index) != celestials_.end()) <<
      "No body at index " << index;
  CHECK(celestials_.find(parent) != celestials_.end()) <<
      "No body at index " << parent;
  celestials_[index]->parent = celestials_[parent].get();
}

bool Plugin::InsertOrKeepVessel(std::string guid, int parent) {
  CHECK(celestials_.find(parent) != celestials_.end()) <<
      "No body at index " << parent;
  return vessels_.insert(
      {guid, std::make_unique<Vessel>(celestials_[parent].get())}).second;
}

void Plugin::CleanupVessels() {
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    if(!it->second->keep) {
      // |std::map::erase| invalidates its parameter so we post-increment.
      vessels_.erase(it++);
    } else {
      it->second->keep = false;
      ++it;
    }
  }
}

}  // namespace ksp_plugin
}  // namespace principia
