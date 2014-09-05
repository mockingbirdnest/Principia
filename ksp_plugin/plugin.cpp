#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace ksp_plugin {

using physics::Body;

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
  auto const inserted = vessels_.insert({guid, std::make_unique<Vessel>()});
  inserted.first->second->parent = celestials_[parent].get();
  return inserted.second;
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
