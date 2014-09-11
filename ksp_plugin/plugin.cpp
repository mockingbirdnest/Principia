#include "ksp_plugin/plugin.hpp"

#include "geometry/permutation.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::Bivector;
using geometry::Permutation;
using physics::Body;
using physics::Trajectory;

// Represents a KSP |CelestialBody|.
struct Plugin::Celestial {
  explicit Celestial(GravitationalParameter const& gravitational_parameter) 
    : body(new Body(gravitational_parameter)) {}
  std::unique_ptr<Body const> const body;
  // The parent body for the 2-body approximation. Not owning, should only
  // be null for the sun.
  Celestial const* parent = nullptr;
  // The past and present trajectory of the body.
  std::unique_ptr<Trajectory<Barycentre>> history;
};

// Represents a KSP |Vessel|.
struct Plugin::Vessel {
  // Constructs a vessel whose parent is initially |*parent|. |parent| should
  // not be null. No transfer of ownership.
  explicit Vessel(Celestial const* parent) : parent(parent) {
    CHECK(parent != nullptr) << "null parent";
  }
  // A massless |Body|.
  std::unique_ptr<Body const> const body = new Body(GravitationalParameter());
  // The parent body for the 2-body approximation. Not owning, should not be
  // null.
  Celestial const* parent;
  // The past and present trajectory of the body.
  std::unique_ptr<Trajectory<Barycentre>> history;
  // Whether to keep the |Vessel| during the next call to |AdvanceTime|.
  bool keep = true;
};

Permutation<WorldSun, AliceSun> const kWorldLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

Permutation<WorldSun, AliceSun> const kSunLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

Rotation<Barycentre, WorldSun> Plugin::PlanetariumRotation() {
  return Rotation<Barycentre, WorldSun>(
      planetarium_rotation_,
      Bivector<double, Barycentre>({0, 1, 0}));
}

Plugin::Plugin(Instant const& initial_time, int const sun_index,
               GravitationalParameter const& sun_gravitational_parameter,
               Angle const& planetarium_rotation)
    : current_time_(initial_time),
      planetarium_rotation_(planetarium_rotation) {
  celestials_.insert(
      {sun_index, std::make_unique<Celestial>(sun_gravitational_parameter)});
  sun_->history = std::make_unique<Trajectory<Barycentre>>(*sun_->body);
  sun_->history->Append(current_time_,
                        {Position<Barycentre>(), Velocity<Barycentre>()});
}

void Plugin::InsertCelestial(
    int const index,
    GravitationalParameter const& gravitational_parameter,
    int const parent,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun>  const& from_parent_velocity) {
  auto const found_parent = celestials_.find(parent);
  CHECK(found_parent != celestials_.end()) << "No body at index " << parent;
  Celestial const* const parent_body = found_parent->second.get();
  auto const inserted = celestials_.insert(
      {index, std::make_unique<Celestial>(gravitational_parameter)});
  CHECK(inserted.second) << "Body already present at index " << index;
  Celestial *const celestial = inserted.first->second.get();
  celestial->history =
      std::make_unique<Trajectory<Barycentre>>(*celestial->body);
  celestial->history->Append(
      current_time_,
      {parent_body->history->last_position() + PlanetariumRotation().Inverse()(
               kSunLookingGlass.Inverse()(from_parent_position)),
       parent_body->history->last_velocity() + PlanetariumRotation().Inverse()(
               kSunLookingGlass.Inverse()(from_parent_velocity))});
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

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
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
