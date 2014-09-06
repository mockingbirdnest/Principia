#pragma once

#include <map>
#include <memory>
#include <string>

#include "geometry/point.hpp"
#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {

// The "universal" inertial reference frame in which we perform the integration.
struct Universe;

// Unity world space. Thank's to KSP madness, occasionally an inertial reference
// frame, occasionally an uniformly rotating one. Handle with caution.
struct World;

// Unity world space, with the y and z axes switched through the looking-glass.
// "We're all mad here. I'm mad. You're mad."
struct AliceWorld;

typedef geometry::Point<Time> Date;

template<typename Frame>
using Displacement = Vector<Length, Frame>;
template<typename Frame>
using VelocityOffset = Vector<Speed, Frame>;

class Plugin {
 public:
  // Insert a new celestial body with index |index| and gravitational parameter
  // |gravitational_parameter| and null parent. No body with index |index|
  // should already have been inserted.
  void InsertCelestial(int index,
                       GravitationalParameter gravitational_parameter);

  // Update the parent of the celestial body with index |index| to the one with
  // index |parent|. Both bodies should already have been inserted.
  void UpdateCelestialHierarchy(int index, int parent);

  // Insert a new vessel with GUID |guid| if it does not already exist, and
  // spares the vessel with GUID |guid| at the next call to |CleanupVessels|.
  // The parent body for the vessel is set to the one with index |index|. It
  // should already have been inserted using |InsertCelestial|.
  // Returns true if a new vessel was inserted.
  bool InsertOrKeepVessel(std::string guid, int parent);

  // Remove vessels whose GUID was not given as a parameter to a call to
  // |InsertOrKeepVessel| since the last call to |CleanupVessels|.
  void CleanupVessels();

  void SetCurrentTime(Date const& t);

  // Arbitrarily sets the current position of the celestial body with index
  // |index|, so that the |SetCelestialStateOffset| and |SetVesselStateOffset|
  // functions may be used on its children. The position of the celestial body
  // should not have been previously set.
  void DeclareSun(int index);

  // Set the position and velocity of the celestial body with index |index|
  // relative to its parent at current time. The parent's position should have
  // been previously set. The position of the celestial body should not have
  // been previously set.
  void SetCelestialStateOffset(int index,
                               Displacement<AliceWorld> from_parent_position,
                               VelocityOffset<AliceWorld> from_parent_velocity);

  // Set the position and velocity of the vessel with GUID |guid| relative to
  // its parent at current time. The parent's position should have been
  // previously set. The position of the vessel should not have been
  // previously set.
  void SetVesselStateOffset(std::string guid,
                            Displacement<AliceWorld> from_parent_position,
                            VelocityOffset<AliceWorld> from_parent_velocity);

 private:
  // Represents a KSP |CelestialBody|.
  struct Celestial {
    // Takes ownership of |body|.
    explicit Celestial(physics::Body const* const body) : body(body) {}
    std::unique_ptr<physics::Body const> const body;
    // The parent body for the 2-body approximation. Not owning, should only
    // be null for the sun.
    Celestial const* parent = nullptr;
    std::unique_ptr<physics::Trajectory<World>>
  };

  // Represents a KSP |Vessel|.
  struct Vessel {
    // A massless |Body|.
    std::unique_ptr<physics::Body const> const body(
        physics::Body(GravitationalParameter()));
    // The parent body for the 2-body approximation. Not owning, should not be
    // null.
    Celestial const* parent;
    // Whether to keep the |Vessel| during the next cleanup.
    bool keep = true;
  };

  std::map<std::string, std::unique_ptr<Vessel>> vessels_;
  std::map<int, std::unique_ptr<Celestial>> celestials_;

  Date current_time_;
};

}  // namespace ksp_plugin
}  // namespace principia
