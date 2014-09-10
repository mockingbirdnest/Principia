#pragma once

#include <map>
#include <memory>
#include <string>

#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::Displacement;
using geometry::Instant;

// Universal time 0, time of game creation.
// Putting the origin here makes the instants we use equal as 8-byte objects to
// the corresponding KSP universal time |double|s.
Instant const kUniversalTimeEpoch;

// Thank's to KSP madness, the reference frame of the celestial body orbited by
// the active vessel, occasionally rotating with its surface, occasionally
// nonrotating.
// The basis is that of Unity's "world space" (this is an indirect basis).
struct World;

// Same as |World| but with the y and z axes switched through the looking-glass.
// "We're all mad here. I'm mad. You're mad."
struct AliceWorld;

// The barycentric reference frame of the solar system.
// The basis is the basis of |World| at |kUniversalTimeEpoch|.
struct Universe; 
// The position of the sun at the instant |initial_time| passed at construction.
Position<Universe> kInitialSunPosition;

// A nonrotating referencence frame whose axes are those of |AliceWorld|.
// Since it is nonrotating (though uniformly accelerating), differences between
// velocities are consistent with those in an inertial reference frame.
// When |AliceWorld| rotates this means the axes are not fixed in the reference
// frame, so this pair is inconsistent across instants. Only perform operations
// between simultaneous quantities, then convert to a consistent (frame, basis)
// pair before use.
struct InconsistentNonRotating;

class Plugin {
 public:
  // Creates a |Plugin|. The current time of that instance is |initial_time|.
  // Inserts a celestial body with an arbitrary position, index |sun_index| and
  // gravitational parameter |sun_gravitational_parameter|.
  // The arguments correspond to KSP's |Planetarium.GetUniversalTime()|,
  // |Planetarium.fetch.Sun.flightGlobalsIndex|,
  // |Planetarium.fetch.Sun.gravParameter|.
  Plugin(Instant const& initial_time, int const sun_index,
         GravitationalParameter const& sun_gravitational_parameter);

  // Insert a new celestial body with index |index| and gravitational parameter
  // |gravitational_parameter|. No body with index |index| should already have
  // been inserted. The parent of the new body is the body at index |parent|,
  // which should already have been inserted. The state of the new body at
  // current time is given by |InconsistentNonRotating| offsets from the parent.
  // For a KSP |CelestialBody| |b|, the arguments correspond to:
  // |b.flightGlobalsIndex|, |b.gravParameter|,
  // |b.orbit.referenceBody.flightGlobalsIndex|, |b.orbit.pos|, |b.orbit.vel|.
  void InsertCelestial(
      int const index,
      GravitationalParameter const& gravitational_parameter,
      int const parent,
      Displacement<InconsistentNonRotating> const& from_parent_position,
      Velocity<InconsistentNonRotating> const& from_parent_velocity);

  // Sets the parent of the celestial body with index |index| to the one with
  // index |parent|. Both bodies should already have been inserted.
  // For a KSP |CelestialBody| |b|, the arguments correspond to
  // |b.flightGlobalsIndex| and |b.orbit.referenceBody.flightGlobalsIndex|.
  void UpdateCelestialHierarchy(int index, int parent);

  // Insert a new vessel with GUID |guid| if it does not already exist, and
  // flags the vessel with GUID |guid| so it is kept when calling |AdvanceTime|.
  // The parent body for the vessel is set to the one with index |index|. It
  // should already have been inserted using |InsertCelestial|. Returns |true|
  // if a new vessel was inserted. In that case, |SetVesselStateOffset| should
  // be called with the same |guid| the next call to |AdvanceTime|, so that the
  // initial state of the new vessel is known.
  // For a KSP |Vessel| |v|, the arguments correspond to |v.id|,
  // |v.orbit.referenceBody|.
  bool InsertOrKeepVessel(std::string guid, int parent);

  // Set the position and velocity of the vessel with GUID |guid| relative to
  // its parent at current time. |SetVesselStateOffset| should only be called
  // once per vessel.
  // For a KSP |Vessel| |v|, the arguments correspond to |v.id|, |v.orbit.pos|,
  // |v.orbit.vel|.
  void SetVesselStateOffset(
      std::string guid,
      Displacement<InconsistentNonRotating> from_parent_position,
      Velocity<InconsistentNonRotating> from_parent_velocity);

  // Simulates the system until time |t|. All vessels that have not been
  // refreshed by calling |InsertOrKeepVessel| since the last call to
  //|AdvanceTime| will be removed.
  void AdvanceTime(Instant const& t);

 private:
  struct Celestial;
  struct Vessel;

  std::map<std::string, std::unique_ptr<Vessel>> vessels_;
  std::map<int, std::unique_ptr<Celestial>> celestials_;

  Instant current_time_;
  Celestial* sun_;  // Not owning.
};

}  // namespace ksp_plugin
}  // namespace principia
