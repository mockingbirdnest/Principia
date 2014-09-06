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
   // Creates a |Plugin|. The current time of that instance is |inital_time|.
   Plugin(Date const& initial_time);

   // The state of a body described as an offset from the state of the celestial
   // body at index |parent|, with displacement from the parent given by
   // |from_parent_position| and velocity relative to the parent given by
   // |from_parent_velocity|.
   // These members wrap the following expressions for a KSP |Orbit|:
   // |Orbit.referenceBody.flightGlobalsIndex|, |Orbit.pos|, |Orbit.vel|.
   struct CelestialRelativeState {
     int parent;
     Displacement<AliceWorld> from_parent_position;
     VelocityOffset<AliceWorld> from_parent_velocity;
   };

  // Insert a new celestial body with index |index| and gravitational parameter
  // |gravitational_parameter|. No body with index |index| should already have
  // been inserted. If |state| is null the celestial is considered as the sun
  // and given an arbitrary position. Otherwise the parent of the new body
  // is the body at index |state->parent|, which should already have been
  // inserted, and the state of the new body is defined by |state|.
  // |InsertCelestial| should only be called once with null |parent|.
  void InsertCelestial(int index,
                       GravitationalParameter gravitational_parameter,
                       CelestialRelativeState* state);

  // Sets the parent of the celestial body with index |index| to the one with
  // index |parent|. Both bodies should already have been inserted.
  void UpdateCelestialHierarchy(int index, int parent);

  // Insert a new vessel with GUID |guid| if it does not already exist, and
  // flags the vessel with GUID |guid| so it is kept when calling |AdvanceTime|.
  // The parent body for the vessel is set to the one with index |index|. It
  // should already have been inserted using |InsertCelestial|. Returns |true|
  // if a new vessel was inserted. In that case, |SetVesselStateOffset| should
  // be called with the same |guid| the next call to |AdvanceTime|, so that the
  // initial state of the new vessel is known.
  bool InsertOrKeepVessel(std::string guid, int parent);


  // Set the position and velocity of the vessel with GUID |guid| relative to
  // its parent at current time. |SetVesselStateOffset| should only be called
  // once per vessel.
  void SetVesselStateOffset(std::string guid,
                            Displacement<AliceWorld> from_parent_position,
                            VelocityOffset<AliceWorld> from_parent_velocity);

  // Simulates the system until time |t|. All vessels that have not been
  // refreshed by calling |InsertOrKeepVessel| since the last call to
  //|AdvanceTime| will be removed. The initial states of the vessels that remain
  // should have been set using |SetVesselStateOffset|.
  void AdvanceTime(Date const& t)

 private:
  struct Celestial;
  struct Vessel;

  std::map<std::string, std::unique_ptr<Vessel>> vessels_;
  std::map<int, std::unique_ptr<Celestial>> celestials_;

  std::unique_ptr<Date> current_time_;
  std::unique_ptr<Celestial> sun_;
};

}  // namespace ksp_plugin
}  // namespace principia
