#pragma once

#include <map>
#include <memory>
#include <string>

#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gtest/gtest.h"
#include "ksp_plugin/monostable.hpp"
#include "physics/body.hpp"
#include "physics/n_body_system.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::Displacement;
using geometry::Instant;
using geometry::Rotation;
using integrators::SPRKIntegrator;
using physics::Body;
using physics::Trajectory;
using physics::NBodySystem;
using quantities::Angle;
using si::Second;

// Universal time 0, time of game creation.
// Putting the origin here makes the instants we use equal to the corresponding
// KSP universal time doubles.
Instant const kUniversalTimeEpoch;

// Thanks to KSP's madness, the reference frame of the celestial body orbited by
// the active vessel, occasionally rotating with its surface, occasionally
// nonrotating.
// The basis is that of Unity's "world space" (this is a left-handed basis).
struct World;

// Same as |World| but with the y and z axes switched through the looking-glass:
// it is a right-handed basis. "We're all mad here. I'm mad. You're mad."
struct AliceWorld;

// The barycentric reference frame of the solar system.
// The basis is the basis of |World| at |kUniversalTimeEpoch|.
// TODO(egg): it *should* be the barycentric frame. For the moment we're using
// the velocity of the sun at the time of construction as our reference.
struct Barycentre;
// The position of the sun at the instant |initial_time| passed at construction.
Position<Barycentre> const kInitialSunPosition;

// A nonrotating referencence frame comoving with the sun with the same axes as
// |AliceWorld|. Since it is nonrotating (though not inertial), differences
// between velocities are consistent with those in an inertial reference frame.
// When |AliceWorld| rotates the axes are not fixed in the reference frame, so
// this (frame, basis) pair is inconsistent across instants. Operations should
// only be performed between simultaneous quantities, then converted to a
// consistent (frame, basis) pair before use.
struct AliceSun;

// Same as above, but with same axes as |World| instead of those of
// |AliceWorld|. The caveats are the same as for |AliceSun|.
struct WorldSun;

// The GUID of a vessel, obtained by |v.id.ToString()| in C#. We use this as a
// key in an |std::map|.
using GUID = std::string;
// The index of a body in |FlightGlobals.Bodies|, obtained by
// |b.flightGlobalsIndex| in C#. We use this as a key in an |std::map|.
using Index = int;

class Plugin {
 public:
  Plugin() = delete;
  Plugin(Plugin const&) = delete;
  Plugin(Plugin&&) = delete;
  virtual ~Plugin() = default;

  // Constructs a |Plugin|. The current time of that instance is |initial_time|.
  // The angle between the axes of |World| and |Barycentre| at |initial_time| is
  // set to |planetarium_rotation|. Inserts a celestial body with an arbitrary
  // position, index |sun_index| and gravitational parameter
  // |sun_gravitational_parameter|.
  // Starts initialization.
  // The arguments correspond to KSP's
  // |Planetarium.GetUniversalTime()|,
  // |Planetarium.fetch.Sun.flightGlobalsIndex|,
  // |Planetarium.fetch.Sun.gravParameter|,
  // |Planetarium.InverseRotAngle|.
  Plugin(Instant const& initial_time,
         Index const sun_index,
         GravitationalParameter const& sun_gravitational_parameter,
         Angle const& planetarium_rotation);

  // Inserts a new celestial body with index |celestial_index| and gravitational
  // parameter |gravitational_parameter|. No body with index |celestial_index|
  // must already have been inserted. The parent of the new body is the body
  // at index |parent_index|, which must already have been inserted. The state
  // of the new body at current time is given by |AliceSun| offsets from the
  // parent. Must only be called during initialization.
  // For a KSP |CelestialBody| |b|, the arguments correspond to:
  // |b.flightGlobalsIndex|,
  // |b.gravParameter|,
  // |b.orbit.referenceBody.flightGlobalsIndex|,
  // |b.orbit.pos|,
  // |b.orbit.vel|.
  virtual void InsertCelestial(
    Index const celestial_index,
    GravitationalParameter const& gravitational_parameter,
    Index const parent_index,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun> const& from_parent_velocity);

  // Ends initialization.
  virtual void EndInitialization();

  // Sets the parent of the celestial body with index |celestial_index| to the
  // one with index |parent_index|. Both bodies must already have been
  // inserted. Must be called after initialization.
  // For a KSP |CelestialBody| |b|, the arguments correspond to
  // |b.flightGlobalsIndex|, |b.orbit.referenceBody.flightGlobalsIndex|.
  virtual void UpdateCelestialHierarchy(Index const celestial_index,
                                        Index const parent_index) const;

  // Inserts a new vessel with GUID |vessel_guid| if it does not already exist,
  // and flags the vessel with GUID |vessel_guid| so it is kept when calling
  // |AdvanceTime|. The parent body for the vessel is set to the one with index
  // |parent_index|. It must already have been inserted using
  // |InsertCelestial|.
  // Returns true if a new vessel was inserted. In that case,
  // |SetVesselStateOffset| must be called with the same GUID before the
  // next call to |AdvanceTime|, |VesselDisplacementFromParent| or
  // |VesselParentRelativeVelocity|, so that the initial state of the new
  // vessel is known. Must be called after initialization.
  // For a KSP |Vessel| |v|, the arguments correspond to
  // |v.id|, |v.orbit.referenceBody.flightGlobalsIndex|.
  virtual bool InsertOrKeepVessel(GUID const& vessel_guid,
                                  Index const parent_index);

  // Set the position and velocity of the vessel with GUID |vessel_guid|
  // relative to its parent at current time. |SetVesselStateOffset| must only
  // be called once per vessel. Must be called after initialization.
  // For a KSP |Vessel| |v|, the arguments correspond to
  // |v.id.ToString()|,
  // |v.orbit.pos|,
  // |v.orbit.vel|.
  virtual void SetVesselStateOffset(
      GUID const& vessel_guid,
      Displacement<AliceSun> const& from_parent_position,
      Velocity<AliceSun> const& from_parent_velocity);

  // Simulates the system until instant |t|. All vessels that have not been
  // refreshed by calling |InsertOrKeepVessel| since the last call to
  // |AdvanceTime| will be removed. Must be called after initialization.
  // |planetarium_rotation| is the value of KSP's |Planetarium.InverseRotAngle|
  // at instant |t|, which provides the rotation between the |World| axes and
  // the |Barycentre| axes (we don't use Planetarium.Rotation since it undergoes
  // truncation to single-precision even though it's a double-precision value).
  // Note that KSP's |Planetarium.InverseRotAngle| is in degrees.
  virtual void AdvanceTime(Instant const& t, Angle const& planetarium_rotation);

  // Returns the position of the vessel with GUID |vessel_guid| relative to its
  // parent at current time. For a KSP |Vessel| |v|, the argument corresponds to
  // |v.id.ToString()|, the return value to |v.orbit.pos|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept. Must
  // be called after initialization.
  virtual Displacement<AliceSun> VesselDisplacementFromParent(
      GUID const& vessel_guid) const;

  // Returns the velocity of the vessel with GUID |vessel_guid| relative to its
  // parent at current time. For a KSP |Vessel| |v|, the argument corresponds to
  // |v.id.ToString()|, the return value to |v.orbit.vel|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept. Must
  // be called after initialization.
  virtual Velocity<AliceSun> VesselParentRelativeVelocity(
      GUID const& vessel_guid) const;

  // Returns the position of the celestial at index |celestial_index| relative
  // to its parent at current time. For a KSP |CelestialBody| |b|, the argument
  // corresponds to |b.flightGlobalsIndex|, the return value to |b.orbit.pos|.
  // A celestial with index |celestial_index| must have been inserted, and it
  // must not be the sun. Must be called after initialization.
  virtual Displacement<AliceSun> CelestialDisplacementFromParent(
      Index const celestial_index) const;

  // Returns the velocity of the celestial at index |celestial_index| relative
  // to its parent at current time. For a KSP |CelestialBody| |b|, the argument
  // corresponds to |b.flightGlobalsIndex|, the return value to |b.orbit.vel|.
  // A celestial with index |celestial_index| must have been inserted, and it
  // must not be the sun. Must be called after initialization.
  virtual Velocity<AliceSun> CelestialParentRelativeVelocity(
      Index const celestial_index) const;

 private:
  // Represents a KSP |CelestialBody|.
  struct Celestial {
    Celestial() = delete;
    Celestial(Celestial const&) = delete;
    Celestial(Celestial&&) = delete;
    ~Celestial() = default;

    explicit Celestial(GravitationalParameter const& gravitational_parameter)
      : body(new Body(gravitational_parameter)) {}
    std::unique_ptr<Body const> const body;
    // The parent body for the 2-body approximation. Not owning, must only
    // be null for the sun.
    Celestial const* parent = nullptr;
    // The past and present trajectory of the body. It ends at |HistoryTime()|.
    std::unique_ptr<Trajectory<Barycentre>> history;
    // A child trajectory of |*history|. It is forked at |history->last_time()|
    // and continues it until |current_time_|. It is computed with a
    // non-constant timestep, which breaks symplecticity. |history| is advanced
    // with a constant timestep as soon as possible, and |prolongation| is then
    // restarted from this new end of |history|.
    // Not owning, not null.
    Trajectory<Barycentre>* prolongation = nullptr;
  };

  // Represents a KSP |Vessel|.
  struct Vessel {
    Vessel() = delete;
    Vessel(Vessel const&) = delete;
    Vessel(Vessel&&) = delete;
    ~Vessel() = default;

    // Constructs a vessel whose parent is initially |*parent|. |parent| must
    // not be null. No transfer of ownership.
    explicit Vessel(Celestial const* parent) : parent(CHECK_NOTNULL(parent)) {}
    // A massless |Body|.
    std::unique_ptr<Body const> const body = new Body(GravitationalParameter());
    // The parent body for the 2-body approximation. Not owning, must not be
    // null.
    Celestial const* parent;
    // The past and present trajectory of the body. It ends at |HistoryTime()|
    // unless |*this| was created after |HistoryTime()|, in which case it ends
    // at |current_time_|.
    std::unique_ptr<Trajectory<Barycentre>> history;
    // A child trajectory of |*history|. It is forked at |history->last_time()|
    // and continues it until |current_time_|. It is computed with a
    // non-constant timestep, which breaks symplecticity. |history| is advanced
    // with a constant timestep as soon as possible, and |prolongation| is then
    // restarted from this new end of |history|.
    // Not owning, is null when the vessel is added and becomes non-null when
    // |history| is next advanced for all vessels and celestials. In the
    // meantime, |history| is advanced with small, non-constant timesteps to
    // catch up with the synchronous constant-timestep integration.
    // |this| is in |new_vessels_| if and only if |prolongation| is null.
    Trajectory<Barycentre>* prolongation = nullptr;
    // Whether to keep the |Vessel| during the next call to |AdvanceTime|.
    bool keep = true;
  };

  using GUIDToOwnedVessel = std::map<GUID, std::unique_ptr<Vessel>>;
  using GUIDToUnownedVessel = std::map<GUID, Vessel* const>;

  // The common last time of the histories of synchronized vessels and
  // celestials.
  Instant const& HistoryTime() const;

  // The rotation between the |World| basis at |current_time_| and the
  // |Barycentre| axes. Since |WorldSun| is not a rotating reference frame,
  // this change of basis is all that's required to convert relative velocities
  // or displacements between simultaneous events.
  Rotation<Barycentre, WorldSun> PlanetariumRotation() const;

  // Utilities for |AdvanceTime|.

  // Remove vessels for which |keep| is false, and sets |keep| to false for the
  // remaining ones.
  void CleanUpVessels();

  // Given |vessel| and an iterator to it in |new_vessels_|, checks that it has
  // been given an initial state, i.e. that its |history| is not null, and that
  // the following are equivalent:
  // * |vessel| is in |new_vessels_|
  // * its |prolongation| is null
  // * its |history->last_time()| is greater than |HistoryTime()|.
  // Also checks that |history->last_time()| is at least |HistoryTime()|.
  void CheckVesselInvariants(
      Vessel const& vessel,
      GUIDToUnownedVessel::iterator const it_in_new_vessels) const;

  // Evolves the histories of the |celestials_| and of the synchronized vessels
  // up to at most |t|. |t| must be large enough that at least one step of
  // size |kΔt| can fit between |current_time_| and |t|.
  void EvolveSynchronizedHistories(Instant const& t);

  // Synchronizes the |new_vessels_| and clears |new_vessels_|.
  void SynchronizeNewHistories();

  // Resets the prolongations of all vessels and celestials to |HistoryTime()|.
  // All vessels and celestials must have a null |prolongation|.
  void ResetProlongations();

  // Evolves the prolongations of all celestials and synchronized vessels, as
  // well as the histories of unsynchronized vessels, up to exactly instant |t|.
  void EvolveProlongationsAndUnsynchronizedHistories(Instant const& t);

  // TODO(egg): Constant time step for now.
  Time const kΔt = 10 * Second;

  GUIDToOwnedVessel vessels_;
  std::map<Index, std::unique_ptr<Celestial>> celestials_;

  // Vessels which have been recently inserted after |HistoryTime()|. For these
  // vessels, |history->last_time > HistoryTime()|. They have a null
  // |prolongation|. The pointers are not owning and not null.
  std::map<GUID, Vessel* const> new_vessels_;

  NBodySystem<Barycentre> solar_system_;
  // The symplectic integrator computing the synchronized histories.
  SPRKIntegrator<Length, Speed> history_integrator_;
  // The integrator computing the prolongations and the histories before they
  // are synchronized.
  SPRKIntegrator<Length, Speed> prolongation_integrator_;

  // Whether initialization is ongoing.
  Monostable initializing;

  Angle planetarium_rotation_;
  // The current in-game universal time.
  Instant current_time_;
  Celestial* sun_;  // Not owning, not null.

  // We want to check |HistoryTime|.
  FRIEND_TEST(PluginTest, AdvanceTime);
};

}  // namespace ksp_plugin
}  // namespace principia
