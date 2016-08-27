
#pragma once

#include <functional>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/dynamic_frame.hpp"

namespace principia {

using geometry::Frame;
using geometry::Instant;
using geometry::Position;
using physics::DynamicFrame;

namespace ksp_plugin {

// Universal time 0, time of game creation.
// Putting the origin here makes the instants we use equal to the corresponding
// KSP universal time doubles.
Instant const UniversalTimeEpoch;

// Thanks to KSP's madness, the reference frame of the celestial body orbited by
// the active vessel, occasionally rotating with its surface, occasionally
// nonrotating.
// The basis is that of Unity's "world space" (this is a left-handed basis).
// The origin is the ineffable origin of Unity's "world space".
using World = Frame<serialization::Frame::PluginTag,
                    serialization::Frame::WORLD, false>;

// Same as |World| but with the y and z axes switched through the looking-glass:
// it is a right-handed basis. "We're all mad here. I'm mad. You're mad."
using AliceWorld = Frame<serialization::Frame::PluginTag,
                         serialization::Frame::ALICE_WORLD, false>;

// The barycentric reference frame of the solar system.
// The basis is the basis of |AliceWorld| at |UniversalTimeEpoch|.
using Barycentric = Frame<serialization::Frame::PluginTag,
                          serialization::Frame::BARYCENTRIC, true>;

// |Barycentric|, with its y and z axes swapped; the basis is left-handed.
using CelestialSphere = Frame<serialization::Frame::PluginTag,
                              serialization::Frame::CELESTIAL_SPHERE, true>;

// The surface frame of a celestial, with the x axis pointing to the origin of
// latitude and longitude, the y axis pointing to the pole with positive
// latitude, and the z axis oriented to form a left-handed basis.
using BodyWorld = Frame<serialization::Frame::PluginTag,
                        serialization::Frame::BODY_WORLD, false>;

// The frame used for trajectory plotting and manœuvre planning.  Its definition
// depends on the choice of a subclass of DynamicFrame.
using Navigation = Frame<serialization::Frame::PluginTag,
                         serialization::Frame::NAVIGATION, false>;

// A nonrotating referencence frame comoving with the sun with the same axes as
// |AliceWorld|. Since it is nonrotating (though not inertial), differences
// between velocities are consistent with those in an inertial reference frame.
// When |AliceWorld| rotates the axes are not fixed in the reference frame, so
// this (frame, basis) pair is inconsistent across instants. Operations should
// only be performed between simultaneous quantities, then converted to a
// consistent (frame, basis) pair before use.
using AliceSun = Frame<serialization::Frame::PluginTag,
                       serialization::Frame::ALICE_SUN, false>;

// Same as above, but with same axes as |World| instead of those of
// |AliceWorld|. The caveats are the same as for |AliceSun|.
using WorldSun = Frame<serialization::Frame::PluginTag,
                       serialization::Frame::WORLD_SUN, false>;

// Convenient instances of types from |physics| for the above frames.
using NavigationFrame = DynamicFrame<Barycentric, Navigation>;
using NavigationManœuvre = Manœuvre<Barycentric, Navigation>;

}  // namespace ksp_plugin
}  // namespace principia
