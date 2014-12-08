#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"
#include "physics/frame.hpp"

using principia::geometry::Instant;
using principia::geometry::Position;
using principia::physics::Frame;

namespace principia {
namespace ksp_plugin {

enum Tag {
  kAliceSun,
  kAliceWorld,
  kBarycentric,
  kRendering,
  kWorld,
  kWorldSun,
};

// Universal time 0, time of game creation.
// Putting the origin here makes the instants we use equal to the corresponding
// KSP universal time doubles.
Instant const kUniversalTimeEpoch;

// Thanks to KSP's madness, the reference frame of the celestial body orbited by
// the active vessel, occasionally rotating with its surface, occasionally
// nonrotating.
// The basis is that of Unity's "world space" (this is a left-handed basis).
// The origin is the ineffable origin of Unity's "world space".
using World = Frame<Tag, Tag::kWorld, false>;

// Same as |World| but with the y and z axes switched through the looking-glass:
// it is a right-handed basis. "We're all mad here. I'm mad. You're mad."
using AliceWorld = Frame<Tag, Tag::kAliceWorld, false>;

// The barycentric reference frame of the solar system.
// The basis is the basis of |World| at |kUniversalTimeEpoch|.
// TODO(egg): it *should* be the barycentric frame. For the moment we're using
// the velocity of the sun at the time of construction as our reference.
// The origin is the position of the sun at the instant |initial_time| passed at
// construction.
using Barycentric = Frame<Tag, Tag::kBarycentric, true>;

// The frame used for rendering.  Its definition depends on the actual factory
// function used to create it, see class Transforms.
using Rendering = Frame<Tag, Tag::kRendering, false>;

// A nonrotating referencence frame comoving with the sun with the same axes as
// |AliceWorld|. Since it is nonrotating (though not inertial), differences
// between velocities are consistent with those in an inertial reference frame.
// When |AliceWorld| rotates the axes are not fixed in the reference frame, so
// this (frame, basis) pair is inconsistent across instants. Operations should
// only be performed between simultaneous quantities, then converted to a
// consistent (frame, basis) pair before use.
using AliceSun = Frame<Tag, Tag::kAliceSun, false>;

// Same as above, but with same axes as |World| instead of those of
// |AliceWorld|. The caveats are the same as for |AliceSun|.
using WorldSun = Frame<Tag, Tag::kWorldSun, false>;

}  // namespace ksp_plugin
}  // namespace principia
