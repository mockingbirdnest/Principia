
#pragma once

#include <functional>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/permutation.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/dynamic_frame.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_frames {

using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::NonInertial;
using geometry::Permutation;
using physics::DynamicFrame;

// Thanks to KSP's madness, the reference frame of the celestial body orbited by
// the active vessel, occasionally rotating with its surface, occasionally
// nonrotating.
// The basis is that of Unity's "world space" (this is a left-handed basis).
// The origin is the ineffable origin of Unity's "world space".
using World = Frame<serialization::Frame::PluginTag,
                    serialization::Frame::WORLD,
                    NonInertial,
                    Handedness::Left>;

// Same as |World| but with the y and z axes switched through the looking-glass:
// it is a right-handed basis. "We're all mad here. I'm mad. You're mad."
using AliceWorld = Frame<serialization::Frame::PluginTag,
                         serialization::Frame::ALICE_WORLD,
                         NonInertial>;

// The barycentric reference frame of the solar system.
using Barycentric = Frame<serialization::Frame::PluginTag,
                          serialization::Frame::BARYCENTRIC,
                          Inertial>;

// The axes are those of |Barycentric|.  The origin is that of |World|.  This
// frame is used for degrees of freedom obtained after the physics simulation of
// the game has run, and before we perform our correction: the origin has no
// physical significance.
using ApparentBubble = Frame<serialization::Frame::PluginTag,
                             serialization::Frame::APPARENT_BUBBLE,
                             NonInertial>;

// |Barycentric|, with its y and z axes swapped; the basis is left-handed.
using CelestialSphere = Frame<serialization::Frame::PluginTag,
                              serialization::Frame::CELESTIAL_SPHERE,
                              Inertial,
                              Handedness::Left>;

// The surface frame of a celestial, with the x axis pointing to the origin of
// latitude and longitude, the y axis pointing to the pole with positive
// latitude, and the z axis oriented to form a left-handed basis.
using BodyWorld = Frame<serialization::Frame::PluginTag,
                        serialization::Frame::BODY_WORLD,
                        NonInertial,
                        Handedness::Left>;

// The frame used for the navball.  Its definition depends on the choice of a
// subclass of FrameField.  This frame is left-handed.
using Navball = Frame<serialization::Frame::PluginTag,
                      serialization::Frame::NAVBALL,
                      NonInertial,
                      Handedness::Left>;

// The frame used for trajectory plotting and manœuvre planning.  Its definition
// depends on the choice of a subclass of DynamicFrame.
using Navigation = Frame<serialization::Frame::PluginTag,
                         serialization::Frame::NAVIGATION,
                         NonInertial>;

// A nonrotating referencence frame comoving with the sun with the same axes as
// |AliceWorld|. Since it is nonrotating (though not inertial), differences
// between velocities are consistent with those in an inertial reference frame.
// When |AliceWorld| rotates the axes are not fixed in the reference frame, so
// this (frame, basis) pair is inconsistent across instants. Operations should
// only be performed between simultaneous quantities, then converted to a
// consistent (frame, basis) pair before use.
using AliceSun = Frame<serialization::Frame::PluginTag,
                       serialization::Frame::ALICE_SUN,
                       NonInertial>;

// Same as above, but with same axes as |World| instead of those of
// |AliceWorld|. The caveats are the same as for |AliceSun|.
using WorldSun = Frame<serialization::Frame::PluginTag,
                       serialization::Frame::WORLD_SUN,
                       NonInertial,
                       Handedness::Left>;

// Used to identify coordinates in the projective plane.
using Camera = Frame<serialization::Frame::PluginTag,
                     serialization::Frame::CAMERA,
                     NonInertial>;

// The frame that defines the orientation of a part.
using RigidPart = Frame<serialization::Frame::PluginTag,
                        serialization::Frame::RIGID_PART,
                        NonInertial,
                        Handedness::Left>;

// The body-centred non-rotating frame for the current main body.
using MainBodyCentred = Frame<serialization::Frame::PluginTag,
                              serialization::Frame::MAIN_BODY_CENTRED,
                              NonInertial>;

// Convenient instances of types from |physics| for the above frames.
using NavigationFrame = DynamicFrame<Barycentric, Navigation>;
using NavigationManœuvre = Manœuvre<Barycentric, Navigation>;

// The map between the vector spaces of |WorldSun| and |AliceSun|.
Permutation<WorldSun, AliceSun> const sun_looking_glass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

}  // namespace internal_frames

using internal_frames::AliceSun;
using internal_frames::AliceWorld;
using internal_frames::ApparentBubble;
using internal_frames::Barycentric;
using internal_frames::BodyWorld;
using internal_frames::Camera;
using internal_frames::CelestialSphere;
using internal_frames::MainBodyCentred;
using internal_frames::Navball;
using internal_frames::Navigation;
using internal_frames::NavigationFrame;
using internal_frames::NavigationManœuvre;
using internal_frames::RigidPart;
using internal_frames::World;
using internal_frames::WorldSun;
using internal_frames::sun_looking_glass;

}  // namespace ksp_plugin
}  // namespace principia
