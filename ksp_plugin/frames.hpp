
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

using geometry::Arbitrary;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::NonRotating;
using geometry::Permutation;
using physics::DynamicFrame;

// Thanks to KSP's madness, the reference frame of the celestial body orbited by
// the active vessel, occasionally rotating with its surface, occasionally
// nonrotating.
// The basis is that of Unity's "world space".  The origin is the ineffable
// origin of Unity's "world space".
using World = Frame<serialization::Frame::PluginTag,
                    Arbitrary,
                    Handedness::Left,
                    serialization::Frame::WORLD>;

// Same as |World| but with the y and z axes switched through the looking-glass:
// it is a right-handed basis. "We're all mad here. I'm mad. You're mad."
using AliceWorld = Frame<serialization::Frame::PluginTag,
                         Arbitrary,
                         Handedness::Right,
                         serialization::Frame::ALICE_WORLD>;

// The barycentric reference frame of the solar system.
using Barycentric = Frame<serialization::Frame::PluginTag,
                          Inertial,
                          Handedness::Right,
                          serialization::Frame::BARYCENTRIC>;

// The |Apparent...| frames are used for data obtained after the physics
// simulation of the game has run, and before we perform our correction.

// |World| coordinates from the game, but before the correction.
using ApparentWorld = Frame<serialization::Frame::PluginTag,
                            NonRotating,
                            Handedness::Right,
                            serialization::Frame::APPARENT_WORLD>;

// The axes are those of |Barycentric|.  The origin is that of |ApparentWorld|,
// and should not be depended upon.
using Apparent = Frame<serialization::Frame::PluginTag,
                       NonRotating,
                       Handedness::Right,
                       serialization::Frame::APPARENT>;

// |Barycentric|, with its y and z axes swapped.
using CelestialSphere = Frame<serialization::Frame::PluginTag,
                              Inertial,
                              Handedness::Left,
                              serialization::Frame::CELESTIAL_SPHERE>;

// The surface frame of a celestial, with the x axis pointing to the origin of
// latitude and longitude, the y axis pointing to the pole with positive
// latitude, and the z axis oriented to form a left-handed basis.
using BodyWorld = Frame<serialization::Frame::PluginTag,
                        Arbitrary,
                        Handedness::Left,
                        serialization::Frame::BODY_WORLD>;

// The frame used for the navball.  Its definition depends on the choice of a
// subclass of FrameField.
using Navball = Frame<serialization::Frame::PluginTag,
                      Arbitrary,
                      Handedness::Left,
                      serialization::Frame::NAVBALL>;

// The frame used for trajectory plotting and manœuvre planning.  Its definition
// depends on the choice of a subclass of DynamicFrame.
using Navigation = Frame<serialization::Frame::PluginTag,
                         Arbitrary,
                         Handedness::Right,
                         serialization::Frame::NAVIGATION>;

// The plotting frame, but with the y and z axes swapped compared to
// |Navigation|.  This frame defines the camera horizontal, and its angular
// velocity defines the angular velocity of the camera (note that the linear
// motion of the camera is defined in-game by following a specific target, which
// may be in motion with respect to |CameraReference|, so the camera is not
// necessarily at rest in that frame).
using CameraReference = Frame<serialization::Frame::PluginTag,
                              Arbitrary,
                              Handedness::Left,
                              serialization::Frame::CAMERA_REFERENCE>;

// A nonrotating referencence frame comoving with the sun with the same axes as
// |AliceWorld|. Since it is nonrotating (though not inertial), differences
// between velocities are consistent with those in an inertial reference frame.
// When |AliceWorld| rotates the axes are not fixed in the reference frame, so
// this (frame, basis) pair is inconsistent across instants. Operations should
// only be performed between simultaneous quantities, then converted to a
// consistent (frame, basis) pair before use.
using AliceSun = Frame<serialization::Frame::PluginTag,
                       Arbitrary,
                       Handedness::Right,
                       serialization::Frame::ALICE_SUN>;

// Same as above, but with same axes as |World| instead of those of
// |AliceWorld|. The caveats are the same as for |AliceSun|.
using WorldSun = Frame<serialization::Frame::PluginTag,
                       Arbitrary,
                       Handedness::Left,
                       serialization::Frame::WORLD_SUN>;

// Used to identify coordinates in the projective plane.  This is *not* the
// OpenGL camera (which is right-handed) but the Unity camera.
using Camera = Frame<serialization::Frame::PluginTag,
                     Arbitrary,
                     Handedness::Left,
                     serialization::Frame::CAMERA>;

// The frame that defines the orientation of a part.
using RigidPart = Frame<serialization::Frame::PluginTag,
                        Arbitrary,
                        Handedness::Left,
                        serialization::Frame::RIGID_PART>;

// The body-centred non-rotating frame for the current main body.
using MainBodyCentred = Frame<serialization::Frame::PluginTag,
                              NonRotating,
                              Handedness::Right,
                              serialization::Frame::MAIN_BODY_CENTRED>;

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
using internal_frames::CameraReference;
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
