#pragma once

#include "physics/n_body_system.hpp"

namespace principia {
namespace testing_utilities {

// ICRF/J2000.0 reference frame.
// The reference epoch is J2000.0
// The xy plane is the plane of the Earth's orbit at the reference epoch.
// The x axis is out along ascending node of instantaneous plane of the Earth's
// orbit and the Earth's mean equator at the reference epoch.
// The z axis is perpendicular to the xy-plane in the directional (+ or -) sense
// of Earth's north pole at the reference epoch.
// The reference frame is direct.
struct ICRFJ2000EclipticFrame;

// A system containing the 18 largest solar system bodies (Pluto and all larger
// bodies) at the time of the launch of Простейший Спутник-1,
// 1957-10-04T19:28:34Z (JD2436116.3115).
physics::NBodySystem<ICRFJ2000EclipticFrame> SolarSystemAtSputnikLaunch();

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_body.hpp"
