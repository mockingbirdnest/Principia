#pragma once

#include <vector>

#include "physics/n_body_system.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {

// ICRF/J2000.0 reference frame.
// The reference epoch is J2000.0.
// The xy plane is the plane of the Earth's orbit at the reference epoch.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at the reference epoch.
// The z axis is perpendicular to the xy-plane in the directional (+ or -) sense
// of Earth's north pole at the reference epoch.
// The reference frame is direct.
struct ICRFJ2000EclipticFrame;

// A system containing the 18 largest solar system bodies (Pluto and all larger
// bodies) at the time of the launch of Простейший Спутник-1,
// 1957-10-04T19:28:34Z (JD2436116.3115).
// The bodies are in decreasing order of mass,
//  0. Sun,
//  1. Jupiter,
//  2. Saturn,
//  3. Neptune,
//  4. Uranus,
//  5. Earth,
//  6. Venus,
//  7. Mars,
//  8. Mercury,
//  9. Ganymede,
// 10. Titan,
// 11. Callisto,
// 12. Io,
// 13. Moon,
// 14. Europa,
// 15. Triton,
// 16. Eris,
// 17. Pluto.
std::unique_ptr<physics::NBodySystem<ICRFJ2000EclipticFrame>>
SolarSystemAtСпутникLaunch();

// Number of days since the JD epoch. JD2436116.3115 is the time of the launch
// of Простейший Спутник-1.
quantities::Time const kСпутникLaunchDate = 2436116.3115 * si::Day;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_body.hpp"
