#pragma once

#include <vector>

#include "physics/n_body_system.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {

// A reference frame with a basis.
// The frame is the International Celestial Reference Frame.
// The basis is defined from the orbit of the Earth at J2000.0 as follows:
// The xy plane is the plane of the Earth's orbit at J2000.0.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at J2000.0.
// The z axis is perpendicular to the xy-plane in the directional (+ or -) sense
// of Earth's north pole at J2000.0.
// The basis is direct and orthonormal.
struct ICRFJ2000Ecliptic;

class SolarSystem {
 public:
  using Bodies = std::vector<std::unique_ptr<physics::Body const> const>;

  // Factory.  The caller gets ownership of the pointers.
  // A solar system at the time of the launch of Простейший Спутник-1,
  // 1957-10-04T19:28:34Z (JD2436116.3115).
  static std::unique_ptr<SolarSystem> AtСпутникLaunch();

  ~SolarSystem() = default;

  // The caller gets ownership of the bodies.  These functions should only be
  // called once.
  Bodies massive_bodies();
  Bodies massless_bodies();

  // This class retains ownership of the trajectories.
  physics::NBodySystem<ICRFJ2000Ecliptic>::Trajectories trajectories() const;

 private:
  // A system containing the 18 largest solar system bodies (Pluto and all
  // larger bodies)  The bodies are in decreasing order of mass,
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
  SolarSystem();

  Bodies massive_bodies_;
  Bodies massless_bodies_;
  std::vector<std::unique_ptr<physics::Trajectory<ICRFJ2000Ecliptic>>>
      trajectories_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_body.hpp"
