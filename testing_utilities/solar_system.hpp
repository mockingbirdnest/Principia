#pragma once

#include <memory>
#include <string>
#include <vector>

#include "geometry/named_quantities.hpp"
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
geometry::Position<ICRFJ2000Ecliptic> const kSolarSystemBarycentre;

class SolarSystem {
 public:
  using Bodies = std::vector<std::unique_ptr<physics::Body const> const>;

  // The indices of the bodies in the |Bodies| vector.
  // The bodies are in decreasing order of mass.
  enum Index : int {
    kSun = 0,
    kJupiter = 1,
    kSaturn = 2,
    kNeptune = 3,
    kUranus = 4,
    kEarth = 5,
    kVenus = 6,
    kMars = 7,
    kMercury = 8,
    kGanymede = 9,
    kTitan = 10,
    kCallisto = 11,
    kIo = 12,
    kMoon = 13,
    kEuropa = 14,
    kTriton = 15,
    kEris = 16,
    kPluto = 17,
    kTitania = 18,
    kOberon = 19,
    kRhea = 20,
    kIapetus = 21,
    kCharon = 22,
    kAriel = 23,
    kUmbriel = 24,
    kDione = 25,
    kTethys = 26,
  };

  // Specifies the accuracy of the modeling.  |kMajorBodiesOnly| has all the
  // bodies with masses greater than or equal to Pluto, and all these bodies are
  // spheric.  |kMinorBodiesAndOblateness| adds some smaller satellites and the
  // oblateness.
  enum class Accuracy {
    kMajorBodiesOnly,
    kMinorBodiesAndOblateness,
  };

  // Factory.  The caller gets ownership of the pointers.
  // A solar system at the time of the launch of Простейший Спутник-1,
  // 1957-10-04T19:28:34Z (JD2436116.31150).
  static std::unique_ptr<SolarSystem> AtСпутник1Launch(Accuracy const accuracy);

  // Factory.  The caller gets ownership of the pointers.
  // A solar system at the time of the launch of Простейший Спутник-2,
  // 1957-11-03T02:30:00Z (JD 2436145.60417)
  static std::unique_ptr<SolarSystem> AtСпутник2Launch(Accuracy const accuracy);

  ~SolarSystem() = default;

  // The caller gets ownership of the bodies.  These functions should only be
  // called once.
  Bodies massive_bodies();
  Bodies massless_bodies();

  // This class retains ownership of the trajectories.
  physics::NBodySystem<ICRFJ2000Ecliptic>::Trajectories trajectories() const;

  // Returns the index of the parent of the body with the given |index|.
  // Because enums are broken in C++ we use ints.  Sigh.
  static int parent(int const index);
  // The name of the body with the given |index|.
  static std::string name(int const index);

 private:
  // A system containing the bodies listed in |Index|.
  SolarSystem();

  Bodies massive_bodies_;
  Bodies massless_bodies_;
  std::vector<std::unique_ptr<physics::Trajectory<ICRFJ2000Ecliptic>>>
      trajectories_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_body.hpp"
