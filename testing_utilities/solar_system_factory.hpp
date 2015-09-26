#pragma once

#include <memory>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using astronomy::ICRFJ2000Equator;
using base::not_null;
using physics::SolarSystem;

namespace testing_utilities {

// A helper class for constructing physics::SolarSystem objects for testing.
class SolarSystemFactory {
 public:
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

  // Specifies the accuracy of the modeling.
  enum class Accuracy {
    // All the bodies with a mass greater than or equal to that of Pluto,
    // modelled as point masses.
    kMajorBodiesOnly,
    // Same as above, with some smaller satellites of the main planets.
    kMinorAndMajorBodies,
    // Same as above, with oblateness for the gas giants.
    kAllBodiesAndOblateness,
  };

  // A solar system at the time of the launch of Простейший Спутник-1.
  static not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>>
  AtСпутник1Launch(Accuracy const accuracy);

  // A solar system at the time of the launch of Простейший Спутник-2.
  static not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>>
  AtСпутник2Launch(Accuracy const accuracy);

  // Returns the index of the parent of the body with the given |index|.
  // Because enums are broken in C++ we use ints.  Sigh.
  static int parent(int const index);
  // The name of the body with the given |index|.
  static std::string name(int const index);
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_factory_body.hpp"
