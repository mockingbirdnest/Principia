
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
    sun = 0,
    jupiter = 1,
    saturn = 2,
    neptune = 3,
    uranus = 4,
    earth = 5,
    venus = 6,
    mars = 7,
    mercury = 8,
    ganymede = 9,
    titan = 10,
    callisto = 11,
    io = 12,
    moon = 13,
    europa = 14,
    triton = 15,
    eris = 16,
    pluto = 17,
    last_major_body = 17,
    titania = 18,
    oberon = 19,
    rhea = 20,
    iapetus = 21,
    charon = 22,
    ariel = 23,
    umbriel = 24,
    dione = 25,
    tethys = 26,
    last_body = 26,
  };

  // Specifies the accuracy of the modeling.
  enum class Accuracy {
    // All the bodies with a mass greater than or equal to that of Pluto,
    // modelled as point masses.
    major_bodies_only,
    // Same as above, with some smaller satellites of the main planets.
    minor_and_major_bodies,
    // Same as above, with oblateness for the gas giants.
    all_bodies_and_oblateness,
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
