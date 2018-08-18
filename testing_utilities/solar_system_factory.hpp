
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
namespace testing_utilities {
namespace internal_solar_system_factory {

using astronomy::ICRS;
using base::not_constructible;
using base::not_null;
using physics::SolarSystem;

// A helper class for constructing physics::SolarSystem objects for testing.
// TODO(egg): should this be a namespace instead?  It contains only static
// things, and there's no reason for it to be a type, it's not used in fancy
// template things.
class SolarSystemFactory : not_constructible {
 public:
  // The bodies are in decreasing order of mass.
  enum Index : int {
    Sun = 0,
    Jupiter = 1,
    Saturn = 2,
    Neptune = 3,
    Uranus = 4,
    Earth = 5,
    Venus = 6,
    Mars = 7,
    Mercury = 8,
    Ganymede = 9,
    Titan = 10,
    Callisto = 11,
    Io = 12,
    Moon = 13,
    Europa = 14,
    Triton = 15,
    Eris = 16,
    Pluto = 17,
    LastMajorBody = 17,
    Titania = 18,
    Oberon = 19,
    Rhea = 20,
    Iapetus = 21,
    Charon = 22,
    Ariel = 23,
    Umbriel = 24,
    Dione = 25,
    Ceres = 26,
    Tethys = 27,
    Vesta = 28,
    Enceladus = 29,
    Miranda = 30,
    Mimas = 31,
    Phobos = 32,
    Deimos = 33,
    LastBody = 33,
  };

  // Specifies the accuracy of the modeling.
  enum class Accuracy {
    // All the bodies with a mass greater than or equal to that of Pluto,
    // modelled as point masses.
    MajorBodiesOnly,
    // Same as above, with some smaller satellites of the main planets.
    MinorAndMajorBodies,
    // Same as above, with oblateness.
    AllBodiesAndOblateness,
  };

  template<typename Frame>
  static void AdjustAccuracy(Accuracy const accuracy,
                             SolarSystem<Frame>& solar_system);

  // A solar system at the time of the launch of Простейший Спутник-1.
  static not_null<std::unique_ptr<SolarSystem<ICRS>>>
  AtСпутник1Launch(Accuracy accuracy);

  // A solar system at the time of the launch of Простейший Спутник-2.
  static not_null<std::unique_ptr<SolarSystem<ICRS>>>
  AtСпутник2Launch(Accuracy accuracy);

  // Returns the index of the parent of the body with the given |index|.
  // Because enums are broken in C++ we use ints.  Sigh.
  static int parent(int index);
  // The name of the body with the given |index|.
  static std::string name(int index);
};

}  // namespace internal_solar_system_factory

using internal_solar_system_factory::SolarSystemFactory;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_factory_body.hpp"
