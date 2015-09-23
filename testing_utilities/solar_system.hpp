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
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using astronomy::ICRFJ2000Ecliptic;
using astronomy::ICRFJ2000Equator;
using astronomy::kEquatorialToEcliptic;
using astronomy::kSolarSystemBarycentre;
using base::not_null;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::MassiveBody;
using quantities::si::ArcMinute;
using quantities::si::ArcSecond;
using quantities::si::Degree;

namespace testing_utilities {

class SolarSystem {
 public:
  using Bodies =
      // NOTE(Norgg) TODO(Egg) Another instance of removing a const from vector
      // here.
      std::vector<not_null<std::unique_ptr<physics::MassiveBody const>>>;

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

  // Factory.  The caller gets ownership of the pointers.
  // A solar system at the time of the launch of Простейший Спутник-1,
  // 1957-10-04T19:28:34Z (JD2436116.31150).
  static not_null<std::unique_ptr<SolarSystem>> AtСпутник1Launch(
      Accuracy const accuracy);

  // Factory.  The caller gets ownership of the pointers.
  // A solar system at the time of the launch of Простейший Спутник-2,
  // 1957-11-03T02:30:00Z (JD 2436145.60417)
  static not_null<std::unique_ptr<SolarSystem>> AtСпутник2Launch(
      Accuracy const accuracy);

  ~SolarSystem() = default;

  // The caller gets ownership of the bodies.  This function should only be
  // called once.
  Bodies massive_bodies();

  // This class retains ownership of the trajectories.
  std::vector<not_null<DiscreteTrajectory<ICRFJ2000Ecliptic>*>>
  trajectories() const;

  std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> initial_state() const;

  Instant const& time() const;

  // Returns the index of the parent of the body with the given |index|.
  // Because enums are broken in C++ we use ints.  Sigh.
  static int parent(int const index);
  // The name of the body with the given |index|.
  static std::string name(int const index);

 private:
  // A system containing the bodies appropriate for the given |accuracy|.
  explicit SolarSystem(Accuracy const accuracy);

  Bodies massive_bodies_;
  std::vector<not_null<std::unique_ptr<
      physics::DiscreteTrajectory<ICRFJ2000Ecliptic>>>> trajectories_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_body.hpp"
