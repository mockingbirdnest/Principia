#pragma once

#include <memory>
#include <string>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/massive_body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using base::not_null;
using geometry::Frame;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using physics::Trajectory;

namespace testing_utilities {

// A reference frame with a basis.
// The frame is the International Celestial Reference Frame.
// The basis is defined from the orbit of the Earth at J2000.0 as follows:
// The xy plane is the plane of the Earth's orbit at J2000.0.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at J2000.0.
// The z axis is perpendicular to the xy-plane in the directional (+ or -) sense
// of Earth's north pole at J2000.0.
// The basis is right-handed and orthonormal.
using ICRFJ2000Ecliptic =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::ICRF_J2000_ECLIPTIC, true>;

// The xy plane is the plane of the Earth's mean equator at J2000.0.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at J2000.0.
// The z axis is along the Earth's mean north pole at J2000.0.
// The basis is right-handed and orthonormal.
// Note that |ICRFJ2000Equator| and |ICRFJ2000Ecliptic| share their x axis.
using ICRFJ2000Equator =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::ICRF_J2000_EQUATOR, true>;

// Rotation around the common x axis mapping equatorial coordinates to ecliptic
// coordinates.  The angle is the one defined by the XVIth General Assembly of
// the International Astronomical Union.
geometry::Rotation<ICRFJ2000Equator, ICRFJ2000Ecliptic> const
    kEquatorialToEcliptic =
        geometry::Rotation<ICRFJ2000Equator, ICRFJ2000Ecliptic>(
            23 * si::Degree + 26 * si::ArcMinute + 21.448 * si::ArcSecond,
            geometry::Bivector<double, ICRFJ2000Equator>({-1, 0, 0}));

geometry::Position<ICRFJ2000Ecliptic> const kSolarSystemBarycentre;

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
  std::vector<not_null<Trajectory<ICRFJ2000Ecliptic>*>> trajectories() const;

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
  std::vector<not_null<std::unique_ptr<physics::Trajectory<ICRFJ2000Ecliptic>>>>
      trajectories_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/solar_system_body.hpp"
