#pragma once

#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_recurrence {

using physics::RotatingBody;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Time;

// All references in this class are to Capderou (2012), Satellites : de Kepler
// au GPS; the notation follows that of the book.
//
// The ground track of a satellite recurs when, in a reference frame fixing the
// plane of the orbit and the equator of the Earth (such a reference frame is
// possible under the assumption that the orbital plane precesses around the
// equator), the following two conditions occur simultaneously:
// — the satellite has gone through an integer number of revolutions Nᴛₒ;
// — the primary has gone through an integer number of revolutions Cᴛₒ.
// Such a recurrence is characterized by the value Nᴛₒ / Cᴛₒ.
// This class represents that rational number; while it exposes many quantities
// that are relevant to mission design, they are all derived from it.
//
// In this class, the word “day” is used to denote a revolution of the primary
// in the reference frame fixing the orbital plane.  Note that this day depends
// on the nodal precession of the satellite: it is not any of the usual
// definitions of the day:
// — for a sun-synchronous orbit (see sections 7.5.4, 7.5.5), this day is the
//   mean solar day;
// — for an orbit with no nodal precession (i.e., a strictly polar orbit, see
//   section 7.1.3), this day is the stellar day.
// These days are counted negatively for a body with retrograde rotation.
class OrbitRecurrence final {
 public:
  // The following conditions must hold:
  //   Cᴛₒ ≠ 0;
  //   sign Cᴛₒ = sign νₒ if νₒ ≠ 0;
  //   |Dᴛₒ / Cᴛₒ| ≤ 1/2;
  //   gcd(Dᴛₒ, Cᴛₒ) = 1.
  OrbitRecurrence(int νₒ, int Dᴛₒ, int Cᴛₒ);

  // Returns the recurrence that mostly matches the given orbital
  // characteristics, limiting the value of |Cᴛₒ|.
  // The Nᴛₒ / Cᴛₒ of the result is the last convergent of the κ obtained from
  // the given arguments whose denominator is less than |max_abs_Cᴛₒ|.
  template<typename Frame>
  static OrbitRecurrence ClosestRecurrence(
      Time const& nodal_period,
      AngularFrequency const& nodal_precession,
      RotatingBody<Frame> const& primary,
      int max_abs_Cᴛₒ);

  // The Capderou recurrence triple [νₒ; Dᴛₒ; Cᴛₒ], see section 11.1.3.
  // This triple expresses the rational Nᴛₒ / Cᴛₒ as an integer and fractional
  // part, Nᴛₒ / Cᴛₒ = νₒ + Dᴛₒ/ Cᴛₒ.

  // While he does mention orbits around bodies whose rotation is retrograde
  // (Venus and Triton), Capderou does not cover ground track recurrence in that
  // case.  We derive the sign from equations (11.4) and (11.5): Nᴛₒ is always
  // positive, and Cᴛₒ is negative for a primary with retrograde rotation.

  // νₒ is the number of orbits per day, rounded to the nearest integer.
  // Note that νₒ < 0 if the rotation of the primary is retrograde.
  int νₒ() const;
  int Dᴛₒ() const;
  // Cᴛₒ is the length of the cycle in days.
  // Note that Cᴛₒ < 0 if the rotation of the primary is retrograde.
  int Cᴛₒ() const;

  // The number of nodal periods per cycle Nᴛₒ, see section 11.1.2.
  int number_of_revolutions() const;

  // The equatorial shift Δλᴇ, see section 8.3.2.
  // This is counted eastward; it is negative (westward shift) except when
  // orbiting bodies whose rotation is retrograde.
  Angle equatorial_shift() const;

  // The width δʀ of the base interval, see 11.5.2.  This is always positive (it
  // represents the spacing between consecutive ground tracks on the equator).
  // Note that δʀ = |Δλᴇ|.
  Angle base_interval() const;
  // The grid interval δ, see section 11.5.2.
  // This is always positive (it represents the closest spacing between ground
  // tracks on the equator).
  Angle grid_interval() const;

  // The subcycle Eᴛₒ*, see section 11.5.3.
  // After about Eᴛₒ* days (n = round(Nᴛₒ Eᴛₒ* / Cᴛₒ) revolutions), the ground
  // track passes δ away from the origin.
  // In terms of continued fractions, Eᴛₒ* is the denominator of the penultimate
  // convergent of Nᴛₒ / Cᴛₒ (it is the cycle that is closest to being exact
  // before Cᴛₒ).
  // Note that Eᴛₒ* < 0 if the rotation of the primary is retrograde.
  int subcycle() const;

 private:
  int νₒ_;
  int Dᴛₒ_;
  int Cᴛₒ_;
  int subcycle_;
};

std::ostream& operator<<(std::ostream& out, OrbitRecurrence const& recurrence);

}  // namespace internal_orbit_recurrence

using internal_orbit_recurrence::OrbitRecurrence;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_recurrence_body.hpp"
