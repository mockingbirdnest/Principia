
#pragma once

#include "geometry/epoch.hpp"

namespace principia {

using quantities::si::Day;

namespace geometry {

Instant const jd0  = J2000 - 2451545.0 * Day;
Instant const mjd0 = J2000 - 51544.5 * Day;

inline Instant JulianDate(double const days) {
  return jd0 + days * Day;
}

inline Instant ModifiedJulianDate(double const days) {
  return mjd0 + days * Day;
}

}  // namespace geometry
}  // namespace principia
