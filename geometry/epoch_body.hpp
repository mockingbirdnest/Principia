#pragma once

#include "geometry/epoch.hpp"

namespace principia {

using quantities::si::Day;

namespace geometry {

Instant const kJD0  = kJ2000 - 2451545.0 * Day;
Instant const kMJD0 = kJ2000 - 51544.5 * Day;

inline Instant JulianDate(double const days) {
  return kJD0 + days * Day;
}

inline Instant ModifiedJulianDate(double const days) {
  return kMJD0 + days * Day;
}

}  // namespace geometry
}  // namespace principia
