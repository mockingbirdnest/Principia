#pragma once

#include "geometry/epoch.hpp"

using principia::si::Day;

namespace principia {
namespace geometry {

Instant const kJD0  = kJ2000 - 2451545.0 * Day;
Instant const kMJD0 = kJ2000 - 51544.5 * Day;

Instant JulianDate(double const days) {
  return kJD0 + days * Day;
}

Instant ModifiedJulianDate(double const days) {
  return kMJD0 + days * Day;
}

}  // namespace geometry
}  // namespace principia
