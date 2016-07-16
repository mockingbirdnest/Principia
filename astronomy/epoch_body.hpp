
#pragma once

#include "astronomy/epoch.hpp"
#include "quantities/si.hpp"

namespace principia {

using quantities::si::Day;

namespace astronomy {

constexpr Instant jd0  = J2000 - 2451545.0 * Day;
constexpr Instant mjd0 = J2000 - 51544.5 * Day;

constexpr Instant JulianDate(double const days) {
  return jd0 + days * Day;
}

constexpr Instant ModifiedJulianDate(double const days) {
  return mjd0 + days * Day;
}

}  // namespace astronomy
}  // namespace principia
