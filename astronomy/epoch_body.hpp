
#pragma once

#include "astronomy/epoch.hpp"

namespace principia {

using quantities::si::Day;

namespace astronomy {

constexpr Instant JulianDate(double const days) {
  return J2000 + (days - 2451545.0) * Day;
}

constexpr Instant ModifiedJulianDate(double const days) {
  return J2000 + (days - 51544.5) * Day;
}

}  // namespace astronomy
}  // namespace principia
