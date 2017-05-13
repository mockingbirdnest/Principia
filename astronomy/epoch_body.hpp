
#pragma once

#include "astronomy/epoch.hpp"

namespace principia {
namespace astronomy {
namespace internal_epoch {

using quantities::si::Day;

constexpr Instant JulianDate(double const days) {
  return J2000 + (days - 2451545.0) * Day;
}

constexpr Instant ModifiedJulianDate(double const days) {
  return J2000 + (days - 51544.5) * Day;
}

}  // namespace internal_epoch
}  // namespace astronomy
}  // namespace principia
