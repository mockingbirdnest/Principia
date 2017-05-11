
#pragma once

#include "astronomy/epoch.hpp"

namespace principia {
namespace astronomy {
namespace internal_epoch {

using quantities::si::Day;


constexpr Instant JulianDate(double const days) {
  return J2000 + (days - 2451545.0) * Day;
}

constexpr double JulianDayNumber(Instant const& t) {
  return 2451545.0 + (t - J2000) / Day;
}

constexpr Instant ModifiedJulianDate(double const days) {
  return J2000 + (days - 51544.5) * Day;
}

constexpr double ModifiedJulianDayNumber(Instant const& t) {
  return 51544.5 + (t - J2000) / Day;
}

}  // namespace internal_epoch
}  // namespace astronomy
}  // namespace principia
