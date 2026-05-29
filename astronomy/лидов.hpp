#pragma once

#include "absl/log/log.h"
#include "geometry/interval.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace astronomy {
namespace _лидов {
namespace internal {

using namespace principia::geometry::_interval;
using namespace principia::journal::_method;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// All functions in this file refer to an orbit perturbed as in the analysis of
// [Лид61].  The parameters c₁ and c₂ are as defined there.

static inline const Angle i_critical = ArcCos(Sqrt(3.0 / 5.0));

// Returns c₁ such that an orbit with these values of c₁ and c₂ has no
// eccentricity-inclination exchange.
inline double ЛидовFrozenLine(double const c₂) {
  CHECK_LE(c₂, 0);
  return 3.0 / 5.0 - 2 * Sqrt(-3.0 / 5.0 * c₂) - c₂;
}

// Returns c₁ such that the upper bound of eccentricity for an orbit with these
// values of c₁ and c₂ is e.
inline double ЛидовMaximalEccentricityLine(double const e, double const c₂) {
  double const e² = Pow<2>(e);
  return 3.0 / 5.0 - c₂ + c₂ / e² - 3 * e² / 5.0;
}

// Returns the range of values of c₂ such that there exists a c₁ such that the
// upper bound of eccentricity for an orbit with these values of c₁ and c₂
// is e.
Interval<double> ЛидовMaximalEccentricityLineC₂Range(double const e) {
  double const e² = Pow<2>(e);
  double const e⁴ = Pow<4>(e);
  return {-3.0 * e⁴ / 5.0, 2.0 * e² / 5.0};
}

// Returns c₁ such that the upper bound of inclination for an orbit with these
// values of c₁ and c₂ is i.
inline double ЛидовMaximalInclinationLine(Angle const i, double const c₂) {
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return c₂ < 0
             ? cos²_i * (5.0 * cos²_i - 5.0 * c₂ - 3.0) / (5.0 * cos²_i - 3.0)
             : (2.0 - 5.0 * c₂) * cos²_i / 2.0;
}

// Returns the range of values of c₂ such that there exists a c₁ such that the
// upper bound of inclination for an orbit with these values of c₁ and c₂
// is i.
Interval<double> ЛидовMaximalInclinationLineC₂Range(Angle const i) {
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return {i > i_critical ? -Pow<2>(1.0 - 5.0 * Cos(2.0 * i)) / 60.0 : 0,
          2.0 / 5.0};
}

// Returns c₁ such that the lower bound of inclination for an orbit with these
// values of c₁ and c₂ is i.
inline double ЛидовMinimalInclinationLine(Angle const i, double const c₂) {
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return cos²_i * (5.0 * cos²_i - 5.0 * c₂ - 3.0) / (5.0 * cos²_i - 3.0);
}

// Returns the range of values of c₂ such that there exists a c₁ such that the
// lower bound of inclination for an orbit with these values of c₁ and c₂
// is i.
Interval<double> ЛидовMinimalInclinationLineC₂Range(Angle const i) {
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return i > i_critical
             ? Interval<double>{cos²_i - 3.0 / 5.0,
                                -Pow<2>(1.0 - 5.0 * Cos(2 * i)) / 60.0}
             : Interval<double>{0, cos²_i - 3.0 / 5.0};
}

// Returns the value of c₁ such that the lower bound of eccentricity for an
// orbit with these values of c₁ and c₂ is e.
double ЛидовMinimalEccentricityLeftLine(double const e, double const c₂) {
  double const e² = Pow<2>(e);
  return 3.0 / 5.0 - c₂ + c₂ / e² - 3.0 * e² / 5.0;
}

// Returns the range of negative values of c₂ such that there exists a c₁ such
// that the lower bound of eccentricity for an orbit with these values of c₁
// and c₂ is e.
Interval<double> ЛидовMinimalEccentricityLeftLineC₂Range(double const e) {
  double const e² = Pow<2>(e);
  double const e⁴ = Pow<4>(e);
  return {-3.0 * e² / 5.0, -3.0 * e⁴ / 5.0};
}

// Returns the positive value of c₂ for which the lower bound of eccentricity is
// e.
inline double ЛидовMinimalEccentricityRightLineC₂(double const e) {
  double const e² = Pow<2>(e);
  return 2.0 * e² / 5;
}

// Returns the maximal possible value of c₁ that can be attained when c₂ has the
// positive value for which the lower bound of eccentricity is e.
inline double ЛидовMinimalEccentricityRightLineC₁Max(double const e) {
  double const e² = Pow<2>(e);
  return 1.0 - e²;
}

}  // namespace internal

using internal::ЛидовFrozenLine;
using internal::ЛидовMaximalEccentricityLine;
using internal::ЛидовMaximalEccentricityLineC₂Range;
using internal::ЛидовMaximalInclinationLine;
using internal::ЛидовMaximalInclinationLineC₂Range;
using internal::ЛидовMinimalEccentricityLeftLine;
using internal::ЛидовMinimalEccentricityLeftLineC₂Range;
using internal::ЛидовMinimalEccentricityRightLineC₁Max;
using internal::ЛидовMinimalEccentricityRightLineC₂;
using internal::ЛидовMinimalInclinationLine;
using internal::ЛидовMinimalInclinationLineC₂Range;

}  // namespace _лидов
}  // namespace astronomy
}  // namespace principia
