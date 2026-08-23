#pragma once

#include <cstdint>

#include "absl/log/log.h"
#include "astronomy/orbital_elements.hpp"
#include "geometry/interval.hpp"
#include "graphics/colours.hpp"
#include "graphics/graph.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace astronomy {
namespace _лидов {
namespace internal {

using namespace principia::astronomy::_orbital_elements;
using namespace principia::geometry::_interval;
using namespace principia::graphics::_colours;
using namespace principia::graphics::_graph;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

enum class ЛидовGrid {
  None,
  MaxEccentricityMinInclination,
  MinEccentricityMaxInclination,
};

Graph<double, double> ЛидовGraph(OrbitalElements const& elements,
                                 std::int64_t width,
                                 std::int64_t height,
                                 RGBA32 background,
                                 RGB24 region_boundary_colour,
                                 RGB24 inclination_colour,
                                 RGB24 eccentricity_colour,
                                 RGB24 лидов_parameter_colour,
                                 ЛидовGrid grid);

// All functions in this file refer to an orbit perturbed as in the analysis of
// [Лид61].  The parameters c₁ and c₂ are as defined there.

// Returns c₁ such that an orbit with these values of c₁ and c₂ has no
// eccentricity-inclination exchange.
double ЛидовFrozenLine(double c₂);

// Returns c₁ such that the upper bound of eccentricity for an orbit with these
// values of c₁ and c₂ is e.
double ЛидовMaximalEccentricityLine(double e, double c₂);

// Returns the range of values of c₂ such that there exists a c₁ such that the
// upper bound of eccentricity for an orbit with these values of c₁ and c₂
// is e.
Interval<double> ЛидовMaximalEccentricityLineC₂Range(double e);

// Returns c₁ such that the upper bound of inclination for an orbit with these
// values of c₁ and c₂ is i.
double ЛидовMaximalInclinationLine(Angle i, double c₂);

// Returns the range of values of c₂ such that there exists a c₁ such that the
// upper bound of inclination for an orbit with these values of c₁ and c₂
// is i.
Interval<double> ЛидовMaximalInclinationLineC₂Range(Angle i);

// Returns c₁ such that the lower bound of inclination for an orbit with these
// values of c₁ and c₂ is i.
double ЛидовMinimalInclinationLine(Angle i, double c₂);

// Returns the range of values of c₂ such that there exists a c₁ such that the
// lower bound of inclination for an orbit with these values of c₁ and c₂
// is i.
Interval<double> ЛидовMinimalInclinationLineC₂Range(Angle cnst i);

// Returns the value of c₁ such that the lower bound of eccentricity for an
// orbit with these values of c₁ and c₂ is e.
double ЛидовMinimalEccentricityLeftLine(double e, double c₂);

// Returns the range of negative values of c₂ such that there exists a c₁ such
// that the lower bound of eccentricity for an orbit with these values of c₁
// and c₂ is e.
Interval<double> ЛидовMinimalEccentricityLeftLineC₂Range(double e);

// Returns the positive value of c₂ for which the lower bound of eccentricity is
// e.
double ЛидовMinimalEccentricityRightLineC₂(double e);

// Returns the maximal possible value of c₁ that can be attained when c₂ has the
// positive value for which the lower bound of eccentricity is e.
double ЛидовMinimalEccentricityRightLineC₁Max(double e);

}  // namespace internal

using internal::ЛидовFrozenLine;
using internal::ЛидовGraph;
using internal::ЛидовGrid;
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
