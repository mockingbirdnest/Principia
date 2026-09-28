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

}  // namespace internal

using internal::ЛидовGraph;
using internal::ЛидовGrid;

}  // namespace _лидов
}  // namespace astronomy
}  // namespace principia
