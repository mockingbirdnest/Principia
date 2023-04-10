#pragma once

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// This namespace contains the other non-SI units listed in the BIPM's
// SI brochure 8, section 4.1, table 8,
// http://www.bipm.org/en/si/si_brochure/chapter4/table8.html.
namespace _bipm {
namespace internal {

constexpr Pressure Bar                 = 1e5 * Pascal;
constexpr Pressure MillimetreOfMercury = 133.322 * Pascal;
constexpr Length   Ångström            = 1e-10 * Metre;
constexpr Length   NauticalMile        = 1852 * Metre;
constexpr Speed    Knot                = 1 * NauticalMile / Hour;
constexpr Area     Barn                = 1e-28 * Pow<2>(Metre);

}  // namespace internal

using internal::Ångström;
using internal::Bar;
using internal::Barn;
using internal::Knot;
using internal::MillimetreOfMercury;
using internal::NauticalMile;

}  // namespace _bipm
}  // namespace quantities
}  // namespace principia
