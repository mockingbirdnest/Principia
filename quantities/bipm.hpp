#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains the other non-SI units listed in the BIPM's
// SI brochure 8, section 4.1, table 8,
// http://www.bipm.org/en/si/si_brochure/chapter4/table8.html.
namespace bipm {

Pressure const Bar                 = 1e5 * si::Pascal;
Pressure const MillimetreOfMercury = 133.322 * si::Pascal;
Length   const Ångström            = 1e-10 * si::Metre;
Length   const NauticalMile        = 1852 * si::Metre;
Speed    const Knot                = 1 * NauticalMile / si::Hour;
Area     const Barn                = 1e-28 * Pow<2>(si::Metre);

}  // namespace bipm
}  // namespace quantities
}  // namespace principia
