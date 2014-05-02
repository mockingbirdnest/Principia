#pragma once

#include "Quantities.hpp"
#include "NamedQuantities.hpp"
#include "Numbers.hpp"
#include "SI.hpp"

namespace principia {
// This namespace contains the other non-SI units listed in the BIPM's 
// SI brochure 8, section 4.1, table 8,
// http://www.bipm.org/en/si/si_brochure/chapter4/table8.html.
namespace BIPM {
Quantities::Pressure const Bar                 = 1e5 * SI::Pascal;
Quantities::Pressure const MillimetreOfMercury = 133.322 * SI::Pascal;
Quantities::Length   const Ångström            = 1e-10 * SI::Metre;
Quantities::Length   const NauticalMile        = 1852 * SI::Metre;
Quantities::Area     const Barn                = 1e-28 * SI::Metre.Pow<2>();
Quantities::Speed    const Knot                = 1 * NauticalMile / SI::Hour;
}  // namespace BIPM
}  // namespace principia
