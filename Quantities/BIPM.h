#pragma once

#include "Quantities.h"
#include "NamedQuantities.h"
#include "Numbers.h"
#include "SI.h"

namespace Principia {
// This namespace contains the other non-SI units listed in the BIPM's 
// SI brochure 8, section 4.1, table 8,
// http://www.bipm.org/en/si/si_brochure/chapter4/table8.html.
namespace BIPM {
Quantities::Pressure Bar                 = 1e5 * SI::Pascal;
Quantities::Pressure MillimetreOfMercury = 133.322 * SI::Pascal;
Quantities::Length   Ångström            = 1e-10 * SI::Metre;
Quantities::Length   NauticalMile        = 1852 * SI::Metre;
Quantities::Area     Barn                = 100 * SI::Femto(SI::Metre).Pow<2>();
Quantities::Speed    Knot                = 1 * NauticalMile / Hour;
}
}
