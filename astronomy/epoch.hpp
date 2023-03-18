#pragma once

#include "geometry/named_quantities.hpp"

// |geometry::Instant| represents instants of Terrestrial Time (TT).  The
// utilities in this file provide its standard epoch and two ways of specifying
// TT dates.

namespace principia {
namespace astronomy {
namespace _epoch {
namespace internal {

using namespace principia::geometry::_named_quantities;

// |J2000| represents to the standard epoch J2000.0.
// According to Resolution B1 (On the Use of Julian Dates) of the XXIIIrd IAU
// general assembly, "it is recommended that JD be specified as SI seconds in
// Terrestrial Time (TT)", see http://goo.gl/oPemRm. J2000.0 is by definition
// JD 2451545.0, i.e., noon on the first of January, 2000 (TT).
// "2000-01-01T12:00:00"_TT
// "2000-01-01T11:59:27,816"_TAI
// "2000-01-01T11:58:55,816"_UTC
constexpr Instant J2000;

}  // namespace internal

using internal::J2000;

}  // namespace _epoch
}  // namespace astronomy
}  // namespace principia

namespace principia::astronomy {
using namespace principia::astronomy::_epoch;
}  // namespace principia::astronomy
