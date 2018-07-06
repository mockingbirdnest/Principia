
#pragma once

#include "base/macros.hpp"
#include <limits>

#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

// |geometry::Instant| represents instants of Terrestrial Time (TT).  The
// utilities in this file provide its standard epoch and two ways of specifying
// TT dates.

namespace principia {
namespace astronomy {
namespace internal_epoch {

using geometry::Instant;
using quantities::Infinity;
using quantities::Time;
using quantities::si::Second;

// |J2000| represents to the standard epoch J2000.0.
// According to Resolution B1 (On the Use of Julian Dates) of the XXIIIrd IAU
// general assembly, "it is recommended that JD be specified as SI seconds in
// Terrestrial Time (TT)", see http://goo.gl/oPemRm. J2000.0 is by definition
// JD 2451545.0, i.e., noon on the first of January, 2000 (TT).
// "2000-01-01T12:00:00"_TT
// "2000-01-01T11:59:27,816"_TAI
// "2000-01-01T11:58:55,816"_UTC
constexpr Instant J2000;

CONSTEXPR_INFINITY Instant InfinitePast = J2000 - Infinity<Time>();
CONSTEXPR_INFINITY Instant InfiniteFuture = J2000 + Infinity<Time>();

}  // namespace internal_epoch

using internal_epoch::InfiniteFuture;
using internal_epoch::InfinitePast;
using internal_epoch::J2000;

}  // namespace astronomy
}  // namespace principia
