#include "numerics/scale_b.hpp"

namespace principia {
namespace numerics {

static_assert(ScaleB(3.0, 2) == 0x3p2);
static_assert(ScaleB(5.0, -2) == 0x5p-2);
static_assert(ScaleB(7.0, 0) == 0x7p0);

static_assert(ScaleB(0b1100 / 0x1p3, -1074) == 0x2p-1074);
// If subnormals are scaled by successive division by two, this assertion fails
// due to double rounding.
static_assert(ScaleB(0b1011 / 0x1p3, -1074) == 0x1p-1074);
static_assert(ScaleB(0b1010 / 0x1p3, -1074) == 0x1p-1074);

// Rescale the largest representable power of two to the smallest representable
// positive power of two.  The ratio between the two is not representable, so
// this test would for an implementation that attempts to represent bᴺ while
// computing x * bᴺ.
static_assert(ScaleB(0x1p1023, -1023 - 1074) == 0x1p-1074);

}  // namespace numerics
}  // namespace principia
