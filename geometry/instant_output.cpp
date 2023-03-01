#include <limits>
#include <string>
#include <string_view>

#include "absl/strings/str_format.h"
#include "astronomy/time_scales.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _point {
namespace internal {

using astronomy::DateTimeAsTT;
using astronomy::J2000;
using astronomy::operator""_TT;
using astronomy::TTSecond;
using quantities::Time;
using quantities::si::Second;

std::ostream& operator<<(std::ostream& os, Instant const& t) {
  Time const from_j2000 = t - J2000;
  // Dates before JD0.5 and after the year 9999 are not supported; Sterbenz’s
  // lemma fails to apply in the second before J2000, we need to print the
  // sign of 0 for J2000 itself, and the second after J2000 has fractions of a
  // second that are too small to be reasonably printed in fixed format.
  if (t >= "JD0.5"_TT && t < "9999-12-31T24:00:00"_TT &&
      !(t > J2000 - 1 * Second && t < J2000 + 1 * Second)) {
    auto const tt_second = TTSecond(t);
    Instant const start_of_second = DateTimeAsTT(tt_second);
    // This subtraction is exact by Sterbenz’s lemma.
    Time const remainder = t - start_of_second;
    // |remainder| being the result of a subtraction of numbers greater than or
    // equal to 1, it is a multiple of 2u ≈ 2×10⁻¹⁶; 16 fractional decimal
    // digits suffice to unambiguously represent it (alternatively, as shown by
    // the static_assert, 17 decimal places are necessary, of which 16 are
    // fractional for numbers in [1, 10[; the integer part is taken care of by
    // |tt_second|).
    static_assert(std::numeric_limits<double>::max_digits10 - 1 == 16);
    return os << tt_second << ","
              << std::string_view(absl::StrFormat("%.16f", remainder / Second))
                     .substr(2)
              << " (TT)";
  }
  // The operator<< on the Time prints the requisite sign.
  return os << "J2000" << from_j2000 << " (TT)";
}

std::string DebugString(const Instant& t) {
  return (std::stringstream() << t).str();
}

}  // namespace internal
}  // namespace _point
}  // namespace geometry
}  // namespace principia
