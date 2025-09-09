#pragma once

#include "numerics/чебышёв_lobatto.hpp"

#include <array>
#include <cstdint>

#include "glog/logging.h"
#include "numerics/fma.hpp"
#include "numerics/sin_cos.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _чебышёв_lobatto {
namespace internal {

using namespace principia::numerics::_fma;
using namespace principia::numerics::_sin_cos;
using namespace principia::quantities::_si;

// Returns the Чебышёв-Lobatto point Cos(πk/N).  Should only be called with odd
// k (values for even k should be obtained by calling the N/2 instantiation).
template<int N>
double ЧебышёвLobattoPoint(std::int64_t const k) {
  static const std::array<double, N>* const points = []() {
    auto* const points = new std::array<double, N>;
    for (std::int64_t k = 1; k <= N; k += 2) {
      // These points are always computed using a correctly-rounded
      // implementation.  They may be evaluated at any time during the
      // execution, and we wouldn't want them to depend on the current
      // configuration of the elementary functions.
      double const πk_over_N = π * k / N;
      (*points)[k] = CanUseHardwareFMA ? Cos<FMAPresence::Present>(πk_over_N)
                                       : Cos<FMAPresence::Absent>(πk_over_N);
    }
    return points;
  }();

  DCHECK_EQ(k % 2, 1);
  return (*points)[k];
}

}  // namespace internal
}  // namespace _чебышёв_lobatto
}  // namespace numerics
}  // namespace principia
