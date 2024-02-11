#pragma once

#include "numerics/чебышёв_lobatto.hpp"

#include <array>
#include <cstdint>

#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _чебышёв_lobatto {
namespace internal {

using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

// Returns the Чебышёв-Lobatto point Cos(πk/N).  Should only be called with odd
// k (values for even k should be obtained by calling the N/2 instantiation).
template<int N>
double ЧебышёвLobattoPoint(std::int64_t const k) {
  static const std::array<double, N>* const points = []() {
    auto* const points = new std::array<double, N>;
    for (std::int64_t k = 1; k <= N; k += 2) {
      (*points)[k] = Cos(π * k * Radian / N);
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
