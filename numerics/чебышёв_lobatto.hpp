#pragma once

#include <cstdint>

namespace principia {
namespace numerics {
namespace _чебышёв_lobatto {
namespace internal {

// Returns the Чебышёв-Lobatto point Cos(πk/N).  Should only be called with odd
// k (values for even k should be obtained for N/2).
template<int N>
double ЧебышёвLobattoPoint(std::int64_t k);

}  // namespace internal

using internal::ЧебышёвLobattoPoint;

}  // namespace _чебышёв_lobatto
}  // namespace numerics
}  // namespace principia

#include "numerics/чебышёв_lobatto_body.hpp"
