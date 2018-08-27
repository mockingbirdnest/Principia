
#pragma once

#include <cstdint>

namespace principia {
namespace numerics {

// n choose k.
constexpr std::int64_t Binomial(std::int64_t n, std::int64_t k);

constexpr std::int64_t FallingFactorial(std::int64_t n, std::int64_t k);

}  // namespace numerics
}  // namespace principia

#include "numerics/combinatorics_body.hpp"
