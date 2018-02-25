
#pragma once

#include <pmmintrin.h>

#include <type_traits>

#include "quantities/traits.hpp"

namespace principia {
namespace quantities {

namespace internal_quantities {
template<typename D>
class Quantity;
}  // namespace internal_quantities

namespace internal_wide {

using internal_quantities::Quantity;

// Fills both halves of the result.
template<typename D>
__m128d ToM128D(Quantity<D> x);
template<typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
__m128d ToM128D(T x);

}  // namespace internal_wide

using internal_wide::ToM128D;

}  // namespace quantities
}  // namespace principia

// Because of circular dependencies, this file doesn't include wide_body.hpp.
// This will be done by quantities.hpp.
