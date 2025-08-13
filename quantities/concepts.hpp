#pragma once

#include <type_traits>

#include "base/traits.hpp"
#include "quantities/cantor.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace _concepts {
namespace internal {

using namespace principia::base::_traits;
using namespace principia::quantities::_cantor;
using namespace principia::quantities::_quantities;

template<typename T>
concept quantity = instance<T, Quantity> || continuum<T>;

// std::integral || std::floating_point rather than
// std::convertible_to<double, T> because
// the latter introduces ambiguities on Sign * Vector.
template<typename T>
concept convertible_to_quantity =
    quantity<std::remove_cvref_t<T>> || countable<T> || continuum<T>;

}  // namespace internal

using internal::convertible_to_quantity;
using internal::quantity;

}  // namespace _concepts
}  // namespace quantities
}  // namespace principia
