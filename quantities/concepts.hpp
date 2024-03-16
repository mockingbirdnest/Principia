#pragma once

#include <concepts>
#include <type_traits>

#include "base/traits.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace _concepts {
namespace internal {

using namespace base::_traits;
using namespace quantities::_quantities;

template<typename T>
concept quantity = std::integral<T> || std::floating_point<T> ||
                   is_instance_of_v<Quantity, std::remove_const_t<T>>;

}  // namespace internal

using internal::quantity;

}  // namespace _concepts
}  // namespace quantities
}  // namespace principia
