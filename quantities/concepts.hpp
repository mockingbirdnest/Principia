#pragma once

#include <concepts>
#include <type_traits>

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace _concepts {
namespace internal {

using namespace quantities::_quantities;

template<typename T>
concept quantity = std::integral<T> || std::floating_point<T> || requires(T q) {
  { Quantity{q} } -> std::same_as<std::remove_const_t<T>>;
};

}  // namespace internal

using internal::quantity;

}  // namespace _concepts
}  // namespace quantities
}  // namespace principia
