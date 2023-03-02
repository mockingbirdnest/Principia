#pragma once

#include <type_traits>

#include "base/not_constructible.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace _traits {
namespace internal {

using namespace principia::base::_not_constructible;

// A type trait for testing if a type is a quantity.
template<typename T>
struct is_quantity : std::is_arithmetic<T>, not_constructible {};
template<typename D>
struct is_quantity<Quantity<D>> : std::true_type, not_constructible {};
template<typename T>
struct is_quantity<T const> : is_quantity<T> {};

template<typename T>
constexpr bool is_quantity_v = is_quantity<T>::value;

}  // namespace internal

using internal::is_quantity;
using internal::is_quantity_v;

}  // namespace _traits
}  // namespace quantities
}  // namespace principia

namespace principia::quantities {
using namespace principia::quantities::_traits;
}  // namespace principia::quantities
