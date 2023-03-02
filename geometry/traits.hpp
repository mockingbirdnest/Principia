#pragma once

#include <type_traits>

#include "base/not_constructible.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _traits {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::quantities::_named_quantities;

// A type trait for testing if a type is a member of a vector space (as opposed
// to an affine space).
template<typename T, typename = void>
struct is_vector : std::false_type, not_constructible {};
template<typename T>
struct is_vector<T, std::enable_if_t<std::is_same_v<T, Difference<T, T>>>>
    : std::true_type, not_constructible {};
template<typename T>
struct is_vector<T const> : is_vector<T> {};

template<typename T>
constexpr bool is_vector_v = is_vector<T>::value;

}  // namespace internal

using internal::is_vector;
using internal::is_vector_v;

}  // namespace _traits
}  // namespace geometry
}  // namespace principia

namespace principia::geometry {
using namespace principia::geometry::_traits;
}  // namespace principia::geometry
