
#pragma once

#include <type_traits>

#include "base/not_constructible.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_traits {

using base::not_constructible;
using quantities::Difference;

// A type trait for testing if a type is a vector.
template<typename T, typename = void>
struct is_vector : std::false_type, not_constructible {};
template<typename T>
struct is_vector<T, std::enable_if_t<std::is_same_v<T, Difference<T, T>>>>
    : std::true_type, not_constructible {};
template<typename T>
struct is_vector<T const> : is_vector<T> {};

template<typename T>
constexpr bool is_vector_v = is_vector<T>::value;

}  // namespace internal_traits

using internal_traits::is_vector;
using internal_traits::is_vector_v;

}  // namespace geometry
}  // namespace principia
