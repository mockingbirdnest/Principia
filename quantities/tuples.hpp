#pragma once

#include <algorithm>
#include <tuple>
#include <utility>

#include "base/not_constructible.hpp"

namespace principia {
namespace quantities {
namespace _tuples {
namespace internal {

using namespace principia::base::_not_constructible;

// A trait for finding if something is a tuple.
// TODO(phl): We might want to use this for pair and array too.
template<typename T>
struct is_tuple : std::false_type, not_constructible {};
template<typename... D>
struct is_tuple<std::tuple<D...>> : std::true_type, not_constructible {};

template<typename T>
constexpr bool is_tuple_v = is_tuple<T>::value;

// This struct has a |Type| member which is a tuple obtained by applying
// |Transform| to corresponding elements in |Tuples| (which may be tuples,
// arrays, or pairs).  If the |Tuples| have different sizes, |void| is passed to
// |Transform| for the missing element types.
template<template<typename...> typename Transform, typename... Tuples>
class ApplyGenerator;

// This struct has a |Type| member which is an n-element tuple of successive
// derivatives of |Value| with respect to |Argument|; the first element is
// |Value|.
template<typename Value, typename Argument, int n,
         typename = std::make_index_sequence<n>>
struct DerivativesGenerator;

}  // namespace internal

using internal::is_tuple;
using internal::is_tuple_v;

template<template<typename...> typename Transform, typename... Tuples>
using Apply = typename internal::ApplyGenerator<Transform, Tuples...>::Type;

template<typename Value, typename Argument, int n>
using Derivatives =
    typename internal::DerivativesGenerator<Value, Argument, n>::Type;

}  // namespace _tuples
}  // namespace quantities
}  // namespace principia

#include "quantities/tuples_body.hpp"
