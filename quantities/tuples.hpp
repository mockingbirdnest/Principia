
#pragma once

#include <algorithm>
#include <tuple>
#include <utility>

namespace principia {
namespace quantities {
namespace internal_tuples {

// A trait for finding if something is a tuple.
template<typename T>
struct is_tuple : std::false_type, not_constructible {};
template<typename... D>
struct is_tuple<std::tuple<D...>> : std::true_type, not_constructible {};

template<typename T>
constexpr bool is_tuple_v = is_tuple<T>::value;

// This struct has a |Type| member which is a tuple obtained by applying
// |Transform| to each element type in |Tuple| (which must be a tuple or an
// array or a pair).
template<template<typename> class Transform,
         typename Tuple,
         typename = std::make_integer_sequence<int, std::tuple_size_v<Tuple>>>
struct ApplyGenerator;

// Same as above, but |Transform| is applied to corresponding pairs of element
// types from |LTuple| and |RTuple|.  If the tuples have different sizes, |void|
// is passed to |Transform| for the missing element types.
template<template<typename, typename> class Transform,
         typename LTuple, typename RTuple,
         typename = std::make_integer_sequence<
             int,
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
struct Apply2Generator;

// This struct has a |Type| member which is an n-element tuple of successive
// derivatives of |Value| with respect to |Argument|; the first element is
// |Value|.
template<typename Value, typename Argument, int n,
         typename = std::make_integer_sequence<int, n>>
struct DerivativesGenerator;

}  // namespace internal_tuples

using internal_tuples::is_tuple;
using internal_tuples::is_tuple_v;

template<template<typename> class Transform, typename Tuple>
using Apply = typename internal_tuples::ApplyGenerator<Transform, Tuple>::Type;

template<template<typename, typename> class Transform,
         typename LTuple, typename RTuple>
using Apply2 =
    typename internal_tuples::Apply2Generator<Transform, LTuple, RTuple>::Type;

template<typename Value, typename Argument, int n>
using Derivatives =
    typename internal_tuples::DerivativesGenerator<Value, Argument, n>::Type;

}  // namespace quantities
}  // namespace principia

#include "quantities/tuples_body.hpp"
