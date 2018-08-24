
#pragma once

#include <algorithm>
#include <tuple>
#include <utility>

namespace principia {
namespace quantities {
namespace internal_tuples {

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
//TODO(phl):index_sequence vs. integer_sequence?
template<template<typename, typename> class Transform,
         typename LTuple, typename RTuple,
         typename = std::make_index_sequence<
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
struct Apply2Generator;

// This struct has a |Type| member which is an n-element tuple of successive
// derivatives of |Value| with respect to |Argument|; the first element is
// |Value|.
template<typename Value, typename Argument, int n,
         typename = std::make_integer_sequence<int, n>>
struct DerivativesGenerator;

}  // namespace internal_tuples

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
