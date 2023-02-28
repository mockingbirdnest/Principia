#pragma once

#include <algorithm>
#include <tuple>
#include <utility>

#include "quantities/tuples.hpp"

namespace principia {
namespace geometry {

// These operators live in a separate (non-internal) namespace to avoid
// polluting the entire universe in cases where they are not useful.
namespace _cartesian_product {
namespace internal {

template<typename RTuple>
constexpr auto operator+(RTuple const& right);

template<typename RTuple>
constexpr auto operator-(RTuple const& right);

template<typename LTuple, typename RTuple>
constexpr auto operator+(LTuple const& left, RTuple const& right);

template<typename LTuple, typename RTuple>
constexpr auto operator-(LTuple const& left, RTuple const& right);

template<typename Scalar, typename Tuple,
         typename = std::enable_if_t<!quantities::is_tuple_v<Scalar>>,
         typename = std::enable_if_t<quantities::is_tuple_v<Tuple>>>
constexpr auto operator*(Scalar const& left, Tuple const& right);

// The extra typename lifts an ambiguity on the definition.
template<typename Tuple, typename Scalar,
         typename = std::enable_if_t<quantities::is_tuple_v<Tuple>>,
         typename = std::enable_if_t<!quantities::is_tuple_v<Scalar>>,
         typename = void>
constexpr auto operator*(Tuple const& left, Scalar const& right);

template<typename Scalar, typename Tuple>
constexpr auto operator/(Tuple const& left, Scalar const& right);

}  // namespace internal

using internal::operator-;
using internal::operator*;
using internal::operator/;
using internal::operator+;

}  // namespace _cartesian_product

namespace _polynomial_ring {
namespace internal {

// The product assumes that the tuple elements are in the monomial basis.
template<typename LTuple, typename RTuple,
         typename = std::enable_if_t<quantities::is_tuple_v<LTuple>>,
         typename = std::enable_if_t<quantities::is_tuple_v<RTuple>>>
constexpr auto operator*(LTuple const& left, RTuple const& right);

template<int exponent, typename Tuple>
constexpr auto Pow(Tuple const& tuple);

}  // namespace internal

using internal::operator*;
using internal::Pow;

}  // namespace _polynomial_ring

namespace _pointwise_inner_product {
namespace internal {

template<typename LTuple, typename RTuple,
         typename = std::enable_if_t<quantities::is_tuple_v<LTuple>>,
         typename = std::enable_if_t<quantities::is_tuple_v<RTuple>>>
constexpr auto PointwiseInnerProduct(LTuple const& left, RTuple const& right);

}  // namespace internal

using internal::PointwiseInnerProduct;

}  // namespace _pointwise_inner_product

}  // namespace geometry
}  // namespace principia

#include "geometry/cartesian_product_body.hpp"
