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
namespace vector_space {

using namespace principia::quantities::_tuples;

template<typename RTuple>
constexpr auto operator+(RTuple const& right);

template<typename RTuple>
constexpr auto operator-(RTuple const& right);

template<typename LTuple, typename RTuple>
constexpr auto operator+(LTuple const& left, RTuple const& right);

template<typename LTuple, typename RTuple>
constexpr auto operator-(LTuple const& left, RTuple const& right);

template<typename Scalar, typename Tuple,
         typename = std::enable_if_t<!is_tuple_v<Scalar>>,
         typename = std::enable_if_t<is_tuple_v<Tuple>>>
constexpr auto operator*(Scalar const& left, Tuple const& right);

// The extra typename lifts an ambiguity on the definition.
template<typename Tuple, typename Scalar,
         typename = std::enable_if_t<is_tuple_v<Tuple>>,
         typename = std::enable_if_t<!is_tuple_v<Scalar>>,
         typename = void>
constexpr auto operator*(Tuple const& left, Scalar const& right);

template<typename Scalar, typename Tuple>
constexpr auto operator/(Tuple const& left, Scalar const& right);

}  // namespace vector_space

namespace polynomial_ring {

using namespace principia::quantities::_tuples;

// The product assumes that the tuple elements are in the monomial basis.
template<typename LTuple, typename RTuple,
         typename = std::enable_if_t<is_tuple_v<LTuple>>,
         typename = std::enable_if_t<is_tuple_v<RTuple>>>
constexpr auto operator*(LTuple const& left, RTuple const& right);

template<int exponent, typename Tuple>
constexpr auto Pow(Tuple const& tuple);

}  // namespace polynomial_ring

namespace pointwise_inner_product {

using namespace principia::quantities::_tuples;

template<typename LTuple, typename RTuple,
         typename = std::enable_if_t<is_tuple_v<LTuple>>,
         typename = std::enable_if_t<is_tuple_v<RTuple>>>
constexpr auto PointwiseInnerProduct(LTuple const& left, RTuple const& right);

}  // namespace pointwise_inner_product

namespace internal {

using namespace principia::quantities::_tuples;

// Tuple wrapper with vector operations always defined on it.
//
// (The vector operations defined in `vector_space` aren't always available).
template<tuple T>
struct VectorTuple {
  friend auto operator<=>(VectorTuple const& left,
                          VectorTuple const& right) = default;

  template<tuple U>
  VectorTuple<T>& operator+=(VectorTuple<U> const& right);
  template<tuple U>
  VectorTuple<T>& operator-=(VectorTuple<U> const& right);

  template<typename Scalar>
  VectorTuple<T>& operator*=(Scalar const& right);
  template<typename Scalar>
  VectorTuple<T>& operator/=(Scalar const& right);

  T tuple;
};

template<tuple T>
constexpr auto operator+(VectorTuple<T> const& right);

template<tuple T>
constexpr auto operator-(VectorTuple<T> const& right);

template<tuple L, tuple R>
constexpr auto operator+(VectorTuple<L> const& left,
                         VectorTuple<R> const& right);

template<tuple L, tuple R>
constexpr auto operator-(VectorTuple<L> const& left,
                         VectorTuple<R> const& right);

template<typename L, tuple R>
constexpr auto operator*(L const& left, VectorTuple<R> const& right);

template<tuple L, typename R>
constexpr auto operator*(VectorTuple<L> const& left, R const& right);

template<tuple L, typename R>
constexpr auto operator/(VectorTuple<L> const& left, R const& right);

}  // namespace internal

using internal::VectorTuple;

}  // namespace _cartesian_product
}  // namespace geometry
}  // namespace principia

#include "geometry/cartesian_product_body.hpp"
