#pragma once

#include <algorithm>
#include <tuple>
#include <type_traits>
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

template<typename Scalar,
         typename Tuple,
         typename = std::enable_if_t<!is_tuple_v<Scalar>>,
         typename = std::enable_if_t<is_tuple_v<Tuple>>>
constexpr auto operator*(Scalar const& left, Tuple const& right);

// The extra typename lifts an ambiguity on the definition.
template<typename Tuple,
         typename Scalar,
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
template<typename LTuple,
         typename RTuple,
         typename = std::enable_if_t<is_tuple_v<LTuple>>,
         typename = std::enable_if_t<is_tuple_v<RTuple>>>
constexpr auto operator*(LTuple const& left, RTuple const& right);

template<int exponent, typename Tuple>
constexpr auto Pow(Tuple const& tuple);

}  // namespace polynomial_ring

namespace pointwise_inner_product {

using namespace principia::quantities::_tuples;

template<typename LTuple,
         typename RTuple,
         typename = std::enable_if_t<is_tuple_v<LTuple>>,
         typename = std::enable_if_t<is_tuple_v<RTuple>>>
constexpr auto PointwiseInnerProduct(LTuple const& left, RTuple const& right);

}  // namespace pointwise_inner_product

namespace internal {

using namespace principia::quantities::_tuples;

template<typename... T>
struct DirectSum {
  constexpr DirectSum() = default;

  // Constructor from elements.
  constexpr DirectSum(T&&... args);

  // Constructor from tuple.
  constexpr explicit DirectSum(std::tuple<T...>&& tuple);

  template<std::size_t I>
  constexpr auto const& get() const;
  template<std::size_t I>
  constexpr auto& get();

  friend auto operator<=>(DirectSum const& left,
                          DirectSum const& right) = default;

  template<typename... U>
  DirectSum<T...>& operator+=(DirectSum<U...> const& right);
  template<typename... U>
  DirectSum<T...>& operator-=(DirectSum<U...> const& right);

  template<typename Scalar>
  DirectSum<T...>& operator*=(Scalar const& right);
  template<typename Scalar>
  DirectSum<T...>& operator/=(Scalar const& right);

  std::tuple<T...> tuple;
};

template<typename... T>
constexpr auto operator+(DirectSum<T...> const& right);

template<typename... T>
constexpr auto operator-(DirectSum<T...> const& right);

template<typename... L, typename... R>
constexpr auto operator+(DirectSum<L...> const& left,
                         DirectSum<R...> const& right);

template<typename... L, typename... R>
constexpr auto operator-(DirectSum<L...> const& left,
                         DirectSum<R...> const& right);

template<typename L, typename... R>
constexpr auto operator*(L const& left, DirectSum<R...> const& right);

template<typename... L, typename R>
constexpr auto operator*(DirectSum<L...> const& left, R const& right);

template<typename... L, typename R>
constexpr auto operator/(DirectSum<L...> const& left, R const& right);

// Helper for getting a DirectSum corresponding to a tuple when you don't
// have access to the pack.
template<tuple Tuple>
struct direct_sum {};

template<typename... T>
struct direct_sum<std::tuple<T...>> {
  typedef DirectSum<T...> type;
};

template<typename T>
using direct_sum_t = direct_sum<T>::type;

}  // namespace internal

using internal::direct_sum;
using internal::direct_sum_t;
using internal::DirectSum;

}  // namespace _cartesian_product
}  // namespace geometry
}  // namespace principia

// Specializations of tuple traits to enable structured bindings.
namespace std {

template<typename... T>
struct std::tuple_size<principia::geometry::_cartesian_product::DirectSum<T...>>
    : public std::integral_constant<std::size_t,
                                    std::tuple_size_v<std::tuple<T...>>> {};

template<std::size_t I, typename... T>
struct std::tuple_element<I, principia::geometry::_cartesian_product::DirectSum<T...>> {
  using type = std::tuple_element_t<I, std::tuple<T...>>;
};

}  // namespace std

#include "geometry/cartesian_product_body.hpp"
