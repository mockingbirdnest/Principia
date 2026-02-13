#pragma once

#include <algorithm>
#include <tuple>
#include <type_traits>
#include <utility>

#include "base/algebra.hpp"
#include "base/tags.hpp"
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

using namespace principia::base::_algebra;
using namespace principia::base::_tags;
using namespace principia::quantities::_tuples;

// The direct sum of a pack of affine types.
//
// Has as much structure as the least powerful of its components (a `DirectSum`
// of affine R-modules is an affine R-module, a `DirectSum` of additive groups
// is an additive group, etc.).
template<affine... T>
struct DirectSum {
  constexpr DirectSum() = default;
  constexpr explicit DirectSum(uninitialized_t);

  // Constructor from elements.
  constexpr DirectSum(T&&... args);  // NOLINT(runtime/explicit)

  // Constructor from tuple.
  constexpr explicit DirectSum(std::tuple<T...>&& tuple);

  template<std::size_t i, typename Self>
  constexpr auto&& get(this Self&& self);

  constexpr auto Norm() const
    requires hilbert<DirectSum<T...>, DirectSum<T...>>;
  constexpr auto Norm²() const
    requires hilbert<DirectSum<T...>, DirectSum<T...>>;

  bool operator==(DirectSum const&) const = default;
  bool operator!=(DirectSum const&) const = default;

  template<affine... U>
  DirectSum<T...>& operator+=(DirectSum<U...> const& right);
  template<affine... U>
  DirectSum<T...>& operator-=(DirectSum<U...> const& right);

  template<ring Scalar>
  DirectSum<T...>& operator*=(Scalar const& right);
  template<field Scalar>
  DirectSum<T...>& operator/=(Scalar const& right);

  std::tuple<T...> tuple;
};

template<additive_group... T>
constexpr auto operator+(DirectSum<T...> const& right);

template<additive_group... T>
constexpr auto operator-(DirectSum<T...> const& right);

template<affine... L, affine... R>
constexpr auto operator+(DirectSum<L...> const& left,
                         DirectSum<R...> const& right);

template<affine... L, affine... R>
constexpr auto operator-(DirectSum<L...> const& left,
                         DirectSum<R...> const& right);

template<homogeneous_ring L, homogeneous_module<L>... R>
constexpr auto operator*(L const& left, DirectSum<R...> const& right);

template<homogeneous_ring R, homogeneous_module<R>... L>
constexpr auto operator*(DirectSum<L...> const& left, R const& right);

template<homogeneous_field R, homogeneous_vector_space<R>... L>
constexpr auto operator/(DirectSum<L...> const& left, R const& right);

template<affine... T>
constexpr auto InnerProduct(DirectSum<T...> const& left,
                            DirectSum<T...> const& right);

// Helper for getting a DirectSum corresponding to a tuple when you don't
// have access to the pack.
template<typename Tuple>
struct direct_sum {};

template<affine... T>
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
struct std::
    tuple_element<I, principia::geometry::_cartesian_product::DirectSum<T...>> {
  using type = std::tuple_element_t<I, std::tuple<T...>>;
};

}  // namespace std

#include "geometry/cartesian_product_body.hpp"
