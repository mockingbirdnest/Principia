#pragma once

#include <algorithm>
#include <tuple>
#include <type_traits>
#include <utility>

#include "base/algebra.hpp"
#include "base/not_constructible.hpp"
#include "base/tags.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace geometry {
namespace _direct_sum {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_not_constructible;
using namespace principia::base::_tags;
using namespace principia::quantities::_tuples;

// The direct sum of a pack of affine types.
//
// Has as much structure as the least powerful of its components (a `DirectSum`
// of affine R-modules is an affine R-module, a `DirectSum` of additive groups
// is an additive group, etc.).
template<affine... T>
class DirectSum {
 public:
  constexpr DirectSum() = default;
  constexpr explicit DirectSum(uninitialized_t);

  // Constructor from elements.
  constexpr DirectSum(T&&... args);  // NOLINT(runtime/explicit)

  // Constructor from tuple.
  constexpr explicit DirectSum(std::tuple<T...>&& tuple);

  // Getter for the inner tuple field.
  template<typename Self>
  auto&& tuple(this Self&& self);

  constexpr auto Norm() const
    requires hilbert<DirectSum<T...>, DirectSum<T...>>;
  constexpr auto Norm²() const
    requires hilbert<DirectSum<T...>, DirectSum<T...>>;

  bool operator==(DirectSum const&) const = default;

  DirectSum& operator+=(DirectSum<Difference<T>...> const& right);
  DirectSum& operator-=(DirectSum<Difference<T>...> const& right);

  template<ring Scalar>
  DirectSum& operator*=(Scalar const& right)
    requires(module<T, Scalar> && ...);
  template<field Scalar>
  DirectSum& operator/=(Scalar const& right)
    requires(vector_space<T, Scalar> && ...);

 private:
  std::tuple<T...> tuple_;
};

template<std::size_t i, affine... T>
constexpr auto const& get(DirectSum<T...> const& self);
template<std::size_t i, affine... T>
constexpr auto& get(DirectSum<T...>& self);

template<additive_group... T>
constexpr DirectSum<T...> operator+(DirectSum<T...> const& right);

template<additive_group... T>
constexpr DirectSum<T...> operator-(DirectSum<T...> const& right);

template<affine... T>
constexpr DirectSum<T...> operator+(DirectSum<T...> const& left,
                                    DirectSum<Difference<T>...> const& right);
template<affine... T>
constexpr DirectSum<T...> operator+(DirectSum<Difference<T>...> const& left,
                                    DirectSum<T...> const& right)
  requires(!additive_group<T> || ...);

template<affine... T>
constexpr DirectSum<Difference<T>...> operator-(DirectSum<T...> const& left,
                                                DirectSum<T...> const& right);
template<affine... T>
constexpr DirectSum<T...> operator-(DirectSum<T...> const& left,
                                    DirectSum<Difference<T>...> const& right)
  requires(!additive_group<T> || ...);

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
struct direct_sum;

template<affine... T>
struct direct_sum<std::tuple<T...>> : not_constructible {
  typedef DirectSum<T...> type;
};

template<typename T>
using direct_sum_t = direct_sum<T>::type;

}  // namespace internal

using internal::direct_sum;
using internal::direct_sum_t;
using internal::DirectSum;

}  // namespace _direct_sum
}  // namespace geometry
}  // namespace principia

// Specializations of tuple traits to enable structured bindings.
namespace std {

template<typename... T>
struct std::tuple_size<principia::geometry::_direct_sum::DirectSum<T...>>
    : public std::integral_constant<std::size_t,
                                    std::tuple_size_v<std::tuple<T...>>> {};

template<std::size_t I, typename... T>
struct std::tuple_element<I,
                          principia::geometry::_direct_sum::DirectSum<T...>> {
  using type = std::tuple_element_t<I, std::tuple<T...>>;
};

}  // namespace std

#include "geometry/direct_sum_body.hpp"
