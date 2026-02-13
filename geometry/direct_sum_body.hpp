#include "geometry/direct_sum.hpp"

#include <algorithm>
#include <tuple>
#include <type_traits>

#include "base/for_all_of.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _direct_sum {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::geometry::_hilbert;
using namespace principia::numerics::_elementary_functions;

template<affine... T>
constexpr DirectSum<T...>::DirectSum(uninitialized_t) {}

template<affine... T>
constexpr DirectSum<T...>::DirectSum(T&&... args) : tuple(args...) {}

template<affine... T>
constexpr DirectSum<T...>::DirectSum(std::tuple<T...>&& tuple) : tuple(tuple) {}

template<std::size_t i, affine... T>
constexpr auto const& get(DirectSum<T...> const& self) {
  return std::get<i>(self.tuple);
}
template<std::size_t i, affine... T>
constexpr auto& get(DirectSum<T...>& self) {
  return std::get<i>(self.tuple);
}

template<affine... T>
constexpr auto DirectSum<T...>::Norm() const
  requires hilbert<DirectSum<T...>, DirectSum<T...>>
{
  return Sqrt(Norm²());
}

template<affine... T>
constexpr auto DirectSum<T...>::Norm²() const
  requires hilbert<DirectSum<T...>, DirectSum<T...>>
{
  return InnerProduct(*this, *this);
}

template<affine... T>
DirectSum<T...>& DirectSum<T...>::operator+=(
    DirectSum<Difference<T>...> const& right) {
  *this = *this + right;
  return *this;
}

template<affine... T>
DirectSum<T...>& DirectSum<T...>::operator-=(
    DirectSum<Difference<T>...> const& right) {
  *this = *this - right;
  return *this;
}

template<affine... T>
template<ring Scalar>
DirectSum<T...>& DirectSum<T...>::operator*=(Scalar const& right) {
  *this = *this * right;
  return *this;
}

template<affine... T>
template<field Scalar>
DirectSum<T...>& DirectSum<T...>::operator/=(Scalar const& right) {
  *this = *this / right;
  return *this;
}

template<additive_group... T>
constexpr DirectSum<T...> operator+(DirectSum<T...> const& right) {
  return right;
}

template<additive_group... T>
constexpr DirectSum<T...> operator-(DirectSum<T...> const& right) {
  DirectSum<T...> const zero{};
  return zero - right;
}

template<affine... T>
constexpr DirectSum<T...> operator+(DirectSum<T...> const& left,
                                    DirectSum<Difference<T>...> const& right) {
  DirectSum<T...> sum(uninitialized);
  for_all_of(left, right, sum)
      .loop([](auto const& left, auto const& right, auto& sum) {
        sum = left + right;
      });
  return sum;
}

template<affine... T>
constexpr DirectSum<T...> operator+(DirectSum<Difference<T>...> const& left,
                                    DirectSum<T...> const& right)
  requires(!std::is_same<DirectSum<T...>, DirectSum<Difference<T>...>>::value)
{
  DirectSum<T...> sum(uninitialized);
  for_all_of(left, right, sum)
      .loop([](auto const& left, auto const& right, auto& sum) {
        sum = left + right;
      });
  return sum;
}

template<affine... T>
constexpr DirectSum<Difference<T>...> operator-(DirectSum<T...> const& left,
                                                DirectSum<T...> const& right) {
  DirectSum<Difference<T>...> difference(uninitialized);
  for_all_of(left, right, difference)
      .loop([](auto const& left, auto const& right, auto& difference) {
        difference = left - right;
      });
  return difference;
}

template<affine... T>
constexpr DirectSum<T...> operator-(DirectSum<T...> const& left,
                                    DirectSum<Difference<T>...> const& right)
  requires(!std::is_same<DirectSum<T...>, DirectSum<Difference<T>...>>::value)
{
  DirectSum<T...> difference(uninitialized);
  for_all_of(left, right, difference)
      .loop([](auto const& left, auto const& right, auto& difference) {
        difference = left - right;
      });
  return difference;
}

template<homogeneous_ring L, homogeneous_module<L>... R>
constexpr auto operator*(L const& left, DirectSum<R...> const& right) {
  DirectSum<Product<L, R>...> product(uninitialized);
  for_all_of(right, product).loop([&left](auto const& right, auto& product) {
    product = left * right;
  });
  return product;
}

template<homogeneous_ring R, homogeneous_module<R>... L>
constexpr auto operator*(DirectSum<L...> const& left, R const& right) {
  DirectSum<Product<L, R>...> product(uninitialized);
  for_all_of(left, product).loop([&right](auto const& left, auto& product) {
    product = left * right;
  });
  return product;
}

template<homogeneous_field R, homogeneous_vector_space<R>... L>
constexpr auto operator/(DirectSum<L...> const& left, R const& right) {
  DirectSum<Quotient<L, R>...> quotient(uninitialized);
  for_all_of(left, quotient).loop([&right](auto const& left, auto& quotient) {
    quotient = left / right;
  });
  return quotient;
}

template<affine... T>
constexpr auto InnerProduct(DirectSum<T...> const& left,
                            DirectSum<T...> const& right) {
  using T0 = std::tuple_element_t<0, DirectSum<T...>>;
  decltype(std::declval<T0>() * std::declval<T0>()) product = {};
  for_all_of(left, right)
      .loop([&product](auto const& leftᵢ, auto const& rightᵢ) {
        product += leftᵢ * rightᵢ;
      });

  return product;
}

}  // namespace internal
}  // namespace _direct_sum
}  // namespace geometry
}  // namespace principia
