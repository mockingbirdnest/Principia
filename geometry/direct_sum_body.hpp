#include "geometry/direct_sum.hpp"

#include <algorithm>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "absl/strings/str_join.h"
#include "base/for_all_of.hpp"
#include "numerics/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _direct_sum {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::numerics::_elementary_functions;

template<affine... T>
class Uninitializer;

template<affine T0, affine... Ti>
class Uninitializer<T0, Ti...> {
 public:
  static constexpr std::tuple<T0, Ti...> Make() {
    if constexpr (uninitialized_constructible<T0>) {
      T0 const t0(uninitialized);
      return std::tuple_cat(std::tie(t0), Uninitializer<Ti...>::Make());
    } else {
      T0 t0;
      return std::tuple_cat(std::tie(t0), Uninitializer<Ti...>::Make());
    }
  }
};

template<>
class Uninitializer<> {
 public:
  static constexpr std::tuple<> Make() {
    return {};
  }
};

template<affine... T>
constexpr DirectSum<T...>::DirectSum(uninitialized_t)
    : tuple_(Uninitializer<T...>::Make()) {}

template<affine... T>
constexpr DirectSum<T...>::DirectSum(T const&... args) : tuple_(args...) {}

template<affine... T>
template<typename... Args>
  requires(sizeof...(T) >= 1 && (std::constructible_from<T, Args> && ...))
constexpr DirectSum<T...>::DirectSum(Args&&... args)
    : tuple_(std::forward<Args>(args)...) {}

template<affine... T>
constexpr auto DirectSum<T...>::Norm() const
  requires(hilbert<T> && ...)
{
  return Sqrt(Norm²());
}

template<affine... T>
constexpr auto DirectSum<T...>::Norm²() const
  requires(hilbert<T> && ...)
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
DirectSum<T...>& DirectSum<T...>::operator*=(Scalar const& right)
  requires(module<T, Scalar> && ...)
{
  *this = *this * right;
  return *this;
}

template<affine... T>
template<field Scalar>
DirectSum<T...>& DirectSum<T...>::operator/=(Scalar const& right)
  requires(vector_space<T, Scalar> && ...)
{
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
  requires(!additive_group<T> || ...)
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
  requires(!additive_group<T> || ...)
{
  DirectSum<T...> difference(uninitialized);
  for_all_of(left, right, difference)
      .loop([](auto const& left, auto const& right, auto& difference) {
        difference = left - right;
      });
  return difference;
}

template<typename L, typename... R>
  requires(!is_instance_of_v<DirectSum, L>) && (homogeneous_module<R, L> && ...)
constexpr auto operator*(L const& left, DirectSum<R...> const& right) {
  DirectSum<Product<L, R>...> product(uninitialized);
  for_all_of(right, product).loop([&left](auto const& right, auto& product) {
    product = left * right;
  });
  return product;
}

template<typename... L, typename R>
  requires(!is_instance_of_v<DirectSum, R>) && (homogeneous_module<L, R> && ...)
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
  requires(hilbert<T> && ...)
constexpr auto InnerProduct(DirectSum<T...> const& left,
                            DirectSum<T...> const& right) {
  using T0 = std::tuple_element_t<0, DirectSum<T...>>;
  InnerProductType<T0, T0> product = {};
  for_all_of(left, right)
      .loop([&product](auto const& leftᵢ, auto const& rightᵢ) {
        using geometry::_hilbert::InnerProduct;
        product += InnerProduct(leftᵢ, rightᵢ);
      });

  return product;
}

template<affine... T>
std::string DebugString(DirectSum<T...> const& direct_sum) {
  std::vector<std::string> strings;
  strings.reserve(sizeof...(T));
  for_all_of(direct_sum).loop([&strings](auto const& component) {
    strings.push_back(DebugString(component));
  });
  return absl::StrCat("{", absl::StrJoin(strings, ", "), "}");
}

template<affine... T>
std::ostream& operator<<(std::ostream& out, DirectSum<T...> const& direct_sum) {
  out << "{";
  bool first = true;
  for_all_of(direct_sum).loop([&out, &first](auto const& component) {
    if (!first) {
      out << ", ";
    }
    out << component;
    first = false;
  });
  out << "}";
  return out;
}

}  // namespace internal
}  // namespace _direct_sum
}  // namespace geometry
}  // namespace principia
