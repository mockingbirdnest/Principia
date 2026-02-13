#include "geometry/cartesian_product.hpp"

#include <algorithm>
#include <tuple>
#include <type_traits>

#include "base/for_all_of.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/elementary_functions.hpp"

// The use of FORCE_INLINE in this file is because we want the construction of
// complex polynomials to be reasonably efficient.

namespace principia {
namespace geometry {
namespace _cartesian_product {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::geometry::_hilbert;
using namespace principia::numerics::_elementary_functions;

template<typename LTuple, typename RTuple,
         typename = std::make_index_sequence<
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
class CartesianProductAddition;

template<typename LTuple, typename RTuple, std::size_t... indices>
class CartesianProductAddition<LTuple, RTuple,
                               std::index_sequence<indices...>> {
  // The type of the result of addition, with suitable specializations for the
  // void case of Apply.
  template<typename L, typename R>
  struct TypesGenerator {
    using Sum = base::_algebra::Sum<L, R>;
  };
  template<typename L>
  struct TypesGenerator<L, void> {
    using Sum = L;
  };
  template<typename R>
  struct TypesGenerator<void, R> {
    using Sum = R;
  };

  // Alias for use as the transform in Apply2.
  template<typename L, typename R>
  using Sum = typename TypesGenerator<L, R>::Sum;

 public:
  FORCE_INLINE(static constexpr) Apply<Sum, LTuple, RTuple> Add(
      LTuple const& left,
      RTuple const& right);

 private:
  // Utility for adding elements at the given index.
  template<std::size_t index>
  static constexpr std::tuple_element_t<index, Apply<Sum, LTuple, RTuple>>
  AddElement(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple,
         typename = std::make_index_sequence<
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
class CartesianProductSubtraction;

template<typename LTuple, typename RTuple, std::size_t... indices>
class CartesianProductSubtraction<LTuple, RTuple,
                                  std::index_sequence<indices...>> {
  // The type of the result of subtraction, with suitable specializations for
  // the void case of Apply.
  template<typename L, typename R>
  struct TypesGenerator {
    using Difference = base::_algebra::Difference<L, R>;
  };
  template<typename L>
  struct TypesGenerator<L, void> {
    using Difference = L;
  };
  template<typename R>
  struct TypesGenerator<void, R> {
    using Difference = R;
  };

  // Alias for use as the transform in Apply2.
  template<typename L, typename R>
  using Difference = typename TypesGenerator<L, R>::Difference;

 public:
  FORCE_INLINE(static constexpr) Apply<Difference, LTuple, RTuple> Subtract(
      LTuple const& left,
      RTuple const& right);

 private:
  // Utility for subtracting elements at the given index.
  template<std::size_t index>
  static constexpr std::tuple_element_t<index,
                                        Apply<Difference, LTuple, RTuple>>
  SubtractElement(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, std::size_t... indices>
constexpr auto CartesianProductAddition<
    LTuple, RTuple,
    std::index_sequence<indices...>>::Add(LTuple const& left,
                                          RTuple const& right)
    -> Apply<Sum, LTuple, RTuple> {
  return {AddElement<indices>(left, right)...};
}

template<typename LTuple, typename RTuple, std::size_t... indices>
template<std::size_t index>
constexpr auto CartesianProductAddition<
    LTuple, RTuple,
    std::index_sequence<indices...>>::AddElement(LTuple const& left,
                                                 RTuple const& right)
    -> std::tuple_element_t<index, Apply<Sum, LTuple, RTuple>> {
  if constexpr (index < std::min(std::tuple_size_v<LTuple>,
                                 std::tuple_size_v<RTuple>)) {
    return std::get<index>(left) + std::get<index>(right);
  } else if constexpr (index < std::tuple_size_v<LTuple>) {
    return std::get<index>(left);
  } else {
    return std::get<index>(right);
  }
}

template<typename LTuple, typename RTuple, std::size_t... indices>
constexpr auto CartesianProductSubtraction<
    LTuple, RTuple,
    std::index_sequence<indices...>>::Subtract(LTuple const& left,
                                               RTuple const& right)
    -> Apply<Difference, LTuple, RTuple> {
  return {SubtractElement<indices>(left, right)...};
}

template<typename LTuple, typename RTuple, std::size_t... indices>
template<std::size_t index>
constexpr auto CartesianProductSubtraction<
    LTuple, RTuple,
    std::index_sequence<indices...>>::SubtractElement(LTuple const& left,
                                                      RTuple const& right)
    -> std::tuple_element_t<index, Apply<Difference, LTuple, RTuple>> {
  if constexpr (index < std::min(std::tuple_size_v<LTuple>,
                                 std::tuple_size_v<RTuple>)) {
    return std::get<index>(left) - std::get<index>(right);
  } else if constexpr (index < std::tuple_size_v<LTuple>) {
    return std::get<index>(left);
  } else {
    return -std::get<index>(right);
  }
}

template<typename Scalar, typename Tuple,
         typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
class CartesianProductVectorSpace;

template<typename Scalar, typename Tuple, std::size_t... indices>
class CartesianProductVectorSpace<Scalar, Tuple,
                                  std::index_sequence<indices...>> {
 public:
  template<typename L, typename R>
  using Product = base::_algebra::Product<L, R>;
  template<typename T>
  using ScalarLeftProduct = Product<Scalar, T>;
  template<typename T>
  using ScalarRightProduct = Product<T, Scalar>;
  template<typename T>
  using Quotient = base::_algebra::Quotient<T, Scalar>;

  FORCE_INLINE(static constexpr) Apply<ScalarLeftProduct, Tuple> Multiply(
      Scalar const& left,
      Tuple const& right);
  FORCE_INLINE(static constexpr) Apply<ScalarRightProduct, Tuple> Multiply(
      Tuple const& left,
      Scalar const& right);

  // The templatization SFINAEs out this operator unless the division is
  // defined.
  template<typename T = Tuple>
  FORCE_INLINE(static constexpr)
  Apply<Quotient, T> Divide(Tuple const& left, Scalar const& right);
};

template<typename Scalar, typename Tuple, std::size_t... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::index_sequence<indices...>>::Multiply(Scalar const& left,
                                               Tuple const& right)
    -> Apply<ScalarLeftProduct, Tuple> {
  return {left * std::get<indices>(right)...};
}

template<typename Scalar, typename Tuple, std::size_t... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::index_sequence<indices...>>::Multiply(Tuple const& left,
                                               Scalar const& right)
    -> Apply<ScalarRightProduct, Tuple> {
  return {std::get<indices>(left) * right...};
}

template<typename Scalar, typename Tuple, std::size_t... indices>
template<typename T>
constexpr auto
CartesianProductVectorSpace<Scalar, Tuple, std::index_sequence<indices...>>::
    Divide(Tuple const& left, Scalar const& right) -> Apply<Quotient, T> {
  return {std::get<indices>(left) / right...};
}

// A helper for prepending an element to a tuple.  Somewhat similar to
// std::tuple_cat but conveniently exports the type of the result.
template<typename Element, typename Tuple,
         typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
struct ConsGenerator;
template<typename Element, typename Tuple, std::size_t... indices>
struct ConsGenerator<Element, Tuple, std::index_sequence<indices...>> {
  using Type = std::tuple<Element, std::tuple_element_t<indices, Tuple>...>;
  static constexpr Type Cons(Element const& element, Tuple const& tuple);
};

template<typename Element, typename Tuple, std::size_t... indices>
constexpr auto ConsGenerator<Element, Tuple, std::index_sequence<indices...>>::
Cons(Element const& element, Tuple const& tuple) -> Type {
  return std::make_tuple(element, std::get<indices>(tuple)...);
}

// A helper for extracting the tail of a tuple, i.e., everything except element
// 0.
template<typename Tuple,
         typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
struct TailGenerator;
template<typename Tuple, std::size_t head_index, std::size_t... tail_indices>
struct TailGenerator<Tuple, std::index_sequence<head_index, tail_indices...>> {
  using Type = std::tuple<std::tuple_element_t<tail_indices, Tuple>...>;
  static constexpr Type Tail(Tuple const& tuple);
};

template<typename Tuple, std::size_t head_index, std::size_t... tail_indices>
constexpr auto
TailGenerator<Tuple, std::index_sequence<head_index, tail_indices...>>::
Tail(Tuple const& tuple) -> Type {
  return std::make_tuple(std::get<tail_indices>(tuple)...);
}

// Technically, this is only a ring if `CartesianProductMultiplicativeSpace` is
// the bona fide `CartesianProductVectorSpace`.  If it is not, this implements a
// multiplication-like operation which has the expected properties with respect
// to the cartesian product additive group.
// NOTE(phl): The need for the third, defaulted, parameter to the template
// template parameter is because Clang, as of 10.0.1, doesn't implement
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2016/p0522r0.html, see
// https://godbolt.org/z/9T8M3P.  Visual Studio 16.6.3 does.
template<typename LTuple, typename RTuple,
         template<typename, typename Tuple,
                  typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
         class CartesianProductMultiplicativeSpace,
         int lsize_ = std::tuple_size_v<LTuple>,
         int rsize_ = std::tuple_size_v<RTuple>>
class PolynomialRing {
 public:
  // Right is split into head (index 0) and tail (the rest).  The tail is a
  // polynomial with valuation 1.
  using RHead = std::tuple_element_t<0, RTuple>;
  using RTail = typename TailGenerator<RTuple>::Type;

  // To implement the polynomial multiplication left * right_tail, we need to
  // insert a zero for the lowest degree (because of the valuation 1).  This is
  // the type of that zero.
  using LHead = std::tuple_element_t<0, LTuple>;
  using Zero = typename CartesianProductMultiplicativeSpace<LHead, RTuple>::
      template Product<LHead, RHead>;

  // The hard part: generating the type of the result.
  using LTupleRHeadProduct =
      decltype(CartesianProductMultiplicativeSpace<RHead, LTuple>::Multiply(
          std::declval<LTuple>(),
          std::declval<RHead>()));
  using LTupleRTailProduct =
      decltype(PolynomialRing<LTuple, RTail,
                              CartesianProductMultiplicativeSpace>::Multiply(
          std::declval<LTuple>(),
          std::declval<RTail>()));
  using ZeroLTupleRTailProduct =
      typename ConsGenerator<Zero, LTupleRTailProduct>::Type;

  using Result =
      decltype(CartesianProductAddition<LTupleRHeadProduct,
                                        ZeroLTupleRTailProduct>::
                   Add(std::declval<LTupleRHeadProduct>(),
                       std::declval<ZeroLTupleRTailProduct>()));

  FORCE_INLINE(static constexpr)
  Result Multiply(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple,
         template<typename, typename Tuple,
                  typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
         class CartesianProductMultiplicativeSpace,
         int lsize_>
class PolynomialRing<LTuple, RTuple,
                     CartesianProductMultiplicativeSpace,
                     lsize_, 1> {
  using RHead = std::tuple_element_t<0, RTuple>;
  using Result = decltype(
      CartesianProductMultiplicativeSpace<RHead, LTuple>::Multiply(
          std::declval<LTuple>(),
          std::declval<RHead>()));

 public:
  FORCE_INLINE(static constexpr)
  Result Multiply(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple,
         template<typename, typename Tuple,
                  typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
         class CartesianProductMultiplicativeSpace,
         int lsize_, int rsize_>
constexpr auto PolynomialRing<LTuple, RTuple,
                              CartesianProductMultiplicativeSpace,
                              lsize_, rsize_>::Multiply(
    LTuple const& left,
    RTuple const& right) -> Result {
  using vector_space::operator+;

  auto const right_head = std::get<0>(right);
  auto const right_tail = TailGenerator<RTuple>::Tail(right);

  using RightHeadType = std::remove_const_t<decltype(right_head)>;
  using RightTailType = std::remove_const_t<decltype(right_tail)>;

  auto const left_times_right_head =
      CartesianProductMultiplicativeSpace<RightHeadType, LTuple>::
          Multiply(left, right_head);
  auto const left_times_right_tail =  PolynomialRing<
      LTuple, RightTailType,
      CartesianProductMultiplicativeSpace>::Multiply(left, right_tail);

  return left_times_right_head +
         ConsGenerator<Zero, LTupleRTailProduct>::Cons(Zero{},
                                                       left_times_right_tail);
}

template<typename LTuple, typename RTuple,
         template<typename, typename Tuple,
                  typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
         class CartesianProductMultiplicativeSpace,
         int lsize_>
constexpr auto PolynomialRing<LTuple, RTuple,
                              CartesianProductMultiplicativeSpace,
                              lsize_, 1>::Multiply(
    LTuple const& left,
    RTuple const& right) -> Result {
  auto const right_head = std::get<0>(right);
  return CartesianProductMultiplicativeSpace<
             decltype(right_head), LTuple>::Multiply(left, right_head);
}

template<typename Scalar, typename Tuple,
         typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
class CartesianProductPointwiseMultiplicativeSpace;

template<typename Scalar, typename Tuple, std::size_t... indices>
class CartesianProductPointwiseMultiplicativeSpace<
    Scalar, Tuple, std::index_sequence<indices...>> {
 public:
  template<typename L, typename R>
  using Product = typename Hilbert<L, R>::InnerProductType;
  template<typename T>
  using ScalarLeftProduct = Product<Scalar, T>;
  template<typename T>
  using ScalarRightProduct = Product<T, Scalar>;

  FORCE_INLINE(static constexpr) Apply<ScalarLeftProduct, Tuple> Multiply(
      Scalar const& left,
      Tuple const& right);
  FORCE_INLINE(static constexpr) Apply<ScalarRightProduct, Tuple> Multiply(
      Tuple const& left,
      Scalar const& right);
};

template<typename Scalar, typename Tuple, std::size_t... indices>
constexpr auto CartesianProductPointwiseMultiplicativeSpace<
    Scalar, Tuple,
    std::index_sequence<indices...>>::Multiply(Scalar const& left,
                                               Tuple const& right)
    -> Apply<ScalarLeftProduct, Tuple> {
  return {Hilbert<Scalar, std::tuple_element_t<indices, Tuple>>::InnerProduct(
      left, std::get<indices>(right))...};
}

template<typename Scalar, typename Tuple, std::size_t... indices>
constexpr auto CartesianProductPointwiseMultiplicativeSpace<
    Scalar, Tuple,
    std::index_sequence<indices...>>::Multiply(Tuple const& left,
                                               Scalar const& right)
    -> Apply<ScalarRightProduct, Tuple> {
  return {Hilbert<std::tuple_element_t<indices, Tuple>, Scalar>::InnerProduct(
      std::get<indices>(left), right)...};
}

// DirectSum.

template<affine... T>
constexpr DirectSum<T...>::DirectSum(T&&... args) : tuple(args...) {}

template<affine... T>
constexpr DirectSum<T...>::DirectSum(std::tuple<T...>&& tuple) : tuple(tuple) {}

template<affine... T>
template<std::size_t I>
constexpr auto const& DirectSum<T...>::get() const {
  return std::get<I>(tuple);
}

template<affine... T>
template<std::size_t I>
constexpr auto& DirectSum<T...>::get() {
  return std::get<I>(tuple);
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
template<affine... U>
DirectSum<T...>& DirectSum<T...>::operator+=(DirectSum<U...> const& right) {
  *this = *this + right;
  return *this;
}

template<affine... T>
template<affine... U>
DirectSum<T...>& DirectSum<T...>::operator-=(DirectSum<U...> const& right) {
  *this = *this - right;
  return *this;
}

template<affine... T>
template<typename Scalar>
DirectSum<T...>& DirectSum<T...>::operator*=(Scalar const& right) {
  *this = *this * right;
  return *this;
}

template<affine... T>
template<typename Scalar>
DirectSum<T...>& DirectSum<T...>::operator/=(Scalar const& right) {
  *this = *this / right;
  return *this;
}

template<affine... T>
constexpr auto operator+(DirectSum<T...> const& right) {
  return right;
}

template<affine... T>
constexpr auto operator-(DirectSum<T...> const& right) {
  std::tuple<> zero;
  return DirectSum(
      CartesianProductSubtraction<decltype(zero), std::tuple<T...>>::Subtract(
          zero, right.tuple));
}

template<affine... L, affine... R>
constexpr auto operator+(DirectSum<L...> const& left,
                         DirectSum<R...> const& right) {
  return DirectSum(
      CartesianProductAddition<std::tuple<L...>, std::tuple<R...>>::Add(
          left.tuple, right.tuple));
}

template<affine... L, affine... R>
constexpr auto operator-(DirectSum<L...> const& left,
                         DirectSum<R...> const& right) {
  return DirectSum(
      CartesianProductSubtraction<std::tuple<L...>, std::tuple<R...>>::Subtract(
          left.tuple, right.tuple));
}

template<typename L, affine... R>
constexpr auto operator*(L const& left, DirectSum<R...> const& right) {
  return DirectSum(CartesianProductVectorSpace<L, std::tuple<R...>>::Multiply(
      left, right.tuple));
}

template<affine... L, typename R>
constexpr auto operator*(DirectSum<L...> const& left, R const& right) {
  return DirectSum(CartesianProductVectorSpace<R, std::tuple<L...>>::Multiply(
      left.tuple, right));
}

template<affine... L, typename R>
constexpr auto operator/(DirectSum<L...> const& left, R const& right) {
  return DirectSum(CartesianProductVectorSpace<R, std::tuple<L...>>::Divide(
      left.tuple, right));
}

template<affine... T>
constexpr auto InnerProduct(DirectSum<T...> const& left,
                            DirectSum<T...> const& right) {
  using T0 = std::tuple_element_t<0, DirectSum<T...>>;
  decltype(std::declval<T0>() * std::declval<T0>()) product = {};
  for_all_of(left.tuple, right.tuple)
      .loop([&product](auto const& leftᵢ, auto const& rightᵢ) {
        product += leftᵢ * rightᵢ;
      });

  return product;
}

}  // namespace internal

namespace vector_space {

template<typename RTuple>
FORCE_INLINE(constexpr)
auto operator+(RTuple const& right) {
  return right;
}

template<typename RTuple>
FORCE_INLINE(constexpr)
auto operator-(RTuple const& right) {
  std::tuple<> zero;
  return internal::CartesianProductSubtraction<decltype(zero), RTuple>::
      Subtract(zero, right);
}

template<typename LTuple, typename RTuple>
FORCE_INLINE(constexpr)
auto operator+(LTuple const& left, RTuple const& right) {
  return internal::CartesianProductAddition<LTuple, RTuple>::
      Add(left, right);
}

template<typename LTuple, typename RTuple>
FORCE_INLINE(constexpr)
auto operator-(LTuple const& left, RTuple const& right) {
  return internal::CartesianProductSubtraction<LTuple, RTuple>::
      Subtract(left, right);
}

template<typename Scalar, typename Tuple, typename, typename>
FORCE_INLINE(constexpr)
auto operator*(Scalar const& left, Tuple const& right) {
  return internal::CartesianProductVectorSpace<Scalar, Tuple>::
      Multiply(left, right);
}

template<typename Tuple, typename Scalar, typename, typename, typename>
FORCE_INLINE(constexpr)
auto operator*(Tuple const& left, Scalar const& right) {
  return internal::CartesianProductVectorSpace<Scalar, Tuple>::
      Multiply(left, right);
}

template<typename Scalar, typename Tuple>
FORCE_INLINE(constexpr)
auto operator/(Tuple const& left, Scalar const& right) {
  return internal::CartesianProductVectorSpace<Scalar, Tuple>::
      Divide(left, right);
}

}  // namespace vector_space

namespace polynomial_ring {

template<typename LTuple, typename RTuple, typename, typename>
FORCE_INLINE(constexpr)
auto operator*(LTuple const& left, RTuple const& right) {
  return internal::PolynomialRing<
      LTuple, RTuple,
      internal::CartesianProductVectorSpace>::Multiply(left, right);
}

template<int exponent, typename Tuple>
constexpr auto Pow(Tuple const& tuple) {
  static_assert(exponent > 0, "Cannot raise a tuple to the zero-th power");
  if constexpr (exponent == 1) {
    return tuple;
  } else if constexpr (exponent % 2 == 0) {
    return Pow<exponent / 2>(tuple * tuple);
  } else {
    return Pow<exponent / 2>(tuple * tuple) * tuple;
  }
}

}  // namespace polynomial_ring

namespace pointwise_inner_product {

template<typename LTuple, typename RTuple,
         typename, typename>
constexpr auto PointwiseInnerProduct(LTuple const& left, RTuple const& right) {
  return internal::PolynomialRing<
      LTuple, RTuple,
      internal::
          CartesianProductPointwiseMultiplicativeSpace>::Multiply(left, right);
}

}  // namespace pointwise_inner_product
}  // namespace _cartesian_product
}  // namespace geometry
}  // namespace principia
