#include "geometry/cartesian_product.hpp"

#include <algorithm>
#include <tuple>
#include <type_traits>

#include "base/macros.hpp"
#include "geometry/hilbert.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"

// The use of FORCE_INLINE in this file is because we want the construction of
// complex polynomials to be reasonably efficient.

namespace principia {
namespace geometry {
namespace _cartesian_product {
namespace internal {

using quantities::Apply;

template<typename LTuple, typename RTuple,
         typename = std::make_index_sequence<
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
class CartesianProductAdditiveGroup;

template<typename LTuple, typename RTuple, std::size_t... indices>
class CartesianProductAdditiveGroup<LTuple, RTuple,
                                    std::index_sequence<indices...>> {
  // The types of the result of addition and subtraction, with suitable
  // specializations for the void case of Apply.
  template<typename L, typename R>
  struct TypesGenerator {
    using Sum = quantities::Sum<L, R>;
    using Difference = quantities::Difference<L, R>;
  };
  template<typename L>
  struct TypesGenerator<L, void> {
    using Sum = L;
    using Difference = L;
  };
  template<typename R>
  struct TypesGenerator<void, R> {
    using Sum = R;
    using Difference = R;
  };

  // Aliases for use as the transform in Apply2.
  template<typename L, typename R>
  using Sum = typename TypesGenerator<L, R>::Sum;
  template<typename L, typename R>
  using Difference = typename TypesGenerator<L, R>::Difference;

 public:
  FORCE_INLINE(static constexpr) Apply<Sum, LTuple, RTuple> Add(
      LTuple const& left,
      RTuple const& right);
  FORCE_INLINE(static constexpr) Apply<Difference, LTuple, RTuple> Subtract(
      LTuple const& left,
      RTuple const& right);

 private:
  // Utilities for adding/subtracting elements at the given index.
  template<std::size_t index>
  static constexpr std::tuple_element_t<index, Apply<Sum, LTuple, RTuple>>
  AddElement(LTuple const& left, RTuple const& right);
  template<std::size_t index>
  static constexpr std::tuple_element_t<index,
                                        Apply<Difference, LTuple, RTuple>>
  SubtractElement(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, std::size_t... indices>
constexpr auto CartesianProductAdditiveGroup<
    LTuple, RTuple,
    std::index_sequence<indices...>>::Add(LTuple const& left,
                                          RTuple const& right)
    -> Apply<Sum, LTuple, RTuple> {
  return {AddElement<indices>(left, right)...};
}

template<typename LTuple, typename RTuple, std::size_t... indices>
constexpr auto CartesianProductAdditiveGroup<
    LTuple, RTuple,
    std::index_sequence<indices...>>::Subtract(LTuple const& left,
                                               RTuple const& right)
    -> Apply<Difference, LTuple, RTuple> {
  return {SubtractElement<indices>(left, right)...};
}

template<typename LTuple, typename RTuple, std::size_t... indices>
template<std::size_t index>
constexpr auto CartesianProductAdditiveGroup<
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
template<std::size_t index>
constexpr auto CartesianProductAdditiveGroup<
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
  using Product = quantities::Product<L, R>;
  template<typename T>
  using ScalarLeftProduct = Product<Scalar, T>;
  template<typename T>
  using ScalarRightProduct = Product<T, Scalar>;
  template<typename T>
  using Quotient = quantities::Quotient<T, Scalar>;

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

}  // namespace internal
}  // namespace _cartesian_product

namespace _polynomial_ring {
namespace internal {

using _cartesian_product::internal::CartesianProductAdditiveGroup;

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

// Technically, this is only a ring if |CartesianProductMultiplicativeSpace| is
// the bona fide |CartesianProductVectorSpace|.  If it is not, this implements a
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
      decltype(CartesianProductAdditiveGroup<LTupleRHeadProduct,
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
  using _cartesian_product::operator+;

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

}  // namespace internal
}  // namespace _polynomial_ring

namespace _pointwise_inner_product {
namespace internal {

using quantities::Apply;
using _cartesian_product::internal::CartesianProductAdditiveGroup;
using namespace principia::geometry::_hilbert;

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

}  // namespace internal
}  // namespace _pointwise_inner_product

namespace _cartesian_product {
namespace internal {

template<typename RTuple>
FORCE_INLINE(constexpr) auto operator+(RTuple const& right) {
  return right;
}

template<typename RTuple>
FORCE_INLINE(constexpr)
auto operator-(RTuple const& right) {
  std::tuple<> zero;
  return CartesianProductAdditiveGroup<decltype(zero), RTuple>::Subtract(zero,
                                                                         right);
}

template<typename LTuple, typename RTuple>
FORCE_INLINE(constexpr)
auto operator+(LTuple const& left, RTuple const& right) {
  return CartesianProductAdditiveGroup<LTuple, RTuple>::Add(left, right);
}

template<typename LTuple, typename RTuple>
FORCE_INLINE(constexpr)
auto operator-(LTuple const& left, RTuple const& right) {
  return CartesianProductAdditiveGroup<LTuple, RTuple>::Subtract(left, right);
}

template<typename Scalar, typename Tuple, typename, typename>
FORCE_INLINE(constexpr)
auto operator*(Scalar const& left, Tuple const& right) {
  return CartesianProductVectorSpace<Scalar, Tuple>::Multiply(left, right);
}

template<typename Tuple, typename Scalar, typename, typename, typename>
FORCE_INLINE(constexpr)
auto operator*(Tuple const& left, Scalar const& right) {
  return CartesianProductVectorSpace<Scalar, Tuple>::Multiply(left, right);
}

template<typename Scalar, typename Tuple>
FORCE_INLINE(constexpr)
auto operator/(Tuple const& left, Scalar const& right) {
  return CartesianProductVectorSpace<Scalar, Tuple>::Divide(left, right);
}

}  // namespace internal
}  // namespace _cartesian_product

namespace _polynomial_ring {
namespace internal {

template<typename LTuple, typename RTuple, typename, typename>
FORCE_INLINE(constexpr)
auto operator*(LTuple const& left, RTuple const& right) {
  return _polynomial_ring::internal::PolynomialRing<
      LTuple, RTuple,
      _cartesian_product::internal::CartesianProductVectorSpace>::Multiply(
          left, right);
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

}  // namespace internal
}  // namespace _polynomial_ring

namespace _pointwise_inner_product {
namespace internal {

template<typename LTuple, typename RTuple,
         typename, typename>
constexpr auto PointwiseInnerProduct(LTuple const& left, RTuple const& right) {
  return _polynomial_ring::internal::PolynomialRing<
      LTuple, RTuple,
      _pointwise_inner_product::internal::
          CartesianProductPointwiseMultiplicativeSpace>::Multiply(left, right);
}

}  // namespace internal
}  // namespace _pointwise_inner_product

}  // namespace geometry
}  // namespace principia
