
#include "geometry/cartesian_product.hpp"

#include <algorithm>
#include <tuple>
#include <type_traits>

#include "base/macros.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"

// The use of FORCE_INLINE in this file is because we want the construction of
// complex polynomials to be reasonably efficient.

namespace principia {
namespace geometry {
namespace internal_cartesian_product {

using quantities::Apply;

template<typename LTuple, typename RTuple,
         typename = std::make_index_sequence<
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
struct CartesianProductAdditiveGroup;

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
struct CartesianProductVectorSpace;

template<typename Scalar, typename Tuple, std::size_t... indices>
class CartesianProductVectorSpace<Scalar, Tuple,
                                  std::index_sequence<indices...>> {
  template<typename T>
  using ScalarLeftProduct = quantities::Product<Scalar, T>;
  template<typename T>
  using ScalarRightProduct = quantities::Product<T, Scalar>;
  template<typename T>
  using Quotient = quantities::Quotient<T, Scalar>;

 public:
  FORCE_INLINE(static constexpr) Apply<ScalarLeftProduct, Tuple> Multiply(
      Scalar const& left,
      Tuple const& right);
  FORCE_INLINE(static constexpr) Apply<ScalarRightProduct, Tuple> Multiply(
      Tuple const& left,
      Scalar const& right);

  FORCE_INLINE(static constexpr) Apply<Quotient, Tuple> Divide(
      Tuple const& left,
      Scalar const& right);
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
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::index_sequence<indices...>>::Divide(Tuple const& left,
                                             Scalar const& right)
    -> Apply<Quotient, Tuple> {
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

template<typename LTuple, typename RTuple,
         int lsize_ = std::tuple_size_v<LTuple>,
         int rsize_ = std::tuple_size_v<RTuple>>
class PolynomialRing {
  // Right is split into head (index 0) and tail (the rest).  The tail is a
  // polynomial with valuation 1.
  using RHead = std::tuple_element_t<0, RTuple>;
  using RTail = typename TailGenerator<RTuple>::Type;

  // To implement the polynomial multiplication left * right_tail, we need to
  // insert a zero for the lowest degree (because of the valuation 1).  This is
  // the type of that zero.
  using LHead = std::tuple_element_t<0, LTuple>;
  using Zero = quantities::Product<LHead, RHead>;

  // The hard part: generating the type of the result.
  using LTupleRHeadProduct =
      decltype(CartesianProductVectorSpace<RHead, LTuple>::Multiply(
          std::declval<LTuple>(),
          std::declval<RHead>()));
  using LTupleRTailProduct =
      decltype(PolynomialRing<LTuple, RTail>::Multiply(
          std::declval<LTuple>(),
          std::declval<RTail>()));
  using ZeroLTupleRTailProduct =
      typename ConsGenerator<Zero, LTupleRTailProduct>::Type;

  using Result =
      decltype(CartesianProductAdditiveGroup<LTupleRHeadProduct,
                                             ZeroLTupleRTailProduct>::
                   Add(std::declval<LTupleRHeadProduct>(),
                       std::declval<ZeroLTupleRTailProduct>()));

 public:
  FORCE_INLINE(static constexpr)
  Result Multiply(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, int lsize_>
class PolynomialRing<LTuple, RTuple, lsize_, 1> {
  using RHead = std::tuple_element_t<0, RTuple>;
  using Result = decltype(CartesianProductVectorSpace<RHead, LTuple>::Multiply(
      std::declval<LTuple>(),
      std::declval<RHead>()));

 public:
  FORCE_INLINE(static constexpr)
  Result Multiply(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, int lsize_, int rsize_>
constexpr auto PolynomialRing<LTuple, RTuple, lsize_, rsize_>::Multiply(
    LTuple const& left,
    RTuple const& right) -> Result {
  using cartesian_product::operator+;
  using cartesian_product::operator*;
  using polynomial_ring::operator*;

  auto const right_head = std::get<0>(right);
  auto const right_tail = TailGenerator<RTuple>::Tail(right);

  return left * right_head +
         ConsGenerator<Zero, LTupleRTailProduct>::Cons(
             Zero{}, left * right_tail);
}

template<typename LTuple, typename RTuple, int lsize_>
constexpr auto PolynomialRing<LTuple, RTuple, lsize_, 1>::Multiply(
    LTuple const& left,
    RTuple const& right) -> Result {
  using cartesian_product::operator*;

  auto const right_head = std::get<0>(right);
  return left * right_head;
}

}  // namespace internal_cartesian_product

namespace cartesian_product {

template<typename LTuple, typename RTuple>
FORCE_INLINE(constexpr)
auto operator+(LTuple const& left, RTuple const& right) {
  return internal_cartesian_product::
      CartesianProductAdditiveGroup<LTuple, RTuple>::Add(left, right);
}

template<typename LTuple, typename RTuple>
FORCE_INLINE(constexpr)
auto operator-(LTuple const& left, RTuple const& right) {
  return internal_cartesian_product::
      CartesianProductAdditiveGroup<LTuple, RTuple>::Subtract(left, right);
}

template<typename Scalar, typename Tuple, typename, typename>
FORCE_INLINE(constexpr) auto operator*(Scalar const& left, Tuple const& right) {
  return internal_cartesian_product::
      CartesianProductVectorSpace<Scalar, Tuple>::Multiply(left, right);
}

template<typename Tuple, typename Scalar, typename, typename, typename>
FORCE_INLINE(constexpr) auto operator*(Tuple const& left, Scalar const& right) {
  return internal_cartesian_product::
      CartesianProductVectorSpace<Scalar, Tuple>::Multiply(left, right);
}

template<typename Scalar, typename Tuple>
FORCE_INLINE(constexpr) auto operator/(Tuple const& left, Scalar const& right) {
  return internal_cartesian_product::
      CartesianProductVectorSpace<Scalar, Tuple>::Divide(left, right);
}

}  // namespace cartesian_product

namespace polynomial_ring {

template<typename LTuple, typename RTuple, typename, typename>
FORCE_INLINE(constexpr)
auto operator*(LTuple const& left, RTuple const& right) {
  return internal_cartesian_product::
      PolynomialRing<LTuple, RTuple>::Multiply(left, right);
}

}  // namespace polynomial_ring

}  // namespace geometry
}  // namespace principia
