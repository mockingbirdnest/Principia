
#include "geometry/cartesian_product.hpp"

#include <algorithm>
#include <type_traits>

#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace geometry {
namespace internal_cartesian_product {

using quantities::Apply;
using quantities::Apply2;

template<typename LTuple, typename RTuple, std::size_t... indices>
class CartesianProductAdditiveGroup<LTuple, RTuple,
                                     std::index_sequence<indices...>> {
  // The types of the result of addition and subtraction, with suitable
  // specializations for the void case of Apply2.
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
  static constexpr Apply2<Sum, LTuple, RTuple> Add(
      LTuple const& left,
      RTuple const& right);
  static constexpr Apply2<Difference, LTuple, RTuple> Subtract(
      LTuple const& left,
      RTuple const& right);

 private:
  // Utilities for adding/subtracting elements at the given index.
  template<std::size_t index>
  static constexpr std::tuple_element_t<index, Apply2<Sum, LTuple, RTuple>>
  AddElement(LTuple const& left, RTuple const& right);
  template<std::size_t index>
  static constexpr std::tuple_element_t<index,
                                        Apply2<Difference, LTuple, RTuple>>
  SubtractElement(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, std::size_t... indices>
constexpr auto CartesianProductAdditiveGroup<
    LTuple, RTuple,
    std::index_sequence<indices...>>::Add(LTuple const& left,
                                          RTuple const& right)
    -> Apply2<Sum, LTuple, RTuple> {
  return {AddElement<indices>(left, right)...};
}

template<typename LTuple, typename RTuple, std::size_t... indices>
constexpr auto CartesianProductAdditiveGroup<
    LTuple, RTuple,
    std::index_sequence<indices...>>::Subtract(LTuple const& left,
                                               RTuple const& right)
    -> Apply2<Difference, LTuple, RTuple> {
  return {SubtractElement<indices>(left, right)...};
}

template<typename LTuple, typename RTuple, std::size_t... indices>
template<std::size_t index>
constexpr auto
CartesianProductAdditiveGroup<LTuple, RTuple, std::index_sequence<indices...>>::
AddElement(LTuple const& left, RTuple const& right)
    -> std::tuple_element_t<index, Apply2<Sum, LTuple, RTuple>> {
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
constexpr auto
CartesianProductAdditiveGroup<LTuple, RTuple, std::index_sequence<indices...>>::
SubtractElement(LTuple const& left, RTuple const& right)
    -> std::tuple_element_t<index, Apply2<Difference, LTuple, RTuple>> {
  if constexpr (index < std::min(std::tuple_size_v<LTuple>,
                                 std::tuple_size_v<RTuple>)) {
    return std::get<index>(left) - std::get<index>(right);
  } else if constexpr (index < std::tuple_size_v<LTuple>) {
    return std::get<index>(left);
  } else {
    return -std::get<index>(right);
  }
}

template<typename Scalar, typename Tuple, int... indices>
struct CartesianProductVectorSpace<Scalar,
                                   Tuple,
                                   std::integer_sequence<int, indices...>> {
  template<typename T>
  using ScalarLeftProduct = quantities::Product<Scalar, T>;
  template<typename T>
  using ScalarRightProduct = quantities::Product<T, Scalar>;
  template<typename T>
  using Quotient = quantities::Quotient<T, Scalar>;

  static constexpr Apply<ScalarLeftProduct, Tuple> Multiply(
      Scalar const& left,
      Tuple const& right);
  static constexpr Apply<ScalarRightProduct, Tuple> Multiply(
      Tuple const& left,
      Scalar const& right);

  static constexpr Apply<Quotient, Tuple> Divide(
      Tuple const& left,
      Scalar const& right);
};

template<typename Scalar, typename Tuple, int... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::integer_sequence<int, indices...>>::Multiply(Scalar const& left,
                                                      Tuple const& right)
    -> Apply<ScalarLeftProduct, Tuple> {
  return {left * std::get<indices>(right)...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::integer_sequence<int, indices...>>::Multiply(Tuple const& left,
                                                      Scalar const& right)
    -> Apply<ScalarRightProduct, Tuple> {
  return {std::get<indices>(left) * right...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::integer_sequence<int, indices...>>::Divide(Tuple const& left,
                                                    Scalar const& right)
    -> Apply<Quotient, Tuple> {
  return {std::get<indices>(left) / right...};
}

template<typename Element, typename Tuple,
         typename = std::make_integer_sequence<int, std::tuple_size_v<Tuple>>>
struct ConsGenerator;

template<typename Element, typename Tuple, int... indices>
struct ConsGenerator<Element, Tuple,
                     std::integer_sequence<int, indices...>> {
  using Type = std::tuple<Element, std::tuple_element_t<indices, Tuple>...>;
  static constexpr Type Cons(Element const& element, Tuple const& tuple);
};

template<typename Element, typename Tuple, int... indices>
constexpr auto
ConsGenerator<Element, Tuple, std::integer_sequence<int, indices...>>::
Cons(Element const& element, Tuple const& tuple) -> Type {
  return std::make_tuple(element, std::get<indices>(tuple)...);
}

template<typename Tuple,
         typename = std::make_integer_sequence<int, std::tuple_size_v<Tuple>>>
struct TailGenerator;

template<typename Tuple, int head_index, int... tail_indices>
struct TailGenerator<Tuple,
                     std::integer_sequence<int, head_index, tail_indices...>> {
  using Type = std::tuple<std::tuple_element_t<tail_indices, Tuple>...>;
  static constexpr Type Tail(Tuple const& tuple);
};

template<typename Tuple, int head_index, int... tail_indices>
constexpr auto
TailGenerator<Tuple, std::integer_sequence<int, head_index, tail_indices...>>::
Tail(Tuple const& tuple) -> Type {
  return std::make_tuple(std::get<tail_indices>(tuple)...);
}

template<typename LTuple, typename RTuple,
         int lsize_ = std::tuple_size_v<LTuple>,
         int rsize_ = std::tuple_size_v<RTuple>>
struct CartesianProductAlgebra {
  // Right is split into head (index 0) and tail (the rest).  The tail is a
  // polynomial with valuation 1.
  using RHead = std::tuple_element_t<0, RTuple>;
  using RTail = typename TailGenerator<RTuple>::Type;

  // To implement the polynomial multiplication left * right_tail, we need to
  // insert a zero for the lowest degree (because of the valuation 1).  This is
  // the type of that zero.
  using LHead = std::tuple_element_t<0, LTuple>;
  using Zero = quantities::Product<LHead, RHead>;

  using LTupleRHeadProduct =
      decltype(CartesianProductVectorSpace<RHead, LTuple>::Multiply(
          std::declval<LTuple>(),
          std::declval<RHead>()));
  using LTupleRTailProduct =
      decltype(CartesianProductAlgebra<LTuple, RTail>::Mult(
          std::declval<LTuple>(),
          std::declval<RTail>()));
  using ZeroLTupleRTailProduct =
      typename ConsGenerator<Zero, LTupleRTailProduct>::Type;

  using Result =
      decltype(
          CartesianProductAdditiveGroup<LTupleRHeadProduct,
                                        ZeroLTupleRTailProduct>::
              Add(std::declval<LTupleRHeadProduct>(),
                  std::declval<ZeroLTupleRTailProduct>()));
  static constexpr Result Mult(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, int lsize_>
struct CartesianProductAlgebra<LTuple, RTuple, lsize_, 1> {
  using RHead = std::tuple_element_t<0, RTuple>;
  using Result = decltype(CartesianProductVectorSpace<RHead, LTuple>::Multiply(
      std::declval<LTuple>(),
      std::declval<RHead>()));
  static constexpr Result Mult(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, int lsize_, int rsize_>
constexpr auto CartesianProductAlgebra<LTuple, RTuple, lsize_, rsize_>::Mult(
    LTuple const& left,
    RTuple const& right) -> Result {
  RHead const right_head = std::get<0>(right);
  auto const right_tail = TailGenerator<RTuple>::Tail(right);

  return CartesianProductAdditiveGroup<LTupleRHeadProduct,
                                       ZeroLTupleRTailProduct>::Add(
             CartesianProductVectorSpace<RHead, LTuple>::Multiply(left,
                                                                  right_head),
             ConsGenerator<Zero, LTupleRTailProduct>::Cons(
                 Zero{}, 
                 CartesianProductAlgebra<LTuple, RTail>::Mult(left,
                                                              right_tail)));
}

template<typename LTuple, typename RTuple, int lsize_>
constexpr auto CartesianProductAlgebra<LTuple, RTuple, lsize_, 1>::Mult(
    LTuple const& left,
    RTuple const& right) -> Result {
  RHead const right_head = std::get<0>(right);
  return CartesianProductVectorSpace<RHead, LTuple>::Multiply(left, right_head);
}

}  // namespace internal_cartesian_product
}  // namespace geometry
}  // namespace principia
