
#include "geometry/cartesian_product.hpp"

#include <algorithm>

#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace geometry {
namespace internal_cartesian_product {

using quantities::Apply;
using quantities::Apply2;

template<typename LTuple, typename RTuple, int... indices>
struct CartesianProductAdditiveGroup<LTuple, RTuple,
                                     std::integer_sequence<int, indices...>> {
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

  static constexpr Apply2<Sum, LTuple, RTuple> Add(
      LTuple const& left,
      RTuple const& right);
  static constexpr Apply2<Difference, LTuple, RTuple> Subtract(
      LTuple const& left,
      RTuple const& right);
};

template<typename LTuple, typename RTuple, int... indices>
constexpr auto CartesianProductAdditiveGroup<
    LTuple, RTuple,
    std::integer_sequence<int, indices...>>::Add(LTuple const& left,
                                                 RTuple const& right)
    -> Apply2<Sum, LTuple, RTuple> {
  return {(
      indices < std::min(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)
          ? std::get<indices>(left) + std::get<indices>(right)
          : indices < std::tuple_size_v<LTuple>
                ? std::get<indices>(left)
                : std::get<indices>(right))...};
}

template<typename LTuple, typename RTuple, int... indices>
constexpr auto CartesianProductAdditiveGroup<
    LTuple, RTuple,
    std::integer_sequence<int, indices...>>::Subtract(LTuple const& left,
                                                      RTuple const& right)
    -> Apply2<Difference, LTuple, RTuple> {
  return {
      (indices < std::min(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)
           ? std::get<indices>(left) - std::get<indices>(right)
           : indices < std::tuple_size_v<LTuple>
                 ? std::get<indices>(left)
                 : -std::get<indices>(right))...};
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

  static constexpr typename Apply<ScalarLeftProduct, Tuple> Multiply(
      Scalar const& left,
      Tuple const& right);
  static constexpr typename Apply<ScalarRightProduct, Tuple> Multiply(
      Tuple const& left,
      Scalar const& right);

  static constexpr typename Apply<Quotient, Tuple> Divide(
      Tuple const& left,
      Scalar const& right);
};

template<typename Scalar, typename Tuple, int... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::integer_sequence<int, indices...>>::Multiply(Scalar const& left,
                                                      Tuple const& right)
    -> typename Apply<ScalarLeftProduct, Tuple> {
  return {left * std::get<indices>(right)...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::integer_sequence<int, indices...>>::Multiply(Tuple const& left,
                                                      Scalar const& right)
    -> typename Apply<ScalarRightProduct, Tuple> {
  return {std::get<indices>(left) * right...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr auto CartesianProductVectorSpace<
    Scalar, Tuple,
    std::integer_sequence<int, indices...>>::Divide(Tuple const& left,
                                                    Scalar const& right)
    -> typename Apply<Quotient, Tuple> {
  return {std::get<indices>(left) / right...};
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

template<typename LTuple, typename RTuple, int lsize_, int rsize_>
struct CartesianProductAlgebra {
  using RHead = std::element_type_t<0, RTuple>;
  using RTail = typename TailGenerator<LTuple>::Type;
  using LTupleRHeadProduct =
      decltype(CartesianProductVectorSpace<LTuple, RHead>::Multiply(
          std::declval<LTuple>(),
          std::declval<RHead>()));
  using LTupleRTailProduct =
      decltype(CartesianProductVectorSpace<LTuple, RTail>::Multiply(
          std::declval<LTuple>(),
          std::declval<RTail>()));
  using Result = int;
  constexpr Result Mult(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, int lsize_>
struct CartesianProductAlgebra<LTuple, RTuple, lsize_, 1> {
  using RHead = std::element_type_t<0, RTuple>;
  using LTupleRHeadProduct =
      decltype(CartesianProductVectorSpace<LTuple, RHead>::Multiply(
          std::declval<LTuple>(),
          std::declval<RHead>()));
  constexpr LTupleRHeadProduct Mult(LTuple const& left, RTuple const& right);
};

template<typename LTuple, typename RTuple, int lsize_, int rsize_>
constexpr auto CartesianProductAlgebra<LTuple, RTuple, lsize_, rsize_>::Mult(
    LTuple const& left,
    RTuple const& right) -> Result {
  // Right is split into head (index 0) and tail (the rest).  The tail is a
  // polynomial with valuation 1.
  RHead const right_head = std::get<0>(right);
  auto const right_tail = TailGenerator<RTuple>::Tail(right);

  // To implement the polynomial multiplication left * right_tail, we need to
  // insert a zero for the lowest degree (because of the valuation 1).  This is
  // the type of that zero.
  using Zero = quantities::Product<std::tuple_element_t<0, LTuple>,
                                   std::tuple_element_t<0, RTuple>>;
  return CartesianProductVectorSpace<LTuple, RHead>::Multiply(left,
                                                                 right_head) +
         std::tuple_cat({Zero{}}, Mult(left, right_tail));
}

template<typename LTuple, typename RTuple, int lsize_>
constexpr auto CartesianProductAlgebra<LTuple, RTuple, lsize_, 1>::Mult(
    LTuple const& left,
    RTuple const& right) -> LTupleRHeadProduct {
  RHead const right_head = std::get<0>(right);
  return CartesianProductVectorSpace<LTuple, RHead>::Multiply(left, right_head);
}

auto a = std::make_tuple(1.0, true);
auto b = std::make_tuple("foo", 2);
auto c = CartesianProductAlgebra<decltype(a), decltype(b)>::Mult(a, b);

}  // namespace internal_cartesian_product
}  // namespace geometry
}  // namespace principia
