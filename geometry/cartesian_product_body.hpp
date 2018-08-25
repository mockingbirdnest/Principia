
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

template<typename LTuple, typename RTuple,
         typename = std::make_integer_sequence<
             int,
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
struct CartesianProductAdditiveGroup;

template<typename LTuple, typename RTuple, int... indices>
struct CartesianProductAdditiveGroup<LTuple, RTuple,
                                     std::integer_sequence<int, indices...>> {
  // The types of the result of addition and subtraction, with suitable
  // checks for the void case of Apply2.
  template<typename L, typename R>
  using Sum =
      std::conditional_t<std::is_void_v<L>, R,
                         std::conditional_t<std::is_void_v<R>, L,
                                            quantities::Sum<L, R>>>;
  template<typename L, typename R>
  using Difference =
      std::conditional_t<std::is_void_v<L>, R,
                         std::conditional_t<std::is_void_v<R>, L,
                                            quantities::Difference<L, R>>>;

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

template<typename Scalar, typename Tuple,
         typename = std::make_integer_sequence<int, std::tuple_size_v<Tuple>>>
struct CartesianProductVectorSpace;

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

}  // namespace internal_cartesian_product

namespace cartesian_product {

template<typename LTuple, typename RTuple>
constexpr auto operator+(LTuple const& left, RTuple const& right) {
  return internal_cartesian_product::
      CartesianProductAdditiveGroup<LTuple, RTuple>::Add(left, right);
}

template<typename LTuple, typename RTuple>
constexpr auto operator-(LTuple const& left, RTuple const& right) {
  return internal_cartesian_product::
      CartesianProductAdditiveGroup<LTuple, RTuple>::Subtract(left, right);
}

template<typename Scalar, typename Tuple, typename>
constexpr auto operator*(Scalar const& left, Tuple const& right) {
  return internal_cartesian_product::
      CartesianProductVectorSpace<Scalar, Tuple>::Multiply(left, right);
}

template<typename Tuple, typename Scalar, typename, typename>
constexpr auto operator*(Tuple const& left, Scalar const& right) {
  return internal_cartesian_product::
      CartesianProductVectorSpace<Scalar, Tuple>::Multiply(left, right);
}

template<typename Scalar, typename Tuple>
constexpr auto operator/(Tuple const& left, Scalar const& right) {
  return internal_cartesian_product::
      CartesianProductVectorSpace<Scalar, Tuple>::Divide(left, right);
}

}  // namespace cartesian_product
}  // namespace geometry
}  // namespace principia
