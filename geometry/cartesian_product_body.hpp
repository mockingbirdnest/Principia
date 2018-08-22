
#include "geometry/cartesian_product.hpp"

#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace geometry {
namespace internal_cartesian_product {

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

}  // namespace internal_cartesian_product
}  // namespace geometry
}  // namespace principia
