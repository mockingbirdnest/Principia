
#pragma once

#include <algorithm>
#include <tuple>
#include <utility>

namespace principia {
namespace geometry {
namespace internal_cartesian_product {

// This struct exports two functions:
//
//   static constexpr ... Add(LTuple const& left, RTuple const& right);
//   static constexpr ... Subtract(LTuple const& left, RTuple const& right);
//
// where the return types are appropriate for the operation.  The two tuples may
// be of different types (including different sizes) but corresponding element
// types must have an addition and a subtraction.  In case of different sizes,
// the smaller tuple is zero-extended, where zero is the neutral element of
// addition on the relevant types.
template<typename LTuple,
         typename RTuple,
         typename = std::make_integer_sequence<
             int,
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
struct CartesianProductAdditiveGroup;

}  // namespace internal_cartesian_product

using internal_cartesian_product::CartesianProductAdditiveGroup;

}  // namespace geometry
}  // namespace principia

#include "geometry/cartesian_product_body.hpp"
