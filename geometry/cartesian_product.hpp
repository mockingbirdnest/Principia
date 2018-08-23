
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
// types must have an addition and a subtraction.
template<typename LTuple,
         typename RTuple,
         typename = std::make_integer_sequence<
             int,
             std::max(std::tuple_size_v<LTuple>, std::tuple_size_v<RTuple>)>>
struct CartesianProductAdditiveGroup;

// This struct exports three functions:
//
//   static constexpr ... Multiply(Scalar const& left, Tuple const& right);
//   static constexpr ... Multiply(Tuple const& left, Scalar const& right);
//   static constexpr ... Divide(Tuple const& left, Scalar const& right);
//
// where the return types are appropriate for the operation.  The element types
// of the tuple must have a structure of vector space.
template<typename Scalar,
         typename Tuple,
         typename = std::make_integer_sequence<int, std::tuple_size_v<Tuple>>>
struct CartesianProductVectorSpace;

//TODO(phl):comment
template<typename LTuple, typename RTuple,
         int lsize_ = std::tuple_size_v<LTuple>,
         int rsize_ = std::tuple_size_v<RTuple>>
struct CartesianProductAlgebra;

}  // namespace internal_cartesian_product

using internal_cartesian_product::CartesianProductAdditiveGroup;
using internal_cartesian_product::CartesianProductAlgebra;
using internal_cartesian_product::CartesianProductVectorSpace;

}  // namespace geometry
}  // namespace principia

#include "geometry/cartesian_product_body.hpp"
