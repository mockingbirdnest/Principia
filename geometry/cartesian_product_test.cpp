#include "geometry/cartesian_product.hpp"

#include <tuple>

#include "base/algebra.hpp"
#include "geometry/point.hpp"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {

using namespace principia::base::_algebra;
using namespace principia::geometry::_cartesian_product;
using namespace principia::geometry::_point;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_named_quantities;

TEST(CartesianProductTest, Concepts) {
  static_assert(real_vector_space<VectorTuple<std::tuple<double>>>);
  static_assert(real_vector_space<VectorTuple<std::tuple<double, double>>>);
  static_assert(real_vector_space<VectorTuple<std::tuple<Length>>>);
  static_assert(real_vector_space<VectorTuple<std::tuple<Length, Mass, Time>>>);

  static_assert(real_affine_space<VectorTuple<std::tuple<Point<Length>>>>);
  static_assert(
      real_affine_space<VectorTuple<std::tuple<Point<Length>, Speed>>>);
}

TEST(CartesianProductTest, FixedVector) {
  FixedVector<VectorTuple<std::tuple<double, double>>, 1>(
      std::array<VectorTuple<std::tuple<double, double>>, 1>{
          VectorTuple<std::tuple<double, double>>{.tuple = {1, 2}}});
}

}  // namespace geometry
}  // namespace principia
