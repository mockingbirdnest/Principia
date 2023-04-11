#include "geometry/traits.hpp"

#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "physics/inertia_tensor.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::geometry::_traits;
using namespace principia::physics::_inertia_tensor;
using namespace principia::quantities::_named_quantities;

TEST(Traits, IsVectorV) {
  static_assert(is_vector_v<int>);
  static_assert(is_vector_v<double>);
  static_assert(is_vector_v<Area>);
  static_assert(is_vector_v<Frequency const>);
  // Not sure if the following is what we want, but at least let's nail it in a
  // test.
  static_assert(!is_vector_v<Entropy&>);
  static_assert(!is_vector_v<float const&>);

  using F = Frame<struct FTag>;
  static_assert(is_vector_v<Displacement<F>>);
  static_assert(is_vector_v<Velocity<F>>);
  static_assert(is_vector_v<AngularVelocity<F>>);
  static_assert(is_vector_v<InertiaTensor<F>>);
  static_assert(is_vector_v<Velocity<F> const>);
  static_assert(!is_vector_v<Instant>);
  static_assert(!is_vector_v<Position<F>>);
  // Same as above.
  static_assert(!is_vector_v<Velocity<F> const&>);
}

}  // namespace geometry
}  // namespace principia
