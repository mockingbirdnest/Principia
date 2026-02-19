#include "geometry/hilbert.hpp"

#include <type_traits>

#include "base/algebra.hpp"
#include "geometry/complexification.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/algebra.hpp"

namespace principia {
namespace geometry {

using namespace principia::base::_algebra;
using namespace principia::geometry::_complexification;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_algebra;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Ne;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

class HilbertTest : public ::testing::Test {};

TEST_F(HilbertTest, Fields) {
  static_assert(hilbert<double>);
  static_assert(real_dimension<double> == 1);
  static_assert(hilbert<Length>);
  static_assert(real_dimension<Length> == 1);
  EXPECT_THAT(InnerProduct(5 * Metre, 5 * Metre), Eq(25 * Pow<2>(Metre)));
  // Not an ordered field.
  static_assert(field<IntegerModulo<2>>);
  static_assert(!std::totally_ordered<IntegerModulo<2>>);
  static_assert(real_dimension<Length> == 1);
  static_assert(!hilbert<IntegerModulo<2>>);
  // Not an ordered field either, but has an explicit InnerProduct (which
  // differs from operator*).
  static_assert(field<Complexification<double>>);
  static_assert(!std::totally_ordered<Complexification<double>>);
  static_assert(hilbert<Complexification<double>>);
  const Complexification<double> i(0, 1);
  const auto z = 5 * Metre + 5 * Metre * i;
  EXPECT_THAT(InnerProduct(z, z), AllOf(Eq(50 * Pow<2>(Metre)), Ne(z * z)));
  EXPECT_THAT(z * z, Eq(50 * Pow<2>(Metre) * i));
}

TEST_F(HilbertTest, Space) {
  // No InnerProduct on R3Element.
  static_assert(real_dimension<R3Element<Length>> == 3);
  static_assert(!hilbert<R3Element<Length>>);

  static_assert(real_dimension<Displacement<World>> == 3);
  static_assert(hilbert<Displacement<World>>);
  static_assert(real_dimension<AngularVelocity<World>> == 3);
  static_assert(hilbert<AngularVelocity<World>>);
  // Not a field, but a one-dimensional hilbert space.
  static_assert(!homogeneous_field<Trivector<double, World>>);
  static_assert(real_dimension<Trivector<double, World>> == 1);
  static_assert(hilbert<Trivector<double, World>>);
}

}  // namespace geometry
}  // namespace principia
