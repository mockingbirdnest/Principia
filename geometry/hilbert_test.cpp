#include "geometry/hilbert.hpp"

#include <type_traits>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

using quantities::Area;
using quantities::Length;
using quantities::Square;
using quantities::si::Metre;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

TEST(HilbertTest, ScalarTypes) {
  using H1 = Hilbert<double, double>;
  static_assert(std::is_same_v<double, H1::InnerProductType>);
  static_assert(std::is_same_v<double, H1::NormType>);

  using H2 = Hilbert<double, Length>;
  static_assert(std::is_same_v<Length, H2::InnerProductType>);
#if 0
  H2::NormType h2;
#endif

  using H3 = Hilbert<Length, Length>;
  static_assert(std::is_same_v<Area, H3::InnerProductType>);
  static_assert(std::is_same_v<Length, H3::NormType>);
}

TEST(HilbertTest, VectorTypes) {
  using H1 = Hilbert<Vector<double, World>, Vector<double, World>>;
  static_assert(std::is_same_v<double, H1::InnerProductType>);
  static_assert(std::is_same_v<double, H1::NormType>);

  using H2 = Hilbert<Bivector<double, World>, Bivector<Length, World>>;
  static_assert(std::is_same_v<Length, H2::InnerProductType>);
#if 0
  H2::NormType h2;
#endif

  using H3 = Hilbert<Trivector<Length, World>, Trivector<Length, World>>;
  static_assert(std::is_same_v<Area, H3::InnerProductType>);
  static_assert(std::is_same_v<Length, H3::NormType>);
}

TEST(HilbertTest, ScalarValues) {
  using H1 = Hilbert<double, double>;
  EXPECT_EQ(6, H1::InnerProduct(2, 3));

  using H2 = Hilbert<double, Length>;
  EXPECT_EQ(6 * Metre, H2::InnerProduct(3, 2 * Metre));

  using H3 = Hilbert<Length, Length>;
  EXPECT_EQ(6 * Metre * Metre, H3::InnerProduct(3 * Metre, 2 * Metre));
}

TEST(HilbertTest, VectorValues) {
  Vector<double, World> v1({1, -2, 3});
  Vector<Length, World> v2({-4 * Metre, 5 * Metre, -6 * Metre});

  using H1 = Hilbert<Vector<double, World>, Vector<double, World>>;
  EXPECT_EQ(14, H1::InnerProduct(v1, v1));

  using H2 = Hilbert<Vector<double, World>, Vector<Length, World>>;
  EXPECT_EQ(-32 * Metre, H2::InnerProduct(v1, v2));

  using H3 = Hilbert<Vector<Length, World>, Vector<Length, World>>;
  EXPECT_EQ(77 * Metre * Metre, H3::InnerProduct(v2, v2));
}

}  // namespace geometry
}  // namespace principia
