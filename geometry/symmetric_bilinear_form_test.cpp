
#include "geometry/symmetric_bilinear_form.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Length;
using quantities::Square;
using quantities::si::Metre;
using ::testing::Eq;

class SymmetricBilinearFormTest : public ::testing::Test {
 protected:
  using World =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;

  SymmetricBilinearFormTest()
      : v1_(Vector<Length, World>({1 * Metre, 2 * Metre, -4 * Metre})),
        v2_(Vector<Length, World>({2 * Metre, -3 * Metre, 5 * Metre})) {}

  Vector<Length, World> v1_;
  Vector<Length, World> v2_;
};

TEST_F(SymmetricBilinearFormTest, UnaryOperators) {
  auto const f1 =
      SymmetricProduct<Length, Length, World>(v1_, v2_) +
      SymmetricBilinearForm<Square<Length>, World>::InnerProductForm();
  auto const f2 =
      SymmetricProduct<Length, Length, World>(-v1_, v2_) -
      SymmetricBilinearForm<Square<Length>, World>::InnerProductForm();
  EXPECT_THAT(f1, Eq(+f1));
  EXPECT_THAT(f2, Eq(-f1));
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
