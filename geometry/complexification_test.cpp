#include "geometry/complexification.hpp"

#include <sstream>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

using ::testing::Eq;
using namespace principia::geometry::_complexification;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_space;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

class ComplexificationTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  Displacement<World> const v1_{{1 * Metre, 4 * Metre, 9 * Metre}};
  Displacement<World> const v2_{{2 * Metre, 3 * Metre, 5 * Metre}};
  Displacement<World> const v3_{{-1 * Metre, 3 * Metre, -9 * Metre}};
  Displacement<World> const v4_{{0 * Metre, 10 * Metre, 5 * Metre}};
  Complexification<double> const i_{0, 1};
};

TEST_F(ComplexificationTest, Addition) {
  {
    Complexification<Displacement<World>> const negation =
        -Complexification(v1_, v2_);
    EXPECT_THAT(negation.real_part(),
                Eq(Displacement<World>({-1 * Metre, -4 * Metre, -9 * Metre})));
    EXPECT_THAT(negation.imaginary_part(),
                Eq(Displacement<World>({-2 * Metre, -3 * Metre, -5 * Metre})));
  }
  {
    Complexification<Displacement<World>> const sum =
        Complexification(v1_, v2_) + Complexification(v3_, v4_);
    EXPECT_THAT(sum.real_part(),
                Eq(Displacement<World>({0 * Metre, 7 * Metre, 0 * Metre})));
    EXPECT_THAT(sum.imaginary_part(),
                Eq(Displacement<World>({2 * Metre, 13 * Metre, 10 * Metre})));
  }
  {
    Complexification<Displacement<World>> const difference =
        Complexification(v1_, v2_) - Complexification(v3_, v4_);
    EXPECT_THAT(difference.real_part(),
                Eq(Displacement<World>({2 * Metre, 1 * Metre, 18 * Metre})));
    EXPECT_THAT(difference.imaginary_part(),
                Eq(Displacement<World>({2 * Metre, -7 * Metre, 0 * Metre})));
  }
}

TEST_F(ComplexificationTest, Multiplication) {
  auto const& v = v1_;
  // Complex * real.
  Complexification<Displacement<World>> iv = i_ * v;
  EXPECT_THAT(iv.real_part(), Eq(Displacement<World>{}));
  EXPECT_THAT(iv.imaginary_part(), Eq(v));
  // Real * complex.
  EXPECT_THAT((2 * iv).real_part(), Eq(Displacement<World>{}));
  EXPECT_THAT((2 * iv).imaginary_part(), Eq(2 * v));
  EXPECT_THAT((iv * 2).real_part(), Eq(Displacement<World>{}));
  EXPECT_THAT((iv * 2).imaginary_part(), Eq(2 * v));
  // Complex * complex.
  Complexification<Displacement<World>> i²v = i_ * iv;
  EXPECT_THAT(i²v.real_part(), Eq(-v));
  EXPECT_THAT(i²v.imaginary_part(), Eq(Displacement<World>{}));

  Complexification<double> i² = i_ * i_;
  EXPECT_THAT(i².real_part(), Eq(-1));
  EXPECT_THAT(i².imaginary_part(), Eq(0));
}

TEST_F(ComplexificationTest, Division) {
  // Complex / real, the easy one.
  EXPECT_THAT((v1_ + i_ * v2_) / 2, Eq(0.5 * (v1_ + i_ * v2_)));
  // Real / complex.
  EXPECT_THAT(1 / i_, Eq(-i_));
  Complexification<Length> const z =
      v1_.coordinates().x + i_ * v2_.coordinates().x;
  // Complex / complex.
  Complexification<Vector<double, World>> quotient = (v1_ + i_ * v2_) / z;
  EXPECT_THAT(quotient, Eq(v1_ / z + i_ * v2_ / z));

  EXPECT_THAT(quotient.real_part().coordinates().x, Eq(1));
  EXPECT_THAT(quotient.imaginary_part().coordinates().x, Eq(0));
}

TEST_F(ComplexificationTest, Norm²) {
  auto const v = v1_ + i_ * v2_;
  auto const vx = v1_.coordinates().x + i_ * v2_.coordinates().x;
  auto const vy = v1_.coordinates().y + i_ * v2_.coordinates().y;
  auto const vz = v1_.coordinates().z + i_ * v2_.coordinates().z;

  EXPECT_THAT(v.Norm²(),
              Eq((vx * vx.Conjugate() +
                  vy * vy.Conjugate() +
                  vz * vz.Conjugate()).real_part()));
}

TEST_F(ComplexificationTest, Logging) {
  EXPECT_THAT((std::stringstream() << (2 * i_ - 1)).str(), Eq("-1 + 2 i"));
}

}  // namespace geometry
}  // namespace principia
