
#include "geometry/r3_element.hpp"

#include <functional>

#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/bipm.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/algebra.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

#define PRINCIPIA_USE_IACA 0
#if PRINCIPIA_USE_IACA
#include "intel/iacaMarks.h"
#endif

namespace principia {

using quantities::Length;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::astronomy::Parsec;
using quantities::bipm::Knot;
using quantities::constants::SpeedOfLight;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Minute;
using quantities::si::Second;
using quantities::uk::Furlong;
using quantities::uk::Mile;
using quantities::uk::Rod;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;
using ::testing::Eq;

namespace geometry {

class R3ElementTest : public testing::Test {
 protected:
  R3Element<Speed> const null_velocity_ = {0 * Knot, 0 * Knot, 0 * Knot};
  R3Element<Speed> const u_ = {3 * Knot, -42 * Parsec / JulianYear, 0 * Knot};
  R3Element<Speed> const v_ = {-π * SpeedOfLight,
                               -e * Kilo(Metre) / Hour,
                               -1 * Knot};
  R3Element<Speed> const w_ = {2 * Mile / Hour,
                               2 * Furlong / Day,
                               2 * Rod / Minute};
  R3Element<Speed> const a_ = {88 * Mile / Hour,
                               300 * Metre / Second,
                               46 * Knot};
};

using R3ElementDeathTest = R3ElementTest;

#ifdef _DEBUG
TEST_F(R3ElementDeathTest, IndexingOperator) {
  // Sorry about the preprocessing, the regex are not powerful enough to work in
  // both cases.
#ifdef PRINCIPIA_COMPILER_MSVC
  static char const const_method[] =  "\\(int\\) const\\:";
  static char const non_const_method[] =  "\\(int\\)\\:";
#elif PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL
  char const const_method[] =  "\\(int\\) const \\[";
  char const non_const_method[] =  "\\(int\\) \\[";
#endif

  EXPECT_DEATH({
    R3Element<Speed> null_velocity = null_velocity_;
    Speed speed = null_velocity[4];
  }, non_const_method);
  EXPECT_DEATH({
    Speed const& speed = null_velocity_[3];
  }, const_method);
}
#endif

#if PRINCIPIA_USE_IACA
// A convenient skeleton for analysing code with IACA.
TEST_F(R3ElementTest, DISABLED_IACA) {
  auto iaca = [](R3Element<Speed> const& a) {
    IACA_VC64_START;
    auto const result = 2 * Second * a;
    IACA_VC64_END;
    return result;
  };
  CHECK_EQ(iaca(a_), iaca(a_));
}
#endif

TEST_F(R3ElementTest, Dumb3Vector) {
  EXPECT_EQ((e * 42) * v_, e * (42 * v_));
  EXPECT_THAT(303.49234547957647 * Metre / Second, AlmostEquals(a_.Norm(), 0));
  EXPECT_THAT(92107.603764694598 * Metre * Metre / Second / Second,
              AlmostEquals(a_.Norm²(), 0));
  testing_utilities::TestEquality(42 * v_, 43 * v_);
  testing_utilities::TestVectorSpace<R3Element<Speed>, double>(
      null_velocity_, u_, v_, w_, 0.0, 1.0, e, 42.0, 0, 1);
  testing_utilities::TestAlternatingBilinearMap(
      Cross<Speed, Speed>, u_, v_, w_, a_, 42.0, 0, 1);
  EXPECT_EQ(Cross(R3Element<double>(1, 0, 0),
                  R3Element<double>(0, 1, 0)),
            R3Element<double>(0, 0, 1));
  testing_utilities::TestSymmetricPositiveDefiniteBilinearMap(
      Dot<Speed, Speed>, u_, v_, w_, a_, 42.0, 0, 1);
}

TEST_F(R3ElementTest, MixedProduct) {
  testing_utilities::TestBilinearMap(
      std::multiplies<>(), 1 * Second, 1 * JulianYear, u_, v_, 42.0, 0, 1);
  testing_utilities::TestBilinearMap(
      std::multiplies<>(), w_, a_,
      -1 * Day, 1 * Parsec / SpeedOfLight, -π, 0, 1);
  Time const t = -3 * Second;
  EXPECT_EQ(t * u_, u_ * t);
  EXPECT_THAT((u_ * t) / t, AlmostEquals(u_, 1));
}

#ifdef _DEBUG
TEST_F(R3ElementDeathTest, OrthogonalizeError) {
  R3Element<Speed> v1 = {1 * Knot, -2 * Knot, 5 * Knot};
  EXPECT_DEATH({
    v1.OrthogonalizationAgainst(null_velocity_);
  }, "Scalar.*!= norm");
}
#endif

TEST_F(R3ElementTest, OrthogonalizeSuccess) {
  R3Element<Length> const v1 = {1 * Metre, -2 * Metre, 5 * Metre};
  R3Element<Length> const v2 = R3Element<Length>(
      {3 * Metre, 4 * Metre, -1 * Metre}).OrthogonalizationAgainst(v1);
  EXPECT_EQ(0 * Metre * Metre, Dot(v1, v2));
  EXPECT_THAT(v2, AlmostEquals(R3Element<Length>({(10.0 / 3.0) * Metre,
                                                  (10.0 / 3.0) * Metre,
                                                  (2.0 / 3.0) * Metre}), 1));
}

TEST_F(R3ElementTest, Normalize) {
  R3Element<Length> const v = {1 * Metre, -2 * Metre, 5 * Metre};
  EXPECT_THAT(Normalize(v),
              AlmostEquals(R3Element<double>({1.0 / sqrt(30.0),
                                              -2.0 / sqrt(30.0),
                                              5.0 / sqrt(30.0)}), 0));
}

TEST_F(R3ElementDeathTest, SerializationError) {
  R3Element<Speed> v = {1 * Metre / Second,
                        -2 * Metre / Second,
                        5 * Metre / Second};
  R3Element<double> d = {7, 8, -1};

  EXPECT_DEATH({
    serialization::R3Element message;
    v.WriteToMessage(&message);
    d = R3Element<double>::ReadFromMessage(message);
  }, "has_double");
  EXPECT_DEATH({
    serialization::R3Element message;
    d.WriteToMessage(&message);
    v = R3Element<Speed>::ReadFromMessage(message);
  }, "has_quantity");
}

TEST_F(R3ElementTest, SerializationSuccess) {
  serialization::R3Element message;

  R3Element<Speed> const v1 = {1 * Metre / Second,
                               -2 * Metre / Second,
                               5 * Metre / Second};
  v1.WriteToMessage(&message);
  EXPECT_FALSE(message.x().has_double_());
  EXPECT_TRUE(message.x().has_quantity());
  EXPECT_EQ(0x7C01, message.x().quantity().dimensions());
  EXPECT_EQ(1.0, message.x().quantity().magnitude());
  EXPECT_FALSE(message.y().has_double_());
  EXPECT_TRUE(message.y().has_quantity());
  EXPECT_EQ(0x7C01, message.y().quantity().dimensions());
  EXPECT_EQ(-2.0, message.y().quantity().magnitude());
  EXPECT_FALSE(message.z().has_double_());
  EXPECT_TRUE(message.z().has_quantity());
  EXPECT_EQ(0x7C01, message.z().quantity().dimensions());
  EXPECT_EQ(5.0, message.z().quantity().magnitude());
  R3Element<Speed> const v2 = R3Element<Speed>::ReadFromMessage(message);
  EXPECT_EQ(v1, v2);

  R3Element<double> const d1 = {-1, 2, 5};
  d1.WriteToMessage(&message);
  EXPECT_TRUE(message.x().has_double_());
  EXPECT_FALSE(message.x().has_quantity());
  EXPECT_EQ(-1.0, message.x().double_());
  EXPECT_TRUE(message.y().has_double_());
  EXPECT_FALSE(message.y().has_quantity());
  EXPECT_EQ(2.0, message.y().double_());
  EXPECT_TRUE(message.z().has_double_());
  EXPECT_FALSE(message.z().has_quantity());
  EXPECT_EQ(5.0, message.z().double_());
  R3Element<double> const d2 = R3Element<double>::ReadFromMessage(message);
  EXPECT_EQ(d1, d2);
}

TEST_F(R3ElementTest, SphericalCoordinates) {
  R3Element<Length> x{1 * Metre, 0 * Metre, 0 * Metre};
  R3Element<Length> y{0 * Metre, 1 * Metre, 0 * Metre};
  R3Element<Length> z{0 * Metre, 0 * Metre, 1 * Metre};
  R3Element<Length> v{1 * Metre, -1 * Metre, -Sqrt(6) * Metre};

  EXPECT_EQ(0 * Degree, x.ToSpherical().latitude);
  EXPECT_EQ(0 * Degree, x.ToSpherical().longitude);
  EXPECT_EQ(1 * Metre, x.ToSpherical().radius);

  EXPECT_EQ(0 * Degree, y.ToSpherical().latitude);
  EXPECT_EQ(90 * Degree, y.ToSpherical().longitude);
  EXPECT_EQ(1 * Metre, y.ToSpherical().radius);

  EXPECT_EQ(90 * Degree, z.ToSpherical().latitude);
  // The x and y coordinates are both +0, so the angle is +0, see Kahan (1986).
  EXPECT_EQ(0 * Degree, z.ToSpherical().longitude);
  EXPECT_EQ(1 * Metre, z.ToSpherical().radius);

  EXPECT_EQ(-60 * Degree, v.ToSpherical().latitude);
  EXPECT_EQ(-45 * Degree, v.ToSpherical().longitude);
  EXPECT_THAT(v.ToSpherical().radius, AlmostEquals(Sqrt(8) * Metre, 1));

  EXPECT_EQ(x, RadiusLatitudeLongitude(1 * Metre, 0 * Degree, 0 * Degree)
                   .ToCartesian());
  EXPECT_THAT(
      RadiusLatitudeLongitude(1 * Metre, 0 * Degree, 90 * Degree).ToCartesian(),
      Componentwise(VanishesBefore(1 * Metre, 0),
                    Eq(1 * Metre),
                    Eq(0 * Metre)));
  EXPECT_THAT(RadiusLatitudeLongitude(1 * Metre,
                                      90 * Degree,
                                      42 * Degree).ToCartesian(),
              Componentwise(VanishesBefore(1 * Metre, 0),
                            VanishesBefore(1 * Metre, 0),
                            Eq(1 * Metre)));
}

}  // namespace geometry
}  // namespace principia

#undef PRINCIPIA_USE_IACA
