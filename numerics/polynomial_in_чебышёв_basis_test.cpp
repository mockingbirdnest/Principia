#include "numerics/polynomial_in_чебышёв_basis.hpp"

#include "astronomy/frames.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {

using ::testing::ElementsAre;
using ::testing::Lt;
using namespace principia::astronomy::_frames;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_numerics_matchers;

class PolynomialInЧебышёвBasisTest : public ::testing::Test {
 protected:
  PolynomialInЧебышёвBasisTest()
      : t_min_(t0_ - 1 * Second),
        t_max_(t0_ + 3 * Second) {}

  Instant const t0_;
  Instant t_min_;
  Instant t_max_;
};

using PolynomialInЧебышёвBasisDeathTest = PolynomialInЧебышёвBasisTest;

TEST_F(PolynomialInЧебышёвBasisDeathTest, EvaluationErrors) {
#ifdef _DEBUG
  using Series = PolynomialInЧебышёвBasis<double, Instant, 0>;
  EXPECT_DEATH({
    Series p({1}, t_min_, t_max_);
    p(t_min_ - 10 * Second);
  }, ">= -1.1");
  EXPECT_DEATH({
    Series p({1}, t_min_, t_max_);
    p(t_max_ + 10 * Second);
  }, "<= 1.1");
#endif
}

TEST_F(PolynomialInЧебышёвBasisTest, T0) {
  PolynomialInЧебышёвBasis<double, Instant, 0> t0({1}, t_min_, t_max_);
  EXPECT_EQ(1, t0(t0_ + 1 * Second));
  EXPECT_EQ(1, t0(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, T1) {
  PolynomialInЧебышёвBasis<double, Instant, 1> t1({0, 1}, t_min_, t_max_);
  EXPECT_EQ(0, t1(t0_ + 1 * Second));
  EXPECT_EQ(1, t1(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, T2) {
  PolynomialInЧебышёвBasis<double, Instant, 2> t2({0, 0, 1}, t_min_, t_max_);
  EXPECT_EQ(1, t2(t0_ + -1 * Second));
  EXPECT_EQ(-1, t2(t0_ + 1 * Second));
  EXPECT_EQ(1, t2(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, T3) {
  PolynomialInЧебышёвBasis<double, Instant, 3> t3({0, 0, 0, 1}, t_min_, t_max_);
  EXPECT_EQ(-1, t3(t0_ + -1 * Second));
  EXPECT_EQ(0, t3(t0_ + 1 * Second));
  EXPECT_EQ(-1, t3(t0_ + 2 * Second));
  EXPECT_EQ(1, t3(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, X5) {
  PolynomialInЧебышёвBasis<double, Instant, 5> x5(
      {0.0, 10.0 / 16.0, 0, 5.0 / 16.0, 0, 1.0 / 16.0},
      t_min_, t_max_);
  EXPECT_EQ(-1, x5(t0_ + -1 * Second));
  EXPECT_EQ(0, x5(t0_ + 1 * Second));
  EXPECT_EQ(1.0 / 1024.0, x5(t0_ + 1.5 * Second));
  EXPECT_EQ(1.0 / 32.0, x5(t0_ + 2 * Second));
  EXPECT_EQ(1, x5(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, X6) {
  PolynomialInЧебышёвBasis<double, Instant, 6> x6(
      {10.0 / 32.0, 0, 15.0 / 32.0, 0, 6.0 / 32.0, 0, 1.0 / 32.0},
      t_min_, t_max_);
  EXPECT_EQ(1, x6(t0_ + -1 * Second));
  EXPECT_EQ(0, x6(t0_ + 1 * Second));
  EXPECT_EQ(1.0 / 4096.0, x6(t0_ + 1.5 * Second));
  EXPECT_EQ(1.0 / 64.0, x6(t0_ + 2 * Second));
  EXPECT_EQ(1, x6(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, T2Dimension) {
  PolynomialInЧебышёвBasis<Length, Instant, 2> t2(
      {0 * Metre, 0 * Metre, 1 * Metre},
      t_min_, t_max_);
  EXPECT_EQ(1 * Metre, t2(t0_ + -1 * Second));
  EXPECT_EQ(-1 * Metre, t2(t0_ + 1 * Second));
  EXPECT_EQ(1 * Metre, t2(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, T2Double) {
  PolynomialInЧебышёвBasis<Length, double, 2> t2(
      {0 * Metre, 0 * Metre, 1 * Metre},
      -1, 2);
  EXPECT_EQ(1 * Metre, t2(-1));
  EXPECT_EQ(-7.0 / 9.0 * Metre, t2(1));
  EXPECT_EQ(1 * Metre, t2(2));
}

TEST_F(PolynomialInЧебышёвBasisTest, Derivative) {
  PolynomialInЧебышёвBasis<double, Instant, 3> series(
      {-2, 3, 5, 6}, t_min_, t_max_);
  EXPECT_EQ(18.5 / Second, series.EvaluateDerivative(t0_ + -1 * Second));
  EXPECT_EQ(-7.5 / Second, series.EvaluateDerivative(t0_ + 1 * Second));
  EXPECT_EQ(38.5 / Second, series.EvaluateDerivative(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisTest, FrobeniusCompanionMatrix) {
  {
    PolynomialInЧебышёвBasis<double, Instant, 3> series(
        {-2, 3, 5, 6}, t_min_, t_max_);
    auto const matrix = series.FrobeniusCompanionMatrix();
    EXPECT_THAT(matrix,
      AlmostEquals(FixedMatrix<double, 3, 3>(
        { 0.0,             1.0,         0.0,
         1.0 / 2.0,       0.0,   1.0 / 2.0,
         1.0 / 6.0, 1.0 / 4.0, -5.0 / 12.0 }),
        0));
    auto const matrix_schur_decomposition =
        RealSchurDecomposition(matrix, 1e-16);
    EXPECT_THAT(matrix_schur_decomposition.real_eigenvalues,
                ElementsAre(AlmostEquals((1.0 - Sqrt(337.0)) / 24.0, 5),
                            AlmostEquals(-0.5, 0),
                            AlmostEquals((1.0 + Sqrt(337.0)) / 24.0, 2)));
  }
  {
    PolynomialInЧебышёвBasis<double, Instant, 1> series(
        {-2, 3}, t_min_, t_max_);
    auto const matrix = series.FrobeniusCompanionMatrix();
    EXPECT_THAT(matrix,
                AlmostEquals(FixedMatrix<double, 1, 1>({1.0 / 3.0}), 0));
    auto const matrix_schur_decomposition =
        RealSchurDecomposition(matrix, 1e-16);
    EXPECT_THAT(matrix_schur_decomposition.real_eigenvalues,
                ElementsAre(AlmostEquals(1.0 / 3.0, 0)));
  }
}

TEST_F(PolynomialInЧебышёвBasisTest, MayHaveRealRoots) {
  // B₀ path.
  PolynomialInЧебышёвBasis<double, Instant, 3> series1(
      {16, 5, 3, 7}, t_min_, t_max_);
  EXPECT_FALSE(series1.MayHaveRealRoots());
  // An error estimate causes the result to be more pessimistic.
  EXPECT_TRUE(series1.MayHaveRealRoots(/*error_estimate*/1.5));
  // We don't know, but it actually doesn't have zeroes.
  PolynomialInЧебышёвBasis<double, Instant, 3> series2(
      {13, 5, 3, 7}, t_min_, t_max_);
  EXPECT_TRUE(series2.MayHaveRealRoots());
  // We don't know, but it actually has zeroes.
  PolynomialInЧебышёвBasis<double, Instant, 3> series3(
      {4, 5, 3, 7}, t_min_, t_max_);
  EXPECT_TRUE(series3.MayHaveRealRoots());
}

TEST_F(PolynomialInЧебышёвBasisTest, RealRoots) {
  PolynomialInЧебышёвBasis<double, Instant, 3> series(
      {-2, 3, 5, 6}, t_min_, t_max_);
  Instant const r1 = t0_ + (13.0 - Sqrt(337.0)) / 12.0 * Second;
  Instant const r2 = t0_;
  Instant const r3 = t0_ + (13.0 + Sqrt(337.0)) / 12.0 * Second;
  EXPECT_THAT(series.RealRoots(1e-16),
              ElementsAre(AlmostEquals(r1, 19),
                          AbsoluteErrorFrom(r2, Lt(2.3e-16 * Second)),
                          AlmostEquals(r3, 1)));
  EXPECT_THAT(Abs(series(r1)), Lt(8.9e-16));
  EXPECT_THAT(Abs(series(r2)), AlmostEquals(0, 0));
  EXPECT_THAT(Abs(series(r3)), Lt(1.8e-15));
}

TEST_F(PolynomialInЧебышёвBasisTest, X6Vector) {
  using V = Vector<Length, ICRS>;
  // {T3, X5, X6}
  V const c0 = V({0.0 * Metre, 0.0 * Metre, 10.0 / 32.0 * Metre});
  V const c1 = V({0.0 * Metre, 10.0 / 16.0 * Metre, 0.0 * Metre});
  V const c2 = V({0.0 * Metre, 0.0 * Metre, 15.0 / 32.0 * Metre});
  V const c3 = V({1.0 * Metre, 5.0 / 16.0 * Metre, 0.0 * Metre});
  V const c4 = V({0.0 * Metre, 0.0 * Metre, 6.0 / 32.0 * Metre});
  V const c5 = V({0.0 * Metre, 1.0 / 16.0 * Metre, 0 * Metre});
  V const c6 = V({0.0 * Metre, 0.0 * Metre, 1.0 / 32.0 * Metre});
  PolynomialInЧебышёвBasis<Vector<Length, ICRS>, Instant, 6> x6(
      {c0, c1, c2, c3, c4, c5, c6},
      t_min_, t_max_);
  EXPECT_EQ(V({-1 * Metre, -1 * Metre, 1 * Metre}),
            x6(t0_ + -1 * Second));
  EXPECT_EQ(V({0 * Metre, 0 * Metre, 0 * Metre}),
            x6(t0_ + 1 * Second));
  EXPECT_EQ(V({-1 * Metre, 1.0 / 32.0 * Metre, 1 / 64.0 * Metre}),
            x6(t0_ + 2 * Second));
  EXPECT_EQ(V({1 * Metre, 1 * Metre, 1 * Metre}),
            x6(t0_ + 3 * Second));
}

TEST_F(PolynomialInЧебышёвBasisDeathTest, SerializationError) {
  using Series1 = PolynomialInЧебышёвBasis<Speed, Instant, 2>;
  using Series2 = PolynomialInЧебышёвBasis<double, Instant, 2>;
  Series1 v({1 * Metre / Second,
             -2 * Metre / Second,
             5 * Metre / Second},
             t_min_, t_max_);
  Series2 d({7, 8, -1}, t_min_, t_max_);

  EXPECT_DEATH({
    serialization::Polynomial message;
    v.WriteToMessage(&message);
    Series2::ReadFromMessage(message);
  }, "has_double");
  EXPECT_DEATH({
    serialization::Polynomial message;
    d.WriteToMessage(&message);
    Series1::ReadFromMessage(message);
  }, "has_quantity");
}

TEST_F(PolynomialInЧебышёвBasisTest, SerializationSuccess) {
  {
    serialization::Polynomial message;
    PolynomialInЧебышёвBasis<Speed, Instant, 2> const v1(
        {1 * Metre / Second, -2 * Metre / Second, 5 * Metre / Second},
        t_min_, t_max_);
    v1.WriteToMessage(&message);
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInЧебышёвBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInЧебышёвBasis::extension);
    EXPECT_EQ(3, extension.coefficient_size());
    EXPECT_FALSE(extension.coefficient(0).has_double_());
    EXPECT_TRUE(extension.coefficient(0).has_quantity());
    EXPECT_EQ(0x7C01, extension.coefficient(0).quantity().dimensions());
    EXPECT_EQ(1.0, extension.coefficient(0).quantity().magnitude());
    EXPECT_TRUE(extension.has_lower_bound());
    EXPECT_TRUE(extension.lower_bound().has_point());
    EXPECT_TRUE(extension.lower_bound().point().has_scalar());
    EXPECT_TRUE(extension.lower_bound().point().scalar().has_dimensions());
    EXPECT_TRUE(extension.lower_bound().point().scalar().has_magnitude());
    EXPECT_EQ(-1.0, extension.lower_bound().point().scalar().magnitude());
    EXPECT_TRUE(extension.has_upper_bound());
    EXPECT_TRUE(extension.upper_bound().has_point());
    EXPECT_TRUE(extension.upper_bound().point().has_scalar());
    EXPECT_TRUE(extension.upper_bound().point().scalar().has_dimensions());
    EXPECT_TRUE(extension.upper_bound().point().scalar().has_magnitude());
    EXPECT_EQ(3.0, extension.upper_bound().point().scalar().magnitude());
    auto const v2 =
        PolynomialInЧебышёвBasis<Speed, Instant, 2>::ReadFromMessage(message);
    EXPECT_EQ(v1, v2);
  }
  {
    serialization::Polynomial message;
    PolynomialInЧебышёвBasis<double, Instant, 2> const d1(
        {-1, 2, 5}, t_min_, t_max_);
    d1.WriteToMessage(&message);
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInЧебышёвBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInЧебышёвBasis::extension);
    EXPECT_EQ(3, extension.coefficient_size());
    EXPECT_TRUE(extension.coefficient(0).has_double_());
    EXPECT_FALSE(extension.coefficient(0).has_quantity());
    EXPECT_EQ(-1.0, extension.coefficient(0).double_());
    EXPECT_TRUE(extension.has_lower_bound());
    EXPECT_TRUE(extension.lower_bound().has_point());
    EXPECT_TRUE(extension.lower_bound().point().has_scalar());
    EXPECT_TRUE(extension.lower_bound().point().scalar().has_dimensions());
    EXPECT_TRUE(extension.lower_bound().point().scalar().has_magnitude());
    EXPECT_EQ(-1.0, extension.lower_bound().point().scalar().magnitude());
    EXPECT_TRUE(extension.has_upper_bound());
    EXPECT_TRUE(extension.upper_bound().has_point());
    EXPECT_TRUE(extension.upper_bound().point().has_scalar());
    EXPECT_TRUE(extension.upper_bound().point().scalar().has_dimensions());
    EXPECT_TRUE(extension.upper_bound().point().scalar().has_magnitude());
    EXPECT_EQ(3.0, extension.upper_bound().point().scalar().magnitude());
    auto const d2 =
        PolynomialInЧебышёвBasis<double, Instant, 2>::ReadFromMessage(message);
    EXPECT_EQ(d1, d2);
  }
}

}  // namespace numerics
}  // namespace principia
