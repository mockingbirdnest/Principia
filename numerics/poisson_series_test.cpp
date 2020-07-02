
#include "numerics/poisson_series.hpp"

#include <memory>

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using geometry::Instant;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Sin;
using quantities::Time;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;

class PoissonSeriesTest : public ::testing::Test {
 protected:
  using Degree1 = PoissonSeries<double, 1, HornerEvaluator>;

  PoissonSeriesTest()
      : ω0_(0 * Radian / Second),
        ω1_(1 * Radian / Second),
        ω2_(2 * Radian / Second),
        ω3_(-3 * Radian / Second) {
    Degree1::Polynomial pa0({0, 0 / Second}, t0_);
    Degree1::Polynomial psa0({100, 200 / Second}, t0_);
    Degree1::Polynomial pca0({1, 2 / Second}, t0_);
    Degree1::Polynomial pb0({3, 4 / Second}, t0_);

    Degree1::Polynomial psa1({5, 6 / Second}, t0_);
    Degree1::Polynomial pca1({7, 8 / Second}, t0_);
    Degree1::Polynomial psb1({9, 10 / Second}, t0_);
    Degree1::Polynomial pcb1({11, 12 / Second}, t0_);

    Degree1::Polynomial psa2({13, 14 / Second}, t0_);
    Degree1::Polynomial pca2({15, 16 / Second}, t0_);

    Degree1::Polynomial psb3({-17, -18 / Second}, t0_);
    Degree1::Polynomial pcb3({19, 20 / Second}, t0_);

    Degree1::Polynomials psca0{/*sin=*/psa0, /*cos=*/pca0};

    Degree1::Polynomials psca1{/*sin=*/psa1, /*cos=*/pca1};
    Degree1::Polynomials pscb1{/*sin=*/psb1, /*cos=*/pcb1};

    Degree1::Polynomials psca2{/*sin=*/psa2, /*cos=*/pca2};

    Degree1::Polynomials pscb3{/*sin=*/psb3, /*cos=*/pcb3};

    pa_ = std::make_unique<Degree1>(
        pa0,
        Degree1::PolynomialsByAngularFrequency{
            {ω0_, psca0}, {ω1_, psca1}, {ω2_, psca2}});
    pb_ = std::make_unique<Degree1>(
        pb0,
        Degree1::PolynomialsByAngularFrequency{{ω1_, pscb1}, {ω3_, pscb3}});
  }

  Instant const t0_;
  AngularFrequency const ω0_;
  AngularFrequency const ω1_;
  AngularFrequency const ω2_;
  AngularFrequency const ω3_;
  std::unique_ptr<Degree1> pa_;
  std::unique_ptr<Degree1> pb_;
};

TEST_F(PoissonSeriesTest, Evaluate) {
  EXPECT_THAT(pa_->Evaluate(t0_ + 1 * Second),
              AlmostEquals(3 + 11 * Sin(1 * Radian) + 15 * Cos(1 * Radian) +
                               27 * Sin(2 * Radian) + 31 * Cos(2 * Radian),
                           0, 1));
  EXPECT_THAT(pb_->Evaluate(t0_ + 1 * Second),
              AlmostEquals(7 + 19 * Sin(1 * Radian) + 23 * Cos(1 * Radian) +
                               35 * Sin(3 * Radian) + 39 * Cos(3 * Radian),
                           32));
}

TEST_F(PoissonSeriesTest, VectorSpace) {
  {
    auto const identity = +*pa_;
    EXPECT_THAT(identity.Evaluate(t0_ + 1 * Second),
                AlmostEquals(pa_->Evaluate(t0_ + 1 * Second), 0));
  }
  {
    auto const negated = -*pb_;
    EXPECT_THAT(negated.Evaluate(t0_ + 1 * Second),
                AlmostEquals(-pb_->Evaluate(t0_ + 1 * Second), 0));
  }
  {
    auto const sum = *pa_ + *pb_;
    EXPECT_THAT(
        sum.Evaluate(t0_ + 1 * Second),
        AlmostEquals(pa_->Evaluate(t0_ + 1 * Second) +
                     pb_->Evaluate(t0_ + 1 * Second), 1));
  }
  {
    auto const difference = *pa_ - *pb_;
    EXPECT_THAT(
        difference.Evaluate(t0_ + 1 * Second),
        AlmostEquals(pa_->Evaluate(t0_ + 1 * Second) -
                     pb_->Evaluate(t0_ + 1 * Second), 0));
  }
  {
    auto const left_product = 3 * *pa_;
    EXPECT_THAT(left_product.Evaluate(t0_ + 1 * Second),
                AlmostEquals(3 * pa_->Evaluate(t0_ + 1 * Second), 1));
  }
  {
    auto const right_product = *pb_ * 4;
    EXPECT_THAT(right_product.Evaluate(t0_ + 1 * Second),
                AlmostEquals(pb_->Evaluate(t0_ + 1 * Second) * 4, 0));
  }
  {
    auto const quotient = *pb_ / 1.5;
    EXPECT_THAT(quotient.Evaluate(t0_ + 1 * Second),
                AlmostEquals(pb_->Evaluate(t0_ + 1 * Second) / 1.5, 0, 32));
  }
}

TEST_F(PoissonSeriesTest, Algebra) {
  auto const product = *pa_ * *pb_;
  EXPECT_THAT(
      product.Evaluate(t0_ + 1 * Second),
      AlmostEquals(pa_->Evaluate(t0_ + 1 * Second) *
                   pb_->Evaluate(t0_ + 1 * Second), 6, 38));
}

TEST_F(PoissonSeriesTest, Primitive) {
  auto const actual_primitive = pb_->Primitive();

  // The primitive was computed using Mathematica.
  auto const expected_primitive = [=](Time const& t){
    auto const a0 = 3;
    auto const a1 = 4 / Second;
    auto const b0 = 9;
    auto const b1 = 10 / Second;
    auto const c0 = 11;
    auto const c1 = 12 / Second;
    auto const d0 = -17;
    auto const d1 = -18 / Second;
    auto const e0 = 19;
    auto const e1 = 20 / Second;
    return a0 * t + (a1 * t * t) / 2 +
           (c1 * Cos(ω1_ * t) * Radian * Radian) / (ω1_ * ω1_) -
           (b0 * Cos(ω1_ * t) * Radian) / ω1_ -
           (b1 * t * Cos(ω1_ * t) * Radian) / ω1_ +
           (e1 * Cos(ω3_ * t) * Radian * Radian) / (ω3_ * ω3_) -
           (d0 * Cos(ω3_ * t) * Radian) / ω3_ -
           (d1 * t * Cos(ω3_ * t) * Radian) / ω3_ +
           (b1 * Sin(ω1_ * t) * Radian * Radian) / (ω1_ * ω1_) +
           (c0 * Sin(ω1_ * t) * Radian) / ω1_ +
           (c1 * t * Sin(ω1_ * t) * Radian) / ω1_ +
           (d1 * Sin(ω3_ * t) * Radian * Radian) / (ω3_ * ω3_) +
           (e0 * Sin(ω3_ * t) * Radian) / ω3_ +
           (e1 * t * Sin(ω3_ * t) * Radian) / ω3_;
  };

  for (int i = -10; i < 10; ++i) {
    EXPECT_THAT(actual_primitive.Evaluate(t0_ + i * Second),
                AlmostEquals(expected_primitive(i * Second), 0, 6));
  }
}

}  // namespace numerics
}  // namespace principia
