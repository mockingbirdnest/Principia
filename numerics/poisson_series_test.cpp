
#include "numerics/poisson_series.hpp"

#include <memory>

#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Sin;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;

class PoissonSeriesTest : public ::testing::Test {
 protected:
  using Degree1 = PoissonSeries<double, 1, HornerEvaluator>;

  PoissonSeriesTest() {
    Degree1::Polynomial pa0({0, 0 / Second});
    Degree1::Polynomial psa0({100, 200 / Second});
    Degree1::Polynomial pca0({1, 2 / Second});
    Degree1::Polynomial pb0({3, 4 / Second});

    Degree1::Polynomial psa1({5, 6 / Second});
    Degree1::Polynomial pca1({7, 8 / Second});
    Degree1::Polynomial psb1({9, 10 / Second});
    Degree1::Polynomial pcb1({11, 12 / Second});

    Degree1::Polynomial psa2({13, 14 / Second});
    Degree1::Polynomial pca2({15, 16 / Second});

    Degree1::Polynomial psb3({-17, -18 / Second});
    Degree1::Polynomial pcb3({19, 20 / Second});

    Degree1::Polynomials psca0{/*sin=*/psa0, /*cos=*/pca0};

    Degree1::Polynomials psca1{/*sin=*/psa1, /*cos=*/pca1};
    Degree1::Polynomials pscb1{/*sin=*/psb1, /*cos=*/pcb1};

    Degree1::Polynomials psca2{/*sin=*/psa2, /*cos=*/pca2};

    Degree1::Polynomials pscb3{/*sin=*/psb3, /*cos=*/pcb3};

    AngularFrequency ω0 = 0 * Radian / Second;
    AngularFrequency ω1 = 1 * Radian / Second;
    AngularFrequency ω2 = 2 * Radian / Second;
    AngularFrequency ω3 = -3 * Radian / Second;

    pa_ = std::make_unique<Degree1>(pa0,
                                    Degree1::PolynomialsByAngularFrequency{
                                        {ω0, psca0}, {ω1, psca1}, {ω2, psca2}});
    pb_ = std::make_unique<Degree1>(pb0,
                                    Degree1::PolynomialsByAngularFrequency{
                                        {ω1, pscb1}, {ω3, pscb3}});
  }

  std::unique_ptr<Degree1> pa_;
  std::unique_ptr<Degree1> pb_;
};

TEST_F(PoissonSeriesTest, Evaluate) {
  EXPECT_THAT(pa_->Evaluate(1 * Second),
              AlmostEquals(3 + 11 * Sin(1 * Radian) + 15 * Cos(1 * Radian) +
                               27 * Sin(2 * Radian) + 31 * Cos(2 * Radian),
                           0, 1));
  EXPECT_THAT(pb_->Evaluate(1 * Second),
              AlmostEquals(7 + 19 * Sin(1 * Radian) + 23 * Cos(1 * Radian) +
                               35 * Sin(3 * Radian) + 39 * Cos(3 * Radian),
                           32));
}

TEST_F(PoissonSeriesTest, VectorSpace) {
  {
    auto const identity = +*pa_;
    EXPECT_THAT(identity.Evaluate(1 * Second),
                AlmostEquals(pa_->Evaluate(1 * Second), 0));
  }
  {
    auto const negated = -*pb_;
    EXPECT_THAT(negated.Evaluate(1 * Second),
                AlmostEquals(-pb_->Evaluate(1 * Second), 0));
  }
  {
    auto const sum = *pa_ + *pb_;
    EXPECT_THAT(
        sum.Evaluate(1 * Second),
        AlmostEquals(pa_->Evaluate(1 * Second) + pb_->Evaluate(1 * Second), 1));
  }
  {
    auto const difference = *pa_ - *pb_;
    EXPECT_THAT(
        difference.Evaluate(1 * Second),
        AlmostEquals(pa_->Evaluate(1 * Second) - pb_->Evaluate(1 * Second), 0));
  }
  {
    auto const left_product = 3 * *pa_;
    EXPECT_THAT(left_product.Evaluate(1 * Second),
                AlmostEquals(3 * pa_->Evaluate(1 * Second), 1));
  }
  {
    auto const right_product = *pb_ * 4;
    EXPECT_THAT(right_product.Evaluate(1 * Second),
                AlmostEquals(pb_->Evaluate(1 * Second) * 4, 0));
  }
  {
    auto const quotient = *pb_ / 1.5;
    EXPECT_THAT(quotient.Evaluate(1 * Second),
                AlmostEquals(pb_->Evaluate(1 * Second) / 1.5, 0, 32));
  }
}

}  // namespace numerics
}  // namespace principia
