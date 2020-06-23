
#include "numerics/poisson_series.hpp"

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
};

TEST_F(PoissonSeriesTest, VectorSpace) {
  Degree1::Polynomial pa0({1, 2 / Second});
  Degree1::Polynomial pb0({3, 4 / Second});

  Degree1::Polynomial psa1({5, 6 / Second});
  Degree1::Polynomial pca1({7, 8 / Second});
  Degree1::Polynomial psb1({9, 10 / Second});
  Degree1::Polynomial pcb1({11, 12 / Second});

  Degree1::Polynomial psa2({13, 14 / Second});
  Degree1::Polynomial pca2({15, 16 / Second});

  Degree1::Polynomial psb3({17, 18 / Second});
  Degree1::Polynomial pcb3({19, 20 / Second});

  Degree1::Polynomials pa1{.sin = psa1, .cos = pca1};
  Degree1::Polynomials pb1{.sin = psb1, .cos = pcb1};

  Degree1::Polynomials pa2{.sin = psa2, .cos = pca2};

  Degree1::Polynomials pb3{.sin = psb3, .cos = pcb3};

  AngularFrequency ω1 = 1 * Radian / Second;
  AngularFrequency ω2 = 2 * Radian / Second;
  AngularFrequency ω3 = 3 * Radian / Second;

  Degree1 pa(pa0, {{ω1, pa1}, {ω2, pa2}});
  Degree1 pb(pb0, {{ω1, pb1}, {ω3, pb3}});

  EXPECT_THAT(pa.Evaluate(1 * Second),
              AlmostEquals(3 + 11 * Sin(1 * Radian) + 15 * Cos(1 * Radian) +
                               27 * Sin(2 * Radian) + 31 * Cos(2 * Radian),
                           0));
}

}  // namespace numerics
}  // namespace principia
