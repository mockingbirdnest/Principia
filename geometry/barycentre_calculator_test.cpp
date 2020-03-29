
#include "geometry/barycentre_calculator.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using geometry::Frame;
using quantities::Entropy;
using quantities::KinematicViscosity;
using testing_utilities::AlmostEquals;
namespace si = quantities::si;

namespace geometry {

class BarycentreCalculatorTest : public testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  BarycentreCalculatorTest()
      : b1_({1 * si::Unit<Entropy>,
             -2 * si::Unit<Entropy>,
             3 * si::Unit<Entropy>}),
        b2_({-9 * si::Unit<Entropy>,
             8 * si::Unit<Entropy>,
             7 * si::Unit<Entropy>}),
        k1_(4 * si::Unit<KinematicViscosity>),
        k2_(5 * si::Unit<KinematicViscosity>) {}

  Bivector<Entropy, World> b1_;
  Bivector<Entropy, World> b2_;
  KinematicViscosity k1_;
  KinematicViscosity k2_;
};

using BarycentreCalculatorDeathTest = BarycentreCalculatorTest;

TEST_F(BarycentreCalculatorDeathTest, Error) {
  using Calculator = BarycentreCalculator<Bivector<Entropy, World>, double>;
  EXPECT_DEATH({
    Calculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
}

TEST_F(BarycentreCalculatorTest, Bivector) {
  BarycentreCalculator<Bivector<Entropy, World>, KinematicViscosity>
      barycentre_calculator;
  barycentre_calculator.Add(b1_, k1_);
  barycentre_calculator.Add(b2_, k2_);
  EXPECT_THAT(barycentre_calculator.Get(),
              AlmostEquals(
                  Bivector<Entropy, World>({(-41.0 / 9.0) * si::Unit<Entropy>,
                                            (32.0 / 9.0) * si::Unit<Entropy>,
                                            (47.0 / 9.0) * si::Unit<Entropy>}),
                  0));
}

TEST_F(BarycentreCalculatorTest, Scalar) {
  BarycentreCalculator<KinematicViscosity, double> barycentre_calculator;
  barycentre_calculator.Add(k1_, -3);
  barycentre_calculator.Add(k2_, 7);
  EXPECT_THAT(barycentre_calculator.Get(),
              AlmostEquals((23.0 / 4.0) * si::Unit<KinematicViscosity>, 0));
}

}  // namespace geometry
}  // namespace principia
