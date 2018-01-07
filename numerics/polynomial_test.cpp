
#include <tuple>

#include "numerics/polynomial.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using geometry::Frame;
using geometry::Displacement;
using geometry::Vector;
using geometry::Velocity;
using quantities::Acceleration;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;

namespace numerics {

class PolynomialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  using P2 = PolynomialInMonomialBasis<Displacement<World>, Time, 2>;
  using P17 = PolynomialInMonomialBasis<Displacement<World>, Time, 17>;

  PolynomialTest()
      : coefficients_({
            Displacement<World>({0 * Metre,
                                 0 * Metre,
                                 1 * Metre}),
            Velocity<World>({0 * Metre / Second,
                             1 * Metre / Second,
                             0 * Metre / Second}),
            Vector<Acceleration, World>({1 * Metre / Second / Second,
                                         0 * Metre / Second / Second,
                                         0 * Metre / Second / Second})}) {}

  P2::Coefficients const coefficients_;
};

// Check that coefficients can be accessed and have the right type.
TEST_F(PolynomialTest, Coefficients) {
  Displacement<World> const d = std::get<0>(coefficients_);
  Velocity<World> const v = std::get<1>(coefficients_);
  EXPECT_EQ(1 * Metre, d.coordinates().z);
  EXPECT_EQ(1 * Metre / Second, v.coordinates().y);
}

// Check that a polynomial can be constructed and evaluated.
TEST_F(PolynomialTest, Evaluate2) {
  P2 p(coefficients_);
  Displacement<World> const d = p.Evaluate(0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

// Check that a polynomial of high order may be declared even if the quantities
// of the coefficients would not be serializable.
TEST_F(PolynomialTest, Evaluate17) {
  P17::Coefficients const coefficients;
  P17 p(coefficients);
  Displacement<World> const d = p.Evaluate(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}), 0));
  // The following doesn't compile with: "Invalid time exponent".
  //   serialization::Quantity message;
  //   std::get<17>(coefficients).coordinates().x.WriteToMessage(&message);
}

}  // namespace numerics
}  // namespace principia
