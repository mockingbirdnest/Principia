
#include <tuple>

#include "numerics/polynomial.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using geometry::Frame;
using geometry::Displacement;
using geometry::Instant;
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

  using P2V = PolynomialInMonomialBasis<Displacement<World>, Time, 2,
                                        HornerEvaluator>;
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2,
                                        HornerEvaluator>;
  using P17 = PolynomialInMonomialBasis<Displacement<World>, Time, 17,
                                        EstrinEvaluator>;

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

  P2V::Coefficients const coefficients_;
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
  P2V p(coefficients_);
  Displacement<World> const d = p.Evaluate(0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

// Check that a polynomial can be for an affine argument.
TEST_F(PolynomialTest, Evaluate3) {
  Instant const t0 = Instant() + 0.3 * Second;
  P2A p(coefficients_, t0);
  Displacement<World> const d = p.Evaluate(t0 + 0.2 * Second);
  Velocity<World> const v = p.EvaluateDerivative(t0 + 0.2 * Second);
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
