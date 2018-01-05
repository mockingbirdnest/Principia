
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

namespace principia {

using geometry::Frame;
using geometry::Instant;
using geometry::Displacement;
using geometry::Vector;
using geometry::Velocity;
using quantities::Acceleration;
using quantities::si::Metre;
using quantities::si::Second;

namespace numerics {

class PolynomialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  using P2 = PolynomialInMonomialBasis<Displacement<World>, Instant, 2>;

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
  Instant const t0_;
};

// Check that coefficients can be accessed and have the right type.
TEST_F(PolynomialTest, Coefficients) {
  Displacement<World> const d = std::get<0>(coefficients_);
  Velocity<World> const v = std::get<1>(coefficients_);
  EXPECT_EQ(1 * Metre, d.coordinates().z);
  EXPECT_EQ(1 * Metre / Second, v.coordinates().y);
}

// Check that a polynomial can be constructed and evaluated.
TEST_F(PolynomialTest, Evaluate) {
  P2 p(coefficients_, t0_, t0_ + 1 * Second);
  Displacement<World> const d = p.Evaluate(t0_ + 0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(t0_ + 0.5 * Second);
}

}  // namespace numerics
}  // namespace principia
