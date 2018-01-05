
#include <tuple>

#include "numerics/polynomial.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "serialization/geometry.pb.h"

namespace principia {

using geometry::Frame;
using geometry::Instant;
using geometry::Displacement;
using geometry::Velocity;

namespace numerics {

class PolynomialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  using T =
      PolynomialInMonomialBasis<Displacement<World>, Instant, 3>::Coefficients;

  void Test() {
    T t;
    Displacement<World> d = std::get<0>(t);
    Velocity<World> v = std::get<1>(t);
  }
};

}  // namespace numerics
}  // namespace principia
