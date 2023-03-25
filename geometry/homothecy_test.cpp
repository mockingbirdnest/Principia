#include "geometry/homothecy.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace geometry {

using ::testing::Eq;
using namespace principia::geometry::_homothecy;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_componentwise;

class HomothecyTest : public testing::Test {
 protected:
  using R1 = Frame<serialization::Frame::TestTag,
                   Inertial,
                   Handedness::Right,
                   serialization::Frame::TEST1>;
  using R2 = Frame<serialization::Frame::TestTag,
                   Inertial,
                   Handedness::Right,
                   serialization::Frame::TEST2>;
  using R3 = Frame<serialization::Frame::TestTag,
                   Inertial,
                   Handedness::Right,
                   serialization::Frame::TEST3>;

  using AmountHomothecy = Homothecy<Amount, R1, R2>;
  using CurrentHomothecy = Homothecy<Current, R2, R3>;

  HomothecyTest()
      : vector_({1 * Metre, -2 * Metre, 3 * Metre}) {}

  Vector<Length, R1> const vector_;
};

TEST_F(HomothecyTest, Vector) {
  AmountHomothecy h(5 * Mole);
  EXPECT_THAT(h(vector_).coordinates(),
              Componentwise(5 * Mole * Metre,
                            -10 * Mole * Metre,
                            15 * Mole * Metre));
}

TEST_F(HomothecyTest, Inverse) {
  AmountHomothecy h(5 * Mole);
  Vector<Length, R2> const vector({1 * Metre, -2 * Metre, 3 * Metre});
  EXPECT_THAT(h.Inverse()(vector).coordinates(),
              Componentwise(AlmostEquals(0.2 * Metre / Mole, 0),
                            AlmostEquals(-0.4 * Metre / Mole, 0),
                            AlmostEquals(0.6 * Metre / Mole, 1)));
}

TEST_F(HomothecyTest, Identity) {
  using H = Homothecy<double, R1, R2>;
  EXPECT_THAT(H::Identity()(vector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, 3 * Metre));
}

TEST_F(HomothecyTest, Forget) {
  // TBD
}

TEST_F(HomothecyTest, Composition) {
  AmountHomothecy h1(5 * Mole);
  CurrentHomothecy h2(3 * Ampere);

  EXPECT_THAT((h2 * h1)(vector_).coordinates(),
              Componentwise(15 * Metre * Mole * Ampere,
                            -30 * Metre * Mole * Ampere,
                            45 * Metre * Mole * Ampere));
}

TEST_F(HomothecyTest, Serialization) {
  serialization::Homothecy message;

  AmountHomothecy homothecy(5 * Mole);
  homothecy.WriteToMessage(&message);
  EXPECT_THAT(AmountHomothecy::ReadFromMessage(message)(vector_),
              Eq(homothecy(vector_)));
}

TEST_F(HomothecyTest, Output) {
  EXPECT_THAT((std::stringstream{} << AmountHomothecy(5 * Mole)).str(),
              Eq("+5.00000000000000000e+00 mol"));
}

}  // namespace geometry
}  // namespace principia
