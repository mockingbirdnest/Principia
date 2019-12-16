
#include "geometry/signature.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/identity.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using quantities::Length;
using quantities::si::Metre;
using ::testing::Eq;

namespace geometry {

class SignatureTest : public testing::Test {
 protected:
  using R1 = Frame<serialization::Frame::TestTag,
                  Inertial,
                  Handedness::Right,
                   serialization::Frame::TEST1>;
  using R2 = Frame<serialization::Frame::TestTag,
                   Inertial,
                   Handedness::Right,
                   serialization::Frame::TEST2>;
  using L = Frame<enum class LTag, Inertial, Handedness::Left>;

  using PositiveSignature = Signature<R1, R2>;
  using NegativeSignature = Signature<R1, L>;

  SignatureTest()
      : vector_({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}),
        bivector_({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}),
        trivector_(4.0 * Metre) {}

  Vector<Length, R1> vector_;
  Bivector<Length, R1> bivector_;
  Trivector<Length, R1> trivector_;
};

using SignatureDeathTest = SignatureTest;

TEST_F(SignatureTest, Identity) {
  EXPECT_THAT(PositiveSignature::Identity()(vector_).coordinates(),
              Eq(vector_.coordinates()));
  EXPECT_THAT(PositiveSignature::Identity()(bivector_).coordinates(),
              Eq(bivector_.coordinates()));
  EXPECT_THAT(PositiveSignature::Identity()(trivector_).coordinates(),
              Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, CentralReflection) {
  EXPECT_THAT(NegativeSignature::CentralReflection()(vector_).coordinates(),
              Eq(vector_.coordinates()));
  EXPECT_THAT(NegativeSignature::CentralReflection()(bivector_).coordinates(),
              Eq(bivector_.coordinates()));
  EXPECT_THAT(NegativeSignature::CentralReflection()(trivector_).coordinates(),
              Eq(trivector_.coordinates()));
}

}  // namespace geometry
}  // namespace principia
