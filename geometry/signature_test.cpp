
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
#include "testing_utilities/componentwise.hpp"

namespace principia {

using quantities::Length;
using quantities::si::Metre;
using testing_utilities::Componentwise;
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
      : vector_({1 * Metre, 2 * Metre, 3 * Metre}),
        bivector_({1 * Metre, 2 * Metre, 3 * Metre}),
        trivector_(4 * Metre) {}

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
              Eq(-vector_.coordinates()));
  EXPECT_THAT(NegativeSignature::CentralReflection()(bivector_).coordinates(),
              Eq(bivector_.coordinates()));
  EXPECT_THAT(NegativeSignature::CentralReflection()(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, XYPlaneReflection) {
  NegativeSignature const reflection(
      Sign::Positive(), Sign::Positive(), deduce_sign_from_handedness);
  EXPECT_THAT(reflection(vector_).coordinates(),
              Componentwise(1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(reflection(bivector_).coordinates(),
              Componentwise(-1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(reflection(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, YZPlaneReflection) {
  NegativeSignature const reflection(
      deduce_sign_from_handedness, Sign::Positive(), Sign::Positive());
  EXPECT_THAT(reflection(vector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, 3 * Metre));
  EXPECT_THAT(reflection(bivector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, -3 * Metre));
  EXPECT_THAT(reflection(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, XZPlaneReflection) {
  NegativeSignature rotation(
      Sign::Positive(), deduce_sign_from_handedness, Sign::Positive());
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, XAxisRotation) {
  PositiveSignature const rotation(
      Sign::Positive(), Sign::Negative(), deduce_sign_from_handedness);
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(), Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, YAxisRotation) {
  PositiveSignature const rotation(
      Sign::Negative(), Sign::Positive(), deduce_sign_from_handedness);
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(), Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, ZAxisRotation) {
  PositiveSignature const rotation(
      Sign::Negative(), Sign::Negative(), deduce_sign_from_handedness);
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(-1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(-1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(), Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, Inversion) {
  NegativeSignature const reflection(
      Sign::Negative(), Sign::Positive(), deduce_sign_from_handedness);
  PositiveSignature const rotation(
      Sign::Negative(), Sign::Positive(), deduce_sign_from_handedness);
  EXPECT_THAT(rotation.Inverse()(rotation(vector_)), Eq(vector_));
  EXPECT_THAT(reflection.Inverse()(reflection(vector_)), Eq(vector_));
  EXPECT_THAT(rotation.Inverse()(rotation(bivector_)), Eq(bivector_));
  EXPECT_THAT(reflection.Inverse()(reflection(bivector_)), Eq(bivector_));
  EXPECT_THAT(rotation.Inverse()(rotation(trivector_)), Eq(trivector_));
  EXPECT_THAT(reflection.Inverse()(reflection(trivector_)), Eq(trivector_));
}

TEST_F(SignatureTest, Composition) {
  Signature<R2, L> reflection(
      Sign::Negative(), Sign::Positive(), deduce_sign_from_handedness);
  Signature<R1, R2> rotation(
      Sign::Negative(), Sign::Positive(), deduce_sign_from_handedness);

  EXPECT_THAT((reflection * rotation)(vector_),
              Eq(reflection(rotation(vector_))));
  EXPECT_THAT((reflection * rotation)(bivector_),
              Eq(reflection(rotation(bivector_))));
  EXPECT_THAT((reflection * rotation)(trivector_),
              Eq(reflection(rotation(trivector_))));
}

}  // namespace geometry
}  // namespace principia
