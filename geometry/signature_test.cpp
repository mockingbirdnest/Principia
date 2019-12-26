
#include "geometry/signature.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/identity.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/componentwise.hpp"

namespace principia {

using quantities::Length;
using quantities::Pow;
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
        trivector_(4 * Metre),
        form_(SymmetricBilinearForm<Length, R1>(
            R3x3Matrix<Length>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre},
                               {2.0 * Metre, -5.0 * Metre, 6.0 * Metre},
                               {3.0 * Metre, 6.0 * Metre, 4.0 * Metre}))) {}

  Vector<Length, R1> vector_;
  Bivector<Length, R1> bivector_;
  Trivector<Length, R1> trivector_;
  SymmetricBilinearForm<Length, R1> const form_;
};

using SignatureDeathTest = SignatureTest;

TEST_F(SignatureTest, Forget) {
  std::array<PositiveSignature, 4> all_positive_signatures{
      {PositiveSignature::Identity(),
       PositiveSignature(Sign::Positive(),
                         Sign::Negative(),
                         DeduceSignPreservingOrientation{}),
       PositiveSignature(Sign::Negative(),
                         Sign::Positive(),
                         DeduceSignPreservingOrientation{}),
       PositiveSignature(Sign::Negative(),
                         Sign::Negative(),
                         DeduceSignPreservingOrientation{})}};
  std::array<NegativeSignature, 4> all_negative_signatures{
      {NegativeSignature::CentralInversion(),
       NegativeSignature(Sign::Negative(),
                         Sign::Positive(),
                         DeduceSignReversingOrientation{}),
       NegativeSignature(Sign::Positive(),
                         Sign::Negative(),
                         DeduceSignReversingOrientation{}),
       NegativeSignature(Sign::Positive(),
                         Sign::Positive(),
                         DeduceSignReversingOrientation{})}};
  auto const test_forget_for = [this](auto const& signatures) {
    for (auto const& signature : signatures) {
      EXPECT_THAT(signature.Forget()(vector_), Eq(signature(vector_)))
          << signature;
      EXPECT_THAT(signature.Forget()(bivector_), Eq(signature(bivector_)))
          << signature;
      EXPECT_THAT(signature.Forget()(trivector_), Eq(signature(trivector_)))
          << signature;
    }
  };
  test_forget_for(all_positive_signatures);
  test_forget_for(all_negative_signatures);
}

TEST_F(SignatureTest, Identity) {
  EXPECT_THAT(PositiveSignature::Identity()(vector_).coordinates(),
              Eq(vector_.coordinates()));
  EXPECT_THAT(PositiveSignature::Identity()(bivector_).coordinates(),
              Eq(bivector_.coordinates()));
  EXPECT_THAT(PositiveSignature::Identity()(trivector_).coordinates(),
              Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, CentralInversion) {
  EXPECT_THAT(NegativeSignature::CentralInversion()(vector_).coordinates(),
              Eq(-vector_.coordinates()));
  EXPECT_THAT(NegativeSignature::CentralInversion()(bivector_).coordinates(),
              Eq(bivector_.coordinates()));
  EXPECT_THAT(NegativeSignature::CentralInversion()(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, XYPlaneReflection) {
  NegativeSignature const reflection(
      Sign::Positive(), Sign::Positive(), DeduceSignReversingOrientation{});
  EXPECT_THAT(reflection(vector_).coordinates(),
              Componentwise(1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(reflection(bivector_).coordinates(),
              Componentwise(-1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(reflection(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, YZPlaneReflection) {
  NegativeSignature const reflection(
      DeduceSignReversingOrientation{}, Sign::Positive(), Sign::Positive());
  EXPECT_THAT(reflection(vector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, 3 * Metre));
  EXPECT_THAT(reflection(bivector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, -3 * Metre));
  EXPECT_THAT(reflection(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, XZPlaneReflection) {
  NegativeSignature rotation(
      Sign::Positive(), DeduceSignReversingOrientation{}, Sign::Positive());
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(),
              Eq(-trivector_.coordinates()));
}

TEST_F(SignatureTest, XAxisRotation) {
  PositiveSignature const rotation(
      Sign::Positive(), Sign::Negative(), DeduceSignPreservingOrientation{});
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(1 * Metre, -2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(), Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, YAxisRotation) {
  PositiveSignature const rotation(
      Sign::Negative(), Sign::Positive(), DeduceSignPreservingOrientation{});
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(-1 * Metre, 2 * Metre, -3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(), Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, ZAxisRotation) {
  PositiveSignature const rotation(
      Sign::Negative(), Sign::Negative(), DeduceSignPreservingOrientation{});
  EXPECT_THAT(rotation(vector_).coordinates(),
              Componentwise(-1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(rotation(bivector_).coordinates(),
              Componentwise(-1 * Metre, -2 * Metre, 3 * Metre));
  EXPECT_THAT(rotation(trivector_).coordinates(), Eq(trivector_.coordinates()));
}

TEST_F(SignatureTest, Inversion) {
  NegativeSignature const reflection(
      Sign::Negative(), Sign::Positive(), DeduceSignReversingOrientation{});
  PositiveSignature const rotation(
      Sign::Negative(), Sign::Positive(), DeduceSignPreservingOrientation{});
  EXPECT_THAT(rotation.Inverse()(rotation(vector_)), Eq(vector_));
  EXPECT_THAT(reflection.Inverse()(reflection(vector_)), Eq(vector_));
  EXPECT_THAT(rotation.Inverse()(rotation(bivector_)), Eq(bivector_));
  EXPECT_THAT(reflection.Inverse()(reflection(bivector_)), Eq(bivector_));
  EXPECT_THAT(rotation.Inverse()(rotation(trivector_)), Eq(trivector_));
  EXPECT_THAT(reflection.Inverse()(reflection(trivector_)), Eq(trivector_));
}

TEST_F(SignatureTest, Composition) {
  Signature<R2, L> reflection(
      Sign::Negative(), Sign::Positive(), DeduceSignReversingOrientation{});
  Signature<R1, R2> rotation(
      Sign::Negative(), Sign::Positive(), DeduceSignPreservingOrientation{});

  EXPECT_THAT((reflection * rotation)(vector_),
              Eq(reflection(rotation(vector_))));
  EXPECT_THAT((reflection * rotation)(bivector_),
              Eq(reflection(rotation(bivector_))));
  EXPECT_THAT((reflection * rotation)(trivector_),
              Eq(reflection(rotation(trivector_))));
}

TEST_F(SignatureTest, AppliedToSymmetricBilinearForm) {
  Signature<R1, L> reflection(
      Sign::Negative(), Sign::Positive(), DeduceSignReversingOrientation{});
  Vector<Length, L> vector({ 1 * Metre, 2 * Metre, 3 * Metre });

  EXPECT_THAT(form_(vector_, vector_), 115 * Pow<3>(Metre));
  EXPECT_THAT(reflection(form_)(vector, vector), 63 * Pow<3>(Metre));
}

TEST_F(SignatureTest, Serialization) {
  serialization::Signature message;

  PositiveSignature signature(
      Sign::Positive(), Sign::Negative(), DeduceSignPreservingOrientation{});
  signature.WriteToMessage(&message);
  EXPECT_THAT(PositiveSignature::ReadFromMessage(message)(vector_),
              Eq(signature(vector_)));
}

TEST_F(SignatureTest, Output) {
  EXPECT_THAT((std::stringstream{}
               << PositiveSignature(Sign::Positive(),
                                    Sign::Negative(),
                                    DeduceSignPreservingOrientation{}))
                  .str(),
              Eq("{+, -, -}"));
}

}  // namespace geometry
}  // namespace principia
