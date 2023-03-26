#include "geometry/conformal_map.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/rotation.hpp"
#include "geometry/signature.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using ::testing::Eq;
using namespace principia::geometry::_conformal_map;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

class ConformalMapTest : public testing::Test {
 protected:
  using DirectWorld = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST1>;
  using MirrorWorld = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Left,
                      serialization::Frame::TEST2>;
  using DirectConf = ConformalMap<Amount, DirectWorld, DirectWorld>;
  using MirrorConf = ConformalMap<Amount, MirrorWorld, DirectWorld>;
  using DirectOrth = OrthogonalMap<DirectWorld, DirectWorld>;
  using MirrorOrth = OrthogonalMap<MirrorWorld, DirectWorld>;
  using Rot = Rotation<DirectWorld, DirectWorld>;
  using DirectSig = Signature<DirectWorld, DirectWorld>;
  using MirrorSig = Signature<MirrorWorld, DirectWorld>;

  ConformalMapTest()
      : direct_vector_(Vector<Length, DirectWorld>(
            {1.0 * Metre, 2.0 * Metre, 3.0 * Metre})),
        mirror_vector_(Vector<Length, MirrorWorld>(
            {1.0 * Metre, 2.0 * Metre, 3.0 * Metre})),
        conformal_a_(
            MirrorConf(5 * Mole,
                       MirrorOrth(Rot(120 * Degree,
                                      Bivector<double, DirectWorld>({1, 1, 1}))
                                      .quaternion()))),
        conformal_b_(
            DirectConf(5 * Mole,
                       DirectOrth(Rot(90 * Degree,
                                      Bivector<double, DirectWorld>({1, 0, 0}))
                                      .quaternion()))),
        conformal_c_(
            MirrorConf(5 * Mole,
                       MirrorOrth(Rot(90 * Degree,
                                      Bivector<double, DirectWorld>({1, 0, 0}))
                                      .quaternion()))) {}

  Vector<Length, DirectWorld> direct_vector_;
  Vector<Length, MirrorWorld> mirror_vector_;
  MirrorConf conformal_a_;
  DirectConf conformal_b_;
  MirrorConf conformal_c_;
};

using ConformalMapDeathTest = ConformalMapTest;

TEST_F(ConformalMapTest, AppliedToVector) {
  EXPECT_THAT(conformal_a_(mirror_vector_),
              AlmostEquals(Vector<Product<Length, Amount>, DirectWorld>(
                               {15.0 * Metre * Mole,
                                5.0 * Metre * Mole,
                                10.0 * Metre * Mole}),
                  4));
  EXPECT_THAT(conformal_b_(direct_vector_),
              AlmostEquals(Vector<Product<Length, Amount>, DirectWorld>(
                               {5.0 * Metre * Mole,
                                -15.0 * Metre * Mole,
                                10.0 * Metre * Mole}),
                  1,
                  2));
}

TEST_F(ConformalMapTest, Inverse) {
  EXPECT_THAT(conformal_a_.Inverse()(direct_vector_),
              AlmostEquals(Vector<Quotient<Length, Amount>, MirrorWorld>(
                               {-0.2 * Metre / Mole,
                                -0.4 * Metre / Mole,
                                -0.6 * Metre / Mole}),
                           2));
  EXPECT_THAT(conformal_b_.Inverse()(direct_vector_),
              AlmostEquals(Vector<Quotient<Length, Amount>, DirectWorld>(
                               {0.2 * Metre / Mole,
                                0.4 * Metre / Mole,
                                -0.6 * Metre / Mole}),
                  1,
                  2));
}

TEST_F(ConformalMapTest, Composition) {
  auto const conformal_ac = conformal_a_ * conformal_c_.Inverse();
  EXPECT_THAT(conformal_ac(direct_vector_),
              AlmostEquals(Vector<Length, DirectWorld>(
                  R3Element<Length>(-2.0 * Metre,
                                                1.0 * Metre,
                                                3.0 * Metre)), 1, 6));
  EXPECT_THAT((conformal_b_ * conformal_a_).Determinant(),
              Eq(-Pow<6>(5 * Mole)));
  EXPECT_THAT((conformal_a_ * conformal_c_.Inverse()).
              Determinant(), Eq(1));
  EXPECT_THAT((conformal_b_ * conformal_c_).Determinant(),
              Eq(-Pow<6>(5 * Mole)));
}

TEST_F(ConformalMapDeathTest, SerializationError) {
  Identity<DirectWorld, DirectWorld> id;
  EXPECT_DEATH({
    serialization::LinearMap message;
    id.WriteToMessage(&message);
    DirectConf const o = DirectConf::ReadFromMessage(message);
  }, "HasExtension.*ConformalMap");
}

TEST_F(ConformalMapTest, SerializationSuccess) {
  serialization::LinearMap message;
  conformal_a_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_from_frame());
  EXPECT_TRUE(message.has_to_frame());
  EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
            message.to_frame().tag_type_fingerprint());
  EXPECT_NE(message.from_frame().tag(),
            message.to_frame().tag());
  EXPECT_EQ(message.from_frame().is_inertial(),
            message.to_frame().is_inertial());
  EXPECT_TRUE(message.HasExtension(
      serialization::ConformalMap::extension));
  serialization::ConformalMap const& extension =
      message.GetExtension(serialization::ConformalMap::extension);
  EXPECT_THAT(extension.orthogonal_map().quaternion().real_part(),
              AlmostEquals(0.5, 1));
  EXPECT_EQ(
      0.5,
      extension.orthogonal_map().quaternion().imaginary_part().x().double_());
  EXPECT_EQ(
      0.5,
      extension.orthogonal_map().quaternion().imaginary_part().y().double_());
  EXPECT_EQ(
      0.5,
      extension.orthogonal_map().quaternion().imaginary_part().z().double_());
  MirrorConf const o = MirrorConf::ReadFromMessage(message);
  EXPECT_EQ(conformal_a_(mirror_vector_), o(mirror_vector_));
}

TEST_F(ConformalMapTest, Output) {
  std::cout << conformal_a_ << "\n";
}

}  // namespace geometry
}  // namespace principia
