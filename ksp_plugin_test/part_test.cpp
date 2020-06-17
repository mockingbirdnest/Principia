
#include "ksp_plugin/part.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using geometry::Displacement;
using geometry::R3x3Matrix;
using quantities::Force;
using quantities::MomentOfInertia;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using ::testing::_;
using ::testing::MockFunction;
using testing_utilities::AlmostEquals;
using testing_utilities::EqualsProto;

class PartTest : public testing::Test {
 protected:
  PartTest()
      : inertia_tensor_(MakeWaterSphereInertiaTensor(mass_)),
        part_(part_id_,
              "part",
              mass_,
              EccentricPart::origin,
              inertia_tensor_,
              RigidMotion<EccentricPart, Barycentric>::MakeNonRotatingMotion(
                  degrees_of_freedom_),
              /*deletion_callback=*/nullptr) {
    part_.apply_intrinsic_force(intrinsic_force_);
    part_.AppendToHistory(
        astronomy::J2000,
        {Barycentric::origin +
             Displacement<Barycentric>({11 * Metre, 22 * Metre, 33 * Metre}),
         Velocity<Barycentric>(
             {44 * Metre / Second, 55 * Metre / Second, 66 * Metre / Second})});
  }

  DegreesOfFreedom<Barycentric> const degrees_of_freedom_ = {
      Barycentric::origin +
          Displacement<Barycentric>({1 * Metre, 2 * Metre, 3 * Metre}),
      Velocity<Barycentric>(
          {4 * Metre / Second, 5 * Metre / Second, 6 * Metre / Second})};
  PartId const part_id_ = 666;
  Mass const mass_ = 7 * Kilogram;
  Vector<Force, Barycentric> const intrinsic_force_ =
      Vector<Force, Barycentric>({8 * Newton, 9 * Newton, 10 * Newton});
  InertiaTensor<RigidPart> inertia_tensor_;
  Part part_;
};

TEST_F(PartTest, Serialization) {
  MockFunction<int(not_null<not_null<PileUp const*>>)>
      serialization_index_for_pile_up;
  EXPECT_CALL(serialization_index_for_pile_up, Call(_)).Times(0);

  serialization::Part message;
  part_.WriteToMessage(&message,
                       serialization_index_for_pile_up.AsStdFunction());
  EXPECT_EQ(part_id_, message.part_id());
  EXPECT_TRUE(message.has_inertia_tensor());
  EXPECT_TRUE(message.has_mass());
  EXPECT_EQ(7, message.mass().magnitude());
  EXPECT_TRUE(message.has_intrinsic_force());
  EXPECT_TRUE(message.intrinsic_force().has_vector());
  EXPECT_EQ(8, message.intrinsic_force().vector().x().quantity().magnitude());
  EXPECT_EQ(9, message.intrinsic_force().vector().y().quantity().magnitude());
  EXPECT_EQ(10, message.intrinsic_force().vector().z().quantity().magnitude());
  EXPECT_TRUE(message.has_rigid_motion());
  EXPECT_TRUE(message.rigid_motion().has_rigid_transformation());
  EXPECT_TRUE(message.rigid_motion().rigid_transformation().has_to_origin());
  EXPECT_TRUE(message.rigid_motion()
                  .rigid_transformation()
                  .to_origin()
                  .has_multivector());
  EXPECT_TRUE(message.rigid_motion()
                  .rigid_transformation()
                  .to_origin()
                  .multivector()
                  .has_vector());
  EXPECT_EQ(1,
            message.rigid_motion()
                .rigid_transformation()
                .to_origin()
                .multivector()
                .vector()
                .x()
                .quantity()
                .magnitude());
  EXPECT_EQ(2,
            message.rigid_motion()
                .rigid_transformation()
                .to_origin()
                .multivector()
                .vector()
                .y()
                .quantity()
                .magnitude());
  EXPECT_EQ(3,
            message.rigid_motion()
                .rigid_transformation()
                .to_origin()
                .multivector()
                .vector()
                .z()
                .quantity()
                .magnitude());
  EXPECT_TRUE(
      message.rigid_motion().velocity_of_to_frame_origin().has_vector());
  EXPECT_THAT(message.rigid_motion()
                  .velocity_of_to_frame_origin()
                  .vector()
                  .x()
                  .quantity()
                  .magnitude(),
              AlmostEquals(-4, 2));
  EXPECT_THAT(message.rigid_motion()
                  .velocity_of_to_frame_origin()
                  .vector()
                  .y()
                  .quantity()
                  .magnitude(),
              AlmostEquals(-6, 2));
  EXPECT_THAT(message.rigid_motion()
                  .velocity_of_to_frame_origin()
                  .vector()
                  .z()
                  .quantity()
                  .magnitude(),
              AlmostEquals(-5, 2));
  EXPECT_EQ(1, message.prehistory().zfp().timeline_size());
  EXPECT_EQ(1, message.prehistory().children_size());
  EXPECT_EQ(1, message.prehistory().children(0).trajectories_size());
  EXPECT_EQ(
      1,
      message.prehistory().children(0).trajectories(0).zfp().timeline_size());

  auto const p = Part::ReadFromMessage(message, /*deletion_callback=*/nullptr);
  EXPECT_EQ(part_.inertia_tensor(), p->inertia_tensor());
  EXPECT_EQ(part_.intrinsic_force(), p->intrinsic_force());
  EXPECT_EQ(part_.degrees_of_freedom(), p->degrees_of_freedom());

  serialization::Part second_message;
  p->WriteToMessage(&second_message,
                    serialization_index_for_pile_up.AsStdFunction());
  EXPECT_THAT(message, EqualsProto(second_message));
}

}  // namespace internal_part
}  // namespace ksp_plugin
}  // namespace principia
