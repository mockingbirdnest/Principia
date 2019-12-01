
#include "ksp_plugin/part.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "physics/inertia_tensor.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using geometry::Displacement;
using geometry::R3x3Matrix;
using physics::InertiaTensor;
using quantities::Force;
using quantities::MomentOfInertia;
using quantities::SIUnit;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using ::testing::_;
using ::testing::MockFunction;
using testing_utilities::EqualsProto;

class PartTest : public testing::Test {
 protected:
  PartTest()
      : inertia_tensor_(
            InertiaTensor<RigidPart>::MakeWaterSphereInertiaTensor(mass_)),
        part_(part_id_,
              "part",
              inertia_tensor_,
              RigidMotion<RigidPart, Barycentric>::MakeNonRotatingMotion(
                  degrees_of_freedom_),
              /*deletion_callback=*/nullptr) {
    part_.increment_intrinsic_force(intrinsic_force_);
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
  EXPECT_TRUE(message.inertia_tensor().has_mass());
  EXPECT_EQ(7, message.inertia_tensor().mass().magnitude());
  EXPECT_TRUE(message.has_intrinsic_force());
  EXPECT_TRUE(message.intrinsic_force().has_vector());
  EXPECT_EQ(8, message.intrinsic_force().vector().x().quantity().magnitude());
  EXPECT_EQ(9, message.intrinsic_force().vector().y().quantity().magnitude());
  EXPECT_EQ(10, message.intrinsic_force().vector().z().quantity().magnitude());
  EXPECT_TRUE(message.has_degrees_of_freedom());
  EXPECT_TRUE(message.degrees_of_freedom().t1().has_point());
  EXPECT_TRUE(message.degrees_of_freedom().t1().point().has_multivector());
  EXPECT_TRUE(
      message.degrees_of_freedom().t1().point().multivector().has_vector());
  EXPECT_EQ(1,
            message.degrees_of_freedom()
                .t1()
                .point()
                .multivector()
                .vector()
                .x()
                .quantity()
                .magnitude());
  EXPECT_EQ(2,
            message.degrees_of_freedom()
                .t1()
                .point()
                .multivector()
                .vector()
                .y()
                .quantity()
                .magnitude());
  EXPECT_EQ(3,
            message.degrees_of_freedom()
                .t1()
                .point()
                .multivector()
                .vector()
                .z()
                .quantity()
                .magnitude());
  EXPECT_TRUE(message.degrees_of_freedom().t2().has_multivector());
  EXPECT_TRUE(message.degrees_of_freedom().t2().multivector().has_vector());
  EXPECT_EQ(4,
            message.degrees_of_freedom()
                .t2()
                .multivector()
                .vector()
                .x()
                .quantity()
                .magnitude());
  EXPECT_EQ(5,
            message.degrees_of_freedom()
                .t2()
                .multivector()
                .vector()
                .y()
                .quantity()
                .magnitude());
  EXPECT_EQ(6,
            message.degrees_of_freedom()
                .t2()
                .multivector()
                .vector()
                .z()
                .quantity()
                .magnitude());
  EXPECT_EQ(1, message.prehistory().timeline_size());
  EXPECT_EQ(1, message.prehistory().children_size());
  EXPECT_EQ(1, message.prehistory().children(0).trajectories_size());
  EXPECT_EQ(1,
            message.prehistory().children(0).trajectories(0).timeline_size());

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
