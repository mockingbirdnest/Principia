
#include "ksp_plugin/part.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using geometry::Displacement;
using quantities::Force;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;

class PartTest : public testing::Test {
 protected:
  PartTest() : part_(part_id_, mass_, degrees_of_freedom_, /*deletion_callback=*/nullptr) {
    part_.increment_intrinsic_force(intrinsic_force_);
    part_.tail().Append(
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
  Part part_;
};

TEST_F(PartTest, Serialization) {
  serialization::Part message;
  part_.WriteToMessage(&message);
  EXPECT_EQ(part_id_, message.part_id());
  EXPECT_TRUE(message.has_mass());
  EXPECT_EQ(7, message.mass().magnitude());
  EXPECT_TRUE(message.has_intrinsic_force());
  EXPECT_TRUE(message.intrinsic_force().
                  has_vector());
  EXPECT_EQ(8, message.intrinsic_force().
                   vector().x().quantity().magnitude());
  EXPECT_EQ(9, message.intrinsic_force().
                   vector().y().quantity().magnitude());
  EXPECT_EQ(10, message.intrinsic_force().
                    vector().z().quantity().magnitude());
  EXPECT_TRUE(message.has_degrees_of_freedom());
  EXPECT_TRUE(message.degrees_of_freedom().t1().has_point());
  EXPECT_TRUE(message.degrees_of_freedom().t1().point().has_multivector());
  EXPECT_TRUE(message.degrees_of_freedom().t1().
                  point().multivector().has_vector());
  EXPECT_EQ(1, message.degrees_of_freedom().t1().
                   point().multivector().vector().x().quantity().magnitude());
  EXPECT_EQ(2, message.degrees_of_freedom().t1().
                   point().multivector().vector().y().quantity().magnitude());
  EXPECT_EQ(3, message.degrees_of_freedom().t1().
                   point().multivector().vector().z().quantity().magnitude());
  EXPECT_TRUE(message.degrees_of_freedom().t2().has_multivector());
  EXPECT_TRUE(message.degrees_of_freedom().t2().multivector().has_vector());
  EXPECT_EQ(4, message.degrees_of_freedom().t2().
                   multivector().vector().x().quantity().magnitude());
  EXPECT_EQ(5, message.degrees_of_freedom().t2().
                   multivector().vector().y().quantity().magnitude());
  EXPECT_EQ(6, message.degrees_of_freedom().t2().
                   multivector().vector().z().quantity().magnitude());
  EXPECT_EQ(1, message.tail().timeline_size());
  EXPECT_FALSE(message.tail_is_authoritative());

  auto const p = Part::ReadFromMessage(message, /*deletion_callback=*/nullptr);
  EXPECT_EQ(part_.mass(), p->mass());
  EXPECT_EQ(part_.intrinsic_force(), p->intrinsic_force());
  EXPECT_EQ(part_.degrees_of_freedom(), p->degrees_of_freedom());

  serialization::Part second_message;
  p->WriteToMessage(&second_message);
  EXPECT_EQ(message.SerializeAsString(), second_message.SerializeAsString());
}

}  // namespace internal_part
}  // namespace ksp_plugin
}  // namespace principia
