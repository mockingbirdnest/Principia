
#include "ksp_plugin/part.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {

using si::Kilogram;
using si::Metre;
using si::Second;

namespace ksp_plugin {

class PartTest : public testing::Test {
 protected:
  PartTest()
      : degrees_of_freedom_({
            Barycentric::origin +
                Displacement<Barycentric>({1 * Metre, 2 * Metre, 3 * Metre}),
            Velocity<Barycentric>({4 * Metre / Second,
                                   5 * Metre / Second,
                                   6 * Metre / Second})}),
        mass_(7 * Kilogram),
        gravitational_acceleration_to_be_applied_by_ksp_(
            Vector<Acceleration, Barycentric>(
                {8 * Metre / Second / Second,
                 9 * Metre / Second / Second,
                 10 * Metre / Second / Second})),
        part_(degrees_of_freedom_,
              mass_,
              gravitational_acceleration_to_be_applied_by_ksp_) {}

  DegreesOfFreedom<Barycentric> degrees_of_freedom_;
  Mass mass_;
  Vector<Acceleration, Barycentric>
      gravitational_acceleration_to_be_applied_by_ksp_;
  Part<Barycentric> part_;
};

TEST_F(PartTest, Serialization) {
  serialization::Part message;
  part_.WriteToMessage(&message);
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
  EXPECT_TRUE(message.has_mass());
  EXPECT_EQ(7, message.mass().magnitude());
  EXPECT_TRUE(message.has_gravitational_acceleration_to_be_applied_by_ksp());
  EXPECT_TRUE(message.gravitational_acceleration_to_be_applied_by_ksp().
                  has_vector());
  EXPECT_EQ(8, message.gravitational_acceleration_to_be_applied_by_ksp().
                   vector().x().quantity().magnitude());
  EXPECT_EQ(9, message.gravitational_acceleration_to_be_applied_by_ksp().
                   vector().y().quantity().magnitude());
  EXPECT_EQ(10, message.gravitational_acceleration_to_be_applied_by_ksp().
                    vector().z().quantity().magnitude());

  Part<Barycentric> p = Part<Barycentric>::ReadFromMessage(message);
  EXPECT_EQ(part_.degrees_of_freedom(), p.degrees_of_freedom());
  EXPECT_EQ(part_.mass(), p.mass());
  EXPECT_EQ(part_.gravitational_acceleration_to_be_applied_by_ksp(),
            p.gravitational_acceleration_to_be_applied_by_ksp());
}

}  // namespace ksp_plugin
}  // namespace principia
