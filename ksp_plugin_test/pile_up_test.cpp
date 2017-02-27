#include "ksp_plugin/pile_up.hpp"

#include "ksp_plugin/part.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using geometry::Displacement;
using geometry::R3Element;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using quantities::Length;
using quantities::Speed;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using ::testing::Return;
using ::testing::ReturnRef;

class PileUpTest : public testing::Test {
 protected:
  PileUpTest()
      : p1_(part_id1_, mass1_),
        p2_(part_id2_, mass2_) {
    p1_.increment_intrinsic_force(
        Vector<Force, Barycentric>({1 * Newton, 2 * Newton, 3 * Newton}));
    p2_.increment_intrinsic_force(
        Vector<Force, Barycentric>({11 * Newton, 21 * Newton, 31 * Newton}));
    p1_.set_degrees_of_freedom(p1_dof_);
    p2_.set_degrees_of_freedom(p2_dof_);
  }

  using RigidPileUp = PileUp::RigidPileUp;

  PartId const part_id1_ = 111;
  PartId const part_id2_ = 222;
  Mass const mass1_ = 1 * Kilogram;
  Mass const mass2_ = 2 * Kilogram;
  DegreesOfFreedom<Bubble> const p1_dof_ = DegreesOfFreedom<Bubble>(
      Bubble::origin + Displacement<Bubble>({1 * Metre, 2 * Metre, 3 * Metre}),
      Velocity<Bubble>(
          {10 * Metre / Second, 20 * Metre / Second, 30 * Metre / Second}));
  DegreesOfFreedom<Bubble> const p2_dof_ = DegreesOfFreedom<Bubble>(
      Bubble::origin + Displacement<Bubble>({6 * Metre, 5 * Metre, 4 * Metre}),
      Velocity<Bubble>(
          {60 * Metre / Second, 50 * Metre / Second, 40 * Metre / Second}));
  DegreesOfFreedom<Barycentric> const bubble_barycentre =
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin +
              Displacement<Barycentric>({7 * Metre, 8 * Metre, 9 * Metre}),
          Velocity<Barycentric>(
              {70 * Metre / Second, 80 * Metre / Second, 90 * Metre / Second}));

  Part p1_;
  Part p2_;
};

TEST_F(PileUpTest, Lifecycle) {
  PileUp pile_up({&p1_, &p2_}, bubble_barycentre, astronomy::J2000);

  EXPECT_EQ(3 * Kilogram, pile_up.mass_);
  EXPECT_THAT(pile_up.intrinsic_force_,
              AlmostEquals(Vector<Force, Barycentric>(
                  {12 * Newton, 23 * Newton, 34 * Newton}), 0));

  EXPECT_THAT(
      pile_up.psychohistory_.last().degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {34.0 / 3.0 * Metre,
                                          12.0 * Metre,
                                          38.0 / 3.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                     {340.0 / 3.0 * Metre / Second,
                                      120.0 * Metre / Second,
                                      380.0 / 3.0 * Metre / Second}), 1)));

  EXPECT_THAT(
      pile_up.actual_part_degrees_of_freedom_.at(&p1_),
      Componentwise(AlmostEquals(RigidPileUp::origin +
                                     Displacement<RigidPileUp>(
                                         {-10.0 / 3.0 * Metre,
                                          -2.0 * Metre,
                                          -2.0 / 3.0 * Metre}), 1),
                    AlmostEquals(Velocity<RigidPileUp>(
                                     {-100.0 / 3.0 * Metre / Second,
                                      -20.0 * Metre / Second,
                                      -20.0 / 3.0 * Metre / Second}), 3)));
  EXPECT_THAT(
      pile_up.actual_part_degrees_of_freedom_.at(&p2_),
      Componentwise(AlmostEquals(RigidPileUp::origin +
                                     Displacement<RigidPileUp>(
                                         {5.0 / 3.0 * Metre,
                                          1.0 * Metre,
                                          1.0 / 3.0 * Metre}), 3),
                    AlmostEquals(Velocity<RigidPileUp>(
                                     {50.0 / 3.0 * Metre / Second,
                                      10.0 * Metre / Second,
                                      10.0 / 3.0 * Metre / Second}), 5)));

  //pile_up.SetPartApparentDegreesOfFreedom(
  //    &p1_,
  //    DegreesOfFreedom<ApparentBubble>(
  //        ApparentBubble::origin +
  //            Displacement<ApparentBubble>({1 * Metre, 2 * Metre, 3 * Metre}),
  //        Velocity<ApparentBubble>({10 * Metre / Second,
  //                                  20 * Metre / Second,
  //                                  30 * Metre / Second})));
  //pile_up.SetPartApparentDegreesOfFreedom(
  //    &p2_,
  //    DegreesOfFreedom<ApparentBubble>(
  //        ApparentBubble::origin +
  //            Displacement<ApparentBubble>({6 * Metre, 5 * Metre, 4 * Metre}),
  //        Velocity<ApparentBubble>({60 * Metre / Second,
  //                                  50 * Metre / Second,
  //                                  40 * Metre / Second})));

  //pile_up.DeformPileUpIfNeeded();

}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
