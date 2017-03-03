#include "ksp_plugin/pile_up.hpp"

#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"  // For the Default...Parameters.
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/mock_ephemeris.hpp"
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
using physics::MockEphemeris;
using quantities::Length;
using quantities::Speed;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using ::testing::DoAll;
using ::testing::IsEmpty;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::_;

namespace {

// TODO(phl): Move to a common library.
ACTION_P(AppendToDiscreteTrajectory, degrees_of_freedom) {
  arg0->Append(arg2, degrees_of_freedom);
}

}  // namespace

class PileUpTest : public testing::Test {
 protected:
  PileUpTest()
      : p1_(part_id1_, mass1_, p1_dof_, /*deletion_callback=*/nullptr),
        p2_(part_id2_, mass2_, p2_dof_, /*deletion_callback=*/nullptr) {
    p1_.increment_intrinsic_force(
        Vector<Force, Barycentric>({1 * Newton, 2 * Newton, 3 * Newton}));
    p2_.increment_intrinsic_force(
        Vector<Force, Barycentric>({11 * Newton, 21 * Newton, 31 * Newton}));
  }

  using RigidPileUp = PileUp::RigidPileUp;

  PartId const part_id1_ = 111;
  PartId const part_id2_ = 222;
  Mass const mass1_ = 1 * Kilogram;
  Mass const mass2_ = 2 * Kilogram;

  // Centre of mass of |p1_| and |p2_| in |Barycentric|, in SI units:
  //   {13 / 3, 4, 11 / 3} {130 / 3, 40, 110 / 3}
  DegreesOfFreedom<Barycentric> const p1_dof_ = DegreesOfFreedom<Barycentric>(
      Barycentric::origin +
          Displacement<Barycentric>({1 * Metre, 2 * Metre, 3 * Metre}),
      Velocity<Barycentric>(
          {10 * Metre / Second, 20 * Metre / Second, 30 * Metre / Second}));
  DegreesOfFreedom<Barycentric> const p2_dof_ = DegreesOfFreedom<Barycentric>(
      Barycentric::origin +
          Displacement<Barycentric>({6 * Metre, 5 * Metre, 4 * Metre}),
      Velocity<Barycentric>(
          {60 * Metre / Second, 50 * Metre / Second, 40 * Metre / Second}));

  Part p1_;
  Part p2_;
};

// Exercises the entire lifecycle of a |PileUp|.  This is an intrusive test that
// checks the internal state of the class.
TEST_F(PileUpTest, Lifecycle) {
  PileUp pile_up({&p1_, &p2_}, astronomy::J2000);

  EXPECT_EQ(3 * Kilogram, pile_up.mass_);
  EXPECT_THAT(pile_up.intrinsic_force_,
              AlmostEquals(Vector<Force, Barycentric>(
                  {12 * Newton, 23 * Newton, 34 * Newton}), 0));

  EXPECT_THAT(
      pile_up.psychohistory_.last().degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {13.0 / 3.0 * Metre,
                                          4.0 * Metre,
                                          11.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {130.0 / 3.0 * Metre / Second,
                                      40.0 * Metre / Second,
                                      110.0 / 3.0 * Metre / Second}), 0)));

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

  // Centre of mass of |p1_| and |p2_| in |ApparentBubble|, in SI units:
  //   {1 / 9, -1 / 3, -2 / 9} {10 / 9, -10 / 3, -20 / 9}
  pile_up.SetPartApparentDegreesOfFreedom(
      &p1_,
      DegreesOfFreedom<ApparentBubble>(
          ApparentBubble::origin +
              Displacement<ApparentBubble>({-11.0 / 3.0 * Metre,
                                            -1.0 * Metre,
                                            2.0 / 3.0 * Metre}),
          Velocity<ApparentBubble>({-110.0 / 3.0 * Metre / Second,
                                    -10.0 * Metre / Second,
                                    20.0 / 3.0 * Metre / Second})));
  pile_up.SetPartApparentDegreesOfFreedom(
      &p2_,
      DegreesOfFreedom<ApparentBubble>(
          ApparentBubble::origin +
              Displacement<ApparentBubble>({2.0 * Metre,
                                            0.0 * Metre,
                                            -2.0 / 3.0 * Metre}),
          Velocity<ApparentBubble>({20.0 * Metre / Second,
                                    0.0 * Metre / Second,
                                    -20.0 / 3.0 * Metre / Second})));

  EXPECT_THAT(
      pile_up.apparent_part_degrees_of_freedom_.at(&p1_),
      Componentwise(AlmostEquals(ApparentBubble::origin +
                                     Displacement<ApparentBubble>(
                                         {-11.0 / 3.0 * Metre,
                                          -1.0 * Metre,
                                          2.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<ApparentBubble>(
                                     {-110.0 / 3.0 * Metre / Second,
                                      -10.0 * Metre / Second,
                                      20.0 / 3.0 * Metre / Second}), 0)));
  EXPECT_THAT(
      pile_up.apparent_part_degrees_of_freedom_.at(&p2_),
      Componentwise(AlmostEquals(ApparentBubble::origin +
                                     Displacement<ApparentBubble>(
                                         {2.0 * Metre,
                                          0.0 * Metre,
                                          -2.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<ApparentBubble>(
                                     {20.0 * Metre / Second,
                                      0.0 * Metre / Second,
                                      -20.0 / 3.0 * Metre / Second}), 0)));

  pile_up.DeformPileUpIfNeeded();

  EXPECT_THAT(
      pile_up.actual_part_degrees_of_freedom_.at(&p1_),
      Componentwise(AlmostEquals(RigidPileUp::origin +
                                     Displacement<RigidPileUp>(
                                         {-34.0 / 9.0 * Metre,
                                          -2.0 / 3.0 * Metre,
                                          8.0 / 9.0 * Metre}), 1),
                    AlmostEquals(Velocity<RigidPileUp>(
                                     {-340.0 / 9.0 * Metre / Second,
                                      -20.0 / 3.0 * Metre / Second,
                                      80.0 / 9.0 * Metre / Second}), 1)));
  EXPECT_THAT(
      pile_up.actual_part_degrees_of_freedom_.at(&p2_),
      Componentwise(AlmostEquals(RigidPileUp::origin +
                                     Displacement<RigidPileUp>(
                                         {17.0 / 9.0 * Metre,
                                          1.0 / 3.0 * Metre,
                                          -4.0 / 9.0 * Metre}), 0),
                    AlmostEquals(Velocity<RigidPileUp>(
                                     {170.0 / 9.0 * Metre / Second,
                                      10.0 / 3.0 * Metre / Second,
                                      -40.0 / 9.0 * Metre / Second}), 0)));
  EXPECT_THAT(pile_up.apparent_part_degrees_of_freedom_, IsEmpty());

  MockEphemeris<Barycentric> ephemeris;
  EXPECT_CALL(ephemeris, FlowWithAdaptiveStep(_, _, _, _, _, _))
      .WillOnce(DoAll(
          AppendToDiscreteTrajectory(DegreesOfFreedom<Barycentric>(
              Barycentric::origin +
                  Displacement<Barycentric>({1.0 * Metre,
                                             14.0 * Metre,
                                             31.0 / 3.0 * Metre}),
              Velocity<Barycentric>({10.0 * Metre / Second,
                                     140.0 * Metre / Second,
                                     310.0 / 3.0 * Metre / Second}))),
          Return(true)));
  pile_up.AdvanceTime(ephemeris,
                      astronomy::J2000 + 1 * Second,
                      DefaultHistoryParameters(),
                      DefaultProlongationParameters());

  EXPECT_EQ(2, p1_.tail().Size());
  EXPECT_EQ(2, p2_.tail().Size());
  EXPECT_EQ(1, pile_up.psychohistory_.Size());
  EXPECT_THAT(
      pile_up.psychohistory_.last().degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {1.0 * Metre,
                                          14.0 * Metre,
                                          31.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {10.0 * Metre / Second,
                                      140.0 * Metre / Second,
                                      310.0 / 3.0 * Metre / Second}), 0)));
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
