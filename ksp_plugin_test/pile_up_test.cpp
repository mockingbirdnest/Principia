
#include "ksp_plugin/pile_up.hpp"

#include <limits>
#include <map>
#include <vector>

#include "base/status.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/part.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/mock_integrators.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/mock_ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using base::check_not_null;
using base::make_not_null_unique;
using base::Status;
using geometry::Displacement;
using geometry::Position;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::Vector;
using geometry::Velocity;
using integrators::MockFixedStepSizeIntegrator;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::BlanesMoan2002SRKN6B;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using physics::MockEphemeris;
using physics::RigidMotion;
using quantities::Acceleration;
using quantities::Length;
using quantities::MomentOfInertia;
using quantities::Pow;
using quantities::Speed;
using quantities::SIUnit;
using quantities::Time;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Newton;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::EqualsProto;
using ::testing::ByMove;
using ::testing::DoAll;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::MockFunction;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::_;

// A helper class to expose the internal state of a pile-up for testing.
class TestablePileUp : public PileUp {
 public:
  using PileUp::PileUp;
  using PileUp::DeformPileUpIfNeeded;
  using PileUp::AdvanceTime;
  using PileUp::NudgeParts;

  Mass const& mass() const {
    return mass_;
  }

  Vector<Force, Barycentric> const& intrinsic_force() const {
    return intrinsic_force_;
  }

  not_null<DiscreteTrajectory<Barycentric>*> psychohistory() const {
    return psychohistory_;
  }

  PartTo<RigidMotion<RigidPart, NonRotatingPileUp>> const&
  actual_part_rigid_motion() const {
    return actual_part_rigid_motion_;
  }

  PartTo<RigidMotion<RigidPart, ApparentBubble>> const&
  apparent_part_rigid_motion() const {
    return apparent_part_rigid_motion_;
  }
};

class PileUpTest : public testing::Test {
 protected:
  PileUpTest()
      : inertia_tensor1_(MakeWaterSphereInertiaTensor(mass1_)),
        inertia_tensor2_(MakeWaterSphereInertiaTensor(mass2_)),
        p1_(part_id1_,
            "p1",
            mass1_,
            inertia_tensor1_,
            RigidMotion<RigidPart, Barycentric>::MakeNonRotatingMotion(p1_dof_),
            /*deletion_callback=*/nullptr),
        p2_(part_id2_,
            "p2",
            mass2_,
            inertia_tensor2_,
            RigidMotion<RigidPart, Barycentric>::MakeNonRotatingMotion(p2_dof_),
            /*deletion_callback=*/nullptr) {}

  void CheckPreDeformPileUpInvariants(TestablePileUp& pile_up) {
    EXPECT_EQ(3 * Kilogram, pile_up.mass());

    EXPECT_THAT(
        pile_up.psychohistory()->back().degrees_of_freedom,
        Componentwise(AlmostEquals(Barycentric::origin +
                                       Displacement<Barycentric>(
                                           {13.0 / 3.0 * Metre,
                                            4.0 * Metre,
                                            11.0 / 3.0 * Metre}), 0),
                      AlmostEquals(Velocity<Barycentric>(
                                       {130.0 / 3.0 * Metre / Second,
                                        40.0 * Metre / Second,
                                        110.0 / 3.0 * Metre / Second}), 4)));

    EXPECT_THAT(
        pile_up.actual_part_rigid_motion().at(&p1_)({RigidPart::origin,
                                                     RigidPart::unmoving}),
        Componentwise(AlmostEquals(NonRotatingPileUp::origin +
                                       Displacement<NonRotatingPileUp>(
                                           {-10.0 / 3.0 * Metre,
                                            -2.0 * Metre,
                                            -2.0 / 3.0 * Metre}), 1),
                      AlmostEquals(Velocity<NonRotatingPileUp>(
                                       {-100.0 / 3.0 * Metre / Second,
                                        -20.0 * Metre / Second,
                                        -20.0 / 3.0 * Metre / Second}), 7)));
    EXPECT_THAT(
        pile_up.actual_part_rigid_motion().at(&p2_)({RigidPart::origin,
                                                     RigidPart::unmoving}),
        Componentwise(AlmostEquals(NonRotatingPileUp::origin +
                                       Displacement<NonRotatingPileUp>(
                                           {5.0 / 3.0 * Metre,
                                            1.0 * Metre,
                                            1.0 / 3.0 * Metre}), 3),
                      AlmostEquals(Velocity<NonRotatingPileUp>(
                                       {50.0 / 3.0 * Metre / Second,
                                        10.0 * Metre / Second,
                                        10.0 / 3.0 * Metre / Second}), 17)));

    // Centre of mass of |p1_| and |p2_| in |ApparentBubble|, in SI units:
    //   {1 / 9, -1 / 3, -2 / 9} {10 / 9, -10 / 3, -20 / 9}
    DegreesOfFreedom<ApparentBubble> const p1_dof(
        ApparentBubble::origin +
            Displacement<ApparentBubble>({-11.0 / 3.0 * Metre,
                                          -1.0 * Metre,
                                          2.0 / 3.0 * Metre}),
        Velocity<ApparentBubble>({-110.0 / 3.0 * Metre / Second,
                                  -10.0 * Metre / Second,
                                  20.0 / 3.0 * Metre / Second}));
    DegreesOfFreedom<ApparentBubble> const p2_dof(
        ApparentBubble::origin +
            Displacement<ApparentBubble>({2.0 * Metre,
                                          0.0 * Metre,
                                          -2.0 / 3.0 * Metre}),
        Velocity<ApparentBubble>({20.0 * Metre / Second,
                                  0.0 * Metre / Second,
                                  -20.0 / 3.0 * Metre / Second}));
    pile_up.SetPartApparentRigidMotion(
        &p1_,
        RigidMotion<RigidPart, ApparentBubble>::MakeNonRotatingMotion(p1_dof));
    pile_up.SetPartApparentRigidMotion(
        &p2_,
        RigidMotion<RigidPart, ApparentBubble>::MakeNonRotatingMotion(p2_dof));

    EXPECT_THAT(
        pile_up.apparent_part_rigid_motion().at(&p1_)({RigidPart::origin,
                                                       RigidPart::unmoving}),
        Componentwise(AlmostEquals(ApparentBubble::origin +
                                       Displacement<ApparentBubble>(
                                           {-11.0 / 3.0 * Metre,
                                            -1.0 * Metre,
                                            2.0 / 3.0 * Metre}), 0),
                      AlmostEquals(Velocity<ApparentBubble>(
                                       {-110.0 / 3.0 * Metre / Second,
                                        -10.0 * Metre / Second,
                                        20.0 / 3.0 * Metre / Second}), 4)));
    EXPECT_THAT(
        pile_up.apparent_part_rigid_motion().at(&p2_)({RigidPart::origin,
                                                       RigidPart::unmoving}),
        Componentwise(AlmostEquals(ApparentBubble::origin +
                                       Displacement<ApparentBubble>(
                                           {2.0 * Metre,
                                            0.0 * Metre,
                                            -2.0 / 3.0 * Metre}), 0),
                      AlmostEquals(Velocity<ApparentBubble>(
                                       {20.0 * Metre / Second,
                                        0.0 * Metre / Second,
                                        -20.0 / 3.0 * Metre / Second}), 4)));
  }

  void CheckPreAdvanceTimeInvariants(TestablePileUp& pile_up) {
    EXPECT_THAT(
        pile_up.actual_part_rigid_motion().at(&p1_)({RigidPart::origin,
                                                     RigidPart::unmoving}),
        Componentwise(AlmostEquals(NonRotatingPileUp::origin +
                                       Displacement<NonRotatingPileUp>(
                                           {-34.0 / 9.0 * Metre,
                                            -2.0 / 3.0 * Metre,
                                            8.0 / 9.0 * Metre}), 1),
                      AlmostEquals(Velocity<NonRotatingPileUp>(
                                       {-340.0 / 9.0 * Metre / Second,
                                        -20.0 / 3.0 * Metre / Second,
                                        80.0 / 9.0 * Metre / Second}), 8)));
    EXPECT_THAT(
        pile_up.actual_part_rigid_motion().at(&p2_)({RigidPart::origin,
                                                     RigidPart::unmoving}),
        Componentwise(AlmostEquals(NonRotatingPileUp::origin +
                                       Displacement<NonRotatingPileUp>(
                                           {17.0 / 9.0 * Metre,
                                            1.0 / 3.0 * Metre,
                                            -4.0 / 9.0 * Metre}), 0),
                      AlmostEquals(Velocity<NonRotatingPileUp>(
                                       {170.0 / 9.0 * Metre / Second,
                                        10.0 / 3.0 * Metre / Second,
                                        -40.0 / 9.0 * Metre / Second}), 8)));
    EXPECT_THAT(pile_up.apparent_part_rigid_motion(), IsEmpty());
  }

  MockFunction<void()> deletion_callback_;

  PartId const part_id1_ = 111;
  PartId const part_id2_ = 222;
  Mass const mass1_ = 1 * Kilogram;
  Mass const mass2_ = 2 * Kilogram;
  InertiaTensor<RigidPart> inertia_tensor1_;
  InertiaTensor<RigidPart> inertia_tensor2_;

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

// Exercises the entire lifecycle of a |PileUp| that is subject to an intrinsic
// force.
TEST_F(PileUpTest, LifecycleWithIntrinsicForce) {
  MockEphemeris<Barycentric> ephemeris;
  p1_.apply_intrinsic_force(
      Vector<Force, Barycentric>({1 * Newton, 2 * Newton, 3 * Newton}));
  p2_.apply_intrinsic_force(
      Vector<Force, Barycentric>({11 * Newton, 21 * Newton, 31 * Newton}));
  EXPECT_CALL(deletion_callback_, Call()).Times(1);
  TestablePileUp pile_up({&p1_, &p2_},
                         astronomy::J2000,
                         DefaultPsychohistoryParameters(),
                         DefaultHistoryParameters(),
                         &ephemeris,
                         deletion_callback_.AsStdFunction());
  EXPECT_THAT(pile_up.intrinsic_force(),
              AlmostEquals(Vector<Force, Barycentric>(
                  {12 * Newton, 23 * Newton, 34 * Newton}), 0));

  CheckPreDeformPileUpInvariants(pile_up);

  pile_up.DeformPileUpIfNeeded(astronomy::J2000 + 1 * Second);

  CheckPreAdvanceTimeInvariants(pile_up);

  EXPECT_CALL(ephemeris, FlowWithAdaptiveStep(_, _, _, _, _))
      .WillOnce(DoAll(
          AppendToDiscreteTrajectory(DegreesOfFreedom<Barycentric>(
              Barycentric::origin +
                  Displacement<Barycentric>({1.0 * Metre,
                                             14.0 * Metre,
                                             31.0 / 3.0 * Metre}),
              Velocity<Barycentric>({10.0 * Metre / Second,
                                     140.0 * Metre / Second,
                                     310.0 / 3.0 * Metre / Second}))),
          Return(Status::OK)));
  pile_up.AdvanceTime(astronomy::J2000 + 1 * Second);

  EXPECT_EQ(++p1_.history_begin(), p1_.history_end());
  EXPECT_EQ(p1_.psychohistory_begin(), p1_.psychohistory_end());
  EXPECT_THAT(
      p1_.history_begin()->degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {-25.0 / 9.0 * Metre,
                                          40.0 / 3.0 * Metre,
                                          101.0 / 9.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                     {-250.0 / 9.0 * Metre / Second,
                                      400.0 / 3.0 * Metre / Second,
                                      1010.0 / 9.0 * Metre / Second}), 16)));
  EXPECT_EQ(++p2_.history_begin(), p2_.history_end());
  EXPECT_EQ(p2_.psychohistory_begin(), p2_.psychohistory_end());
  EXPECT_THAT(
      p2_.history_begin()->degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {26.0 / 9.0 * Metre,
                                          43.0 / 3.0 * Metre,
                                          89.0 / 9.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {260.0 / 9.0 * Metre / Second,
                                      430.0 / 3.0 * Metre / Second,
                                      890.0 / 9.0 * Metre / Second}), 8)));
  EXPECT_EQ(1, pile_up.psychohistory()->Size());
  EXPECT_THAT(
      pile_up.psychohistory()->back().degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {1.0 * Metre,
                                          14.0 * Metre,
                                          31.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {10.0 * Metre / Second,
                                      140.0 * Metre / Second,
                                      310.0 / 3.0 * Metre / Second}), 0)));

  pile_up.NudgeParts();

  EXPECT_THAT(
      p1_.degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {-25.0 / 9.0 * Metre,
                                          40.0 / 3.0 * Metre,
                                          101.0 / 9.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                     {-250.0 / 9.0 * Metre / Second,
                                      400.0 / 3.0 * Metre / Second,
                                      1010.0 / 9.0 * Metre / Second}), 20)));
  EXPECT_THAT(
      p2_.degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {26.0 / 9.0 * Metre,
                                          43.0 / 3.0 * Metre,
                                          89.0 / 9.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {260.0 / 9.0 * Metre / Second,
                                      430.0 / 3.0 * Metre / Second,
                                      890.0 / 9.0 * Metre / Second}), 12)));
}

// Same as above, but without an intrinsic force.
TEST_F(PileUpTest, LifecycleWithoutIntrinsicForce) {
  MockEphemeris<Barycentric> ephemeris;
  EXPECT_CALL(deletion_callback_, Call()).Times(1);
  TestablePileUp pile_up({&p1_, &p2_},
                         astronomy::J2000,
                         DefaultPsychohistoryParameters(),
                         DefaultHistoryParameters(),
                         &ephemeris,
                         deletion_callback_.AsStdFunction());
  EXPECT_THAT(pile_up.intrinsic_force(),
              AlmostEquals(Vector<Force, Barycentric>(), 0));

  CheckPreDeformPileUpInvariants(pile_up);

  pile_up.DeformPileUpIfNeeded(astronomy::J2000 + 1 * Second);

  CheckPreAdvanceTimeInvariants(pile_up);

  auto history = pile_up.psychohistory()->parent();
  auto instance = make_not_null_unique<MockFixedStepSizeIntegrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation>::MockInstance>();
  EXPECT_CALL(ephemeris,
              NewInstance(ElementsAre(history), _, _))
      .WillOnce(Return(ByMove(std::move(instance))));
  EXPECT_CALL(ephemeris, FlowWithFixedStep(_, _))
      .WillOnce(DoAll(
          AppendToDiscreteTrajectory(
              &history,
              astronomy::J2000 + 0.4 * Second,
              DegreesOfFreedom<Barycentric>(
                  Barycentric::origin +
                      Displacement<Barycentric>(
                          {1.1 * Metre, 14.1 * Metre, 31.1 / 3.0 * Metre}),
                  Velocity<Barycentric>({10.1 * Metre / Second,
                                         140.1 * Metre / Second,
                                         310.1 / 3.0 * Metre / Second}))),
          AppendToDiscreteTrajectory(
              &history,
              astronomy::J2000 + 0.8 * Second,
              DegreesOfFreedom<Barycentric>(
                  Barycentric::origin +
                      Displacement<Barycentric>(
                          {1.2 * Metre, 14.2 * Metre, 31.2 / 3.0 * Metre}),
                  Velocity<Barycentric>({10.2 * Metre / Second,
                                         140.2 * Metre / Second,
                                         310.2 / 3.0 * Metre / Second}))),
          Return(Status::OK)));
  EXPECT_CALL(ephemeris, FlowWithAdaptiveStep(_, _, _, _, _))
      .WillOnce(DoAll(
          AppendToDiscreteTrajectory(DegreesOfFreedom<Barycentric>(
              Barycentric::origin +
                  Displacement<Barycentric>({1.0 * Metre,
                                             14.0 * Metre,
                                             31.0 / 3.0 * Metre}),
              Velocity<Barycentric>({10.0 * Metre / Second,
                                     140.0 * Metre / Second,
                                     310.0 / 3.0 * Metre / Second}))),
          Return(Status::OK)));
  pile_up.AdvanceTime(astronomy::J2000 + 1 * Second);

  EXPECT_EQ(++(++p1_.history_begin()), p1_.history_end());
  EXPECT_EQ(++p1_.psychohistory_begin(), p1_.psychohistory_end());
  EXPECT_THAT(
      p1_.history_begin()->degrees_of_freedom,
      Componentwise(
          AlmostEquals(Barycentric::origin +
                           Displacement<Barycentric>({-24.1 / 9.0 * Metre,
                                                      40.3 / 3.0 * Metre,
                                                      101.3 / 9.0 * Metre}),
                       1),
          AlmostEquals(Velocity<Barycentric>({-249.1 / 9.0 * Metre / Second,
                                              400.3 / 3.0 * Metre / Second,
                                              1010.3 / 9.0 * Metre / Second}),
                       16)));
  EXPECT_THAT(
      (++p1_.history_begin())->degrees_of_freedom,
      Componentwise(
          AlmostEquals(Barycentric::origin +
                           Displacement<Barycentric>({-23.2 / 9.0 * Metre,
                                                      40.6 / 3.0 * Metre,
                                                      101.6 / 9.0 * Metre}),
                       1),
          AlmostEquals(Velocity<Barycentric>({-248.2 / 9.0 * Metre / Second,
                                              400.6 / 3.0 * Metre / Second,
                                              1010.6 / 9.0 * Metre / Second}),
                       17)));
  EXPECT_THAT(
      p1_.psychohistory_begin()->degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {-25.0 / 9.0 * Metre,
                                          40.0 / 3.0 * Metre,
                                          101.0 / 9.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                     {-250.0 / 9.0 * Metre / Second,
                                      400.0 / 3.0 * Metre / Second,
                                      1010.0 / 9.0 * Metre / Second}), 16)));
  EXPECT_EQ(++(++p2_.history_begin()), p2_.history_end());
  EXPECT_EQ(++p2_.psychohistory_begin(), p2_.psychohistory_end());
  EXPECT_THAT(
      p2_.history_begin()->degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {26.9 / 9.0 * Metre,
                                          43.3 / 3.0 * Metre,
                                          89.3 / 9.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                     {260.9 / 9.0 * Metre / Second,
                                      430.3 / 3.0 * Metre / Second,
                                      890.3 / 9.0 * Metre / Second}), 8)));
  EXPECT_THAT(
      (++p2_.history_begin())->degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {27.8 / 9.0 * Metre,
                                          43.6 / 3.0 * Metre,
                                          89.6 / 9.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                     {261.8 / 9.0 * Metre / Second,
                                      430.6 / 3.0 * Metre / Second,
                                      890.6 / 9.0 * Metre / Second}), 8)));
  EXPECT_THAT(
      p2_.psychohistory_begin()->degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {26.0 / 9.0 * Metre,
                                          43.0 / 3.0 * Metre,
                                          89.0 / 9.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {260.0 / 9.0 * Metre / Second,
                                      430.0 / 3.0 * Metre / Second,
                                      890.0 / 9.0 * Metre / Second}), 8)));
  EXPECT_EQ(2, pile_up.psychohistory()->Size());
  EXPECT_THAT(
      pile_up.psychohistory()->front().degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {1.2 * Metre,
                                          14.2 * Metre,
                                          31.2 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {10.2 * Metre / Second,
                                      140.2 * Metre / Second,
                                      310.2 / 3.0 * Metre / Second}), 0)));
  EXPECT_THAT(
      pile_up.psychohistory()->back().degrees_of_freedom,
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {1.0 * Metre,
                                          14.0 * Metre,
                                          31.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {10.0 * Metre / Second,
                                      140.0 * Metre / Second,
                                      310.0 / 3.0 * Metre / Second}), 0)));

  pile_up.NudgeParts();

  EXPECT_THAT(
      p1_.degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {-25.0 / 9.0 * Metre,
                                          40.0 / 3.0 * Metre,
                                          101.0 / 9.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                     {-250.0 / 9.0 * Metre / Second,
                                      400.0 / 3.0 * Metre / Second,
                                      1010.0 / 9.0 * Metre / Second}), 20)));
  EXPECT_THAT(
      p2_.degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                     Displacement<Barycentric>(
                                         {26.0 / 9.0 * Metre,
                                          43.0 / 3.0 * Metre,
                                          89.0 / 9.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                     {260.0 / 9.0 * Metre / Second,
                                      430.0 / 3.0 * Metre / Second,
                                      890.0 / 9.0 * Metre / Second}), 12)));
}

TEST_F(PileUpTest, MidStepIntrinsicForce) {
  // An empty ephemeris; the parameters don't matter, since there are no bodies
  // to integrate.
  // NOTE(egg): ... except we have to put a body because |Ephemeris| doesn't
  // want to be empty.  We put a tiny one very far.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  bodies.emplace_back(make_not_null_unique<MassiveBody>(1 * Kilogram));
  std::vector<DegreesOfFreedom<Barycentric>> initial_state{
      DegreesOfFreedom<Barycentric>{
          Barycentric::origin +
              Displacement<Barycentric>(
                  {std::pow(2, 100) * Metre, 0 * Metre, 0 * Metre}),
          Barycentric::unmoving}};
  Ephemeris<Barycentric> ephemeris{
      std::move(bodies),
      initial_state,
      /*initial_time=*/astronomy::J2000,
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<Barycentric>::FixedStepParameters{
          SymplecticRungeKuttaNyströmIntegrator<BlanesMoan2002SRKN6B,
                                                Position<Barycentric>>(),
          1 * Second}};

  Time const fixed_step = 10 * Second;
  Ephemeris<Barycentric>::FixedStepParameters fixed_parameters{
      SymplecticRungeKuttaNyströmIntegrator<BlanesMoan2002SRKN6B,
                                            Position<Barycentric>>(),
      fixed_step};
  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_parameters{
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Position<Barycentric>>(),
      /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*length_integration_tolerance*/ 1 * Micro(Metre),
      /*speed_integration_tolerance=*/1 * Micro(Metre) / Second};

  EXPECT_CALL(deletion_callback_, Call()).Times(1);
  TestablePileUp pile_up({&p1_}, astronomy::J2000,
                         DefaultPsychohistoryParameters(),
                         DefaultHistoryParameters(),
                         &ephemeris,
                         deletion_callback_.AsStdFunction());
  Velocity<Barycentric> const old_velocity =
      p1_.degrees_of_freedom().velocity();

  pile_up.AdvanceTime(astronomy::J2000 + 1.5 * fixed_step);
  pile_up.NudgeParts();
  EXPECT_THAT(p1_.degrees_of_freedom().velocity(),
              AlmostEquals(old_velocity, 4));

  Vector<Acceleration, Barycentric> const a{{1729 * Metre / Pow<2>(Second),
                                             -168 * Metre / Pow<2>(Second),
                                             504 * Metre / Pow<2>(Second)}};
  p1_.apply_intrinsic_force(p1_.mass() * a);
  pile_up.RecomputeFromParts();
  pile_up.AdvanceTime(astronomy::J2000 + 2 * fixed_step);
  pile_up.NudgeParts();
  EXPECT_THAT(p1_.degrees_of_freedom().velocity(),
              AlmostEquals(old_velocity + 0.5 * fixed_step * a, 1));
}

TEST_F(PileUpTest, Serialization) {
  MockEphemeris<Barycentric> ephemeris;
  p1_.apply_intrinsic_force(
      Vector<Force, Barycentric>({1 * Newton, 2 * Newton, 3 * Newton}));
  p2_.apply_intrinsic_force(
      Vector<Force, Barycentric>({11 * Newton, 21 * Newton, 31 * Newton}));
  EXPECT_CALL(deletion_callback_, Call()).Times(2);
  TestablePileUp pile_up({&p1_, &p2_},
                         astronomy::J2000,
                         DefaultPsychohistoryParameters(),
                         DefaultHistoryParameters(),
                         &ephemeris,
                         deletion_callback_.AsStdFunction());

  serialization::PileUp message;
  pile_up.WriteToMessage(&message);

  EXPECT_EQ(2, message.part_id_size());
  EXPECT_EQ(part_id1_, message.part_id(0));
  EXPECT_EQ(part_id2_, message.part_id(1));
  EXPECT_EQ(1, message.history().timeline_size());
  EXPECT_EQ(2, message.actual_part_rigid_motion().size());
  EXPECT_TRUE(message.apparent_part_rigid_motion().empty());

  auto const part_id_to_part = [this](PartId const part_id) {
    if (part_id == part_id1_) {
      return &p1_;
    }
    if (part_id == part_id2_) {
      return &p2_;
    }
    LOG(FATAL) << "Unexpected part id " << part_id;
    base::noreturn();
  };
  auto const p = PileUp::ReadFromMessage(message,
                                         part_id_to_part,
                                         &ephemeris,
                                         deletion_callback_.AsStdFunction());

  serialization::PileUp second_message;
  p->WriteToMessage(&second_message);
  EXPECT_THAT(message, EqualsProto(second_message));
}

TEST_F(PileUpTest, SerializationCompatibility) {
  MockEphemeris<Barycentric> ephemeris;
  p1_.apply_intrinsic_force(
      Vector<Force, Barycentric>({1 * Newton, 2 * Newton, 3 * Newton}));
  p2_.apply_intrinsic_force(
      Vector<Force, Barycentric>({11 * Newton, 21 * Newton, 31 * Newton}));
  EXPECT_CALL(deletion_callback_, Call()).Times(2);
  TestablePileUp pile_up({&p1_, &p2_},
                         astronomy::J2000,
                         DefaultPsychohistoryParameters(),
                         DefaultHistoryParameters(),
                         &ephemeris,
                         deletion_callback_.AsStdFunction());

  serialization::PileUp message;
  pile_up.WriteToMessage(&message);

  // Clear the children to simulate pre-Cesàro serialization.
  message.mutable_history()->clear_children();
  EXPECT_EQ(1, message.history().timeline_size());

  auto const part_id_to_part = [this](PartId const part_id) {
    if (part_id == part_id1_) {
      return &p1_;
    }
    if (part_id == part_id2_) {
      return &p2_;
    }
    LOG(FATAL) << "Unexpected part id " << part_id;
    base::noreturn();
  };
  auto const p = PileUp::ReadFromMessage(message,
                                         part_id_to_part,
                                         &ephemeris,
                                         deletion_callback_.AsStdFunction());

  EXPECT_CALL(ephemeris, FlowWithAdaptiveStep(_, _, _, _, _))
      .WillOnce(DoAll(
          AppendToDiscreteTrajectory(DegreesOfFreedom<Barycentric>(
              Barycentric::origin +
                  Displacement<Barycentric>({1.0 * Metre,
                                             14.0 * Metre,
                                             31.0 / 3.0 * Metre}),
              Velocity<Barycentric>({10.0 * Metre / Second,
                                     140.0 * Metre / Second,
                                     310.0 / 3.0 * Metre / Second}))),
          Return(Status::OK)));
  p->DeformAndAdvanceTime(astronomy::J2000 + 1 * Second);
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
