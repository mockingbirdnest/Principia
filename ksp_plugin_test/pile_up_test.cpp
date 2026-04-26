#include "ksp_plugin/pile_up.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "astronomy/epoch.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/part.hpp"
#include "numerics/elementary_functions.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/mock_ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/tensors.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::DoAll;
using ::testing::IsEmpty;
using ::testing::MockFunction;
using ::testing::Return;
using ::testing::_;
using namespace principia::astronomy::_epoch;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_identification;
using namespace principia::ksp_plugin::_integrators;
using namespace principia::ksp_plugin::_part;
using namespace principia::ksp_plugin::_pile_up;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory_segment_iterator;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_mock_ephemeris;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_tensors;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_matchers;

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

  DiscreteTrajectorySegmentIterator<Barycentric> psychohistory() const {
    return psychohistory_;
  }

  PartTo<RigidMotion<RigidPart, NonRotatingPileUp>> const&
  actual_part_rigid_motion() const {
    return actual_part_rigid_motion_;
  }

  PartTo<RigidMotion<RigidPart, Apparent>> const&
  apparent_part_rigid_motion() const {
    return apparent_part_rigid_motion_;
  }
};

class PileUpTest : public testing::Test {
 protected:
  using CorrectedPileUp = Frame<struct CorrectedPileUpTag, NonRotating>;
  using Vessel = Frame<struct VesselTag>;

  PileUpTest()
      : inertia_tensor1_(MakeWaterSphereInertiaTensor(mass1_)),
        inertia_tensor2_(MakeWaterSphereInertiaTensor(mass2_)),
        p1_(part_id1_,
            "p1",
            mass1_,
            EccentricPart::origin,
            inertia_tensor1_,
            RigidMotion<EccentricPart, Barycentric>::MakeNonRotatingMotion(
                p1_dof_),
            /*deletion_callback=*/nullptr),
        p2_(part_id2_,
            "p2",
            mass2_,
            EccentricPart::origin,
            inertia_tensor2_,
            RigidMotion<EccentricPart, Barycentric>::MakeNonRotatingMotion(
                p2_dof_),
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

    // Centre of mass of `p1_` and `p2_` in `Apparent`, in SI units:
    //   {1 / 9, -1 / 3, -2 / 9} {10 / 9, -10 / 3, -20 / 9}
    DegreesOfFreedom<Apparent> const p1_dof(
        Apparent::origin +
            Displacement<Apparent>(
                {-11.0 / 3.0 * Metre, -1.0 * Metre, 2.0 / 3.0 * Metre}),
        Velocity<Apparent>({-110.0 / 3.0 * Metre / Second,
                            -10.0 * Metre / Second,
                            20.0 / 3.0 * Metre / Second}));
    DegreesOfFreedom<Apparent> const p2_dof(
        Apparent::origin + Displacement<Apparent>(
                               {2.0 * Metre, 0.0 * Metre, -2.0 / 3.0 * Metre}),
        Velocity<Apparent>({20.0 * Metre / Second,
                            0.0 * Metre / Second,
                            -20.0 / 3.0 * Metre / Second}));
    pile_up.SetPartApparentRigidMotion(
        &p1_, RigidMotion<RigidPart, Apparent>::MakeNonRotatingMotion(p1_dof));
    pile_up.SetPartApparentRigidMotion(
        &p2_, RigidMotion<RigidPart, Apparent>::MakeNonRotatingMotion(p2_dof));

    EXPECT_THAT(
        pile_up.apparent_part_rigid_motion().at(&p1_)(
            {RigidPart::origin, RigidPart::unmoving}),
        Componentwise(
            AlmostEquals(
                Apparent::origin +
                    Displacement<Apparent>(
                        {-11.0 / 3.0 * Metre, -1.0 * Metre, 2.0 / 3.0 * Metre}),
                0),
            AlmostEquals(Velocity<Apparent>({-110.0 / 3.0 * Metre / Second,
                                             -10.0 * Metre / Second,
                                             20.0 / 3.0 * Metre / Second}),
                         4)));
    EXPECT_THAT(
        pile_up.apparent_part_rigid_motion().at(&p2_)(
            {RigidPart::origin, RigidPart::unmoving}),
        Componentwise(
            AlmostEquals(
                Apparent::origin +
                    Displacement<Apparent>(
                        {2.0 * Metre, 0.0 * Metre, -2.0 / 3.0 * Metre}),
                0),
            AlmostEquals(Velocity<Apparent>({20.0 * Metre / Second,
                                             0.0 * Metre / Second,
                                             -20.0 / 3.0 * Metre / Second}),
                         4)));
  }

  void CheckPreAdvanceTimeInvariants(TestablePileUp& pile_up) {
    EXPECT_THAT(pile_up.actual_part_rigid_motion().at(&p1_)(
                    {RigidPart::origin, RigidPart::unmoving}),
                Componentwise(AlmostEquals(NonRotatingPileUp::origin +
                                               Displacement<NonRotatingPileUp>(
                                                   {-34.0 / 9.0 * Metre,
                                                    -2.0 / 3.0 * Metre,
                                                    8.0 / 9.0 * Metre}),
                                           1122),
                              AlmostEquals(Velocity<NonRotatingPileUp>(
                                               {-340.0 / 9.0 * Metre / Second,
                                                -20.0 / 3.0 * Metre / Second,
                                                80.0 / 9.0 * Metre / Second}),
                                           1372)));
    EXPECT_THAT(pile_up.actual_part_rigid_motion().at(&p2_)(
                    {RigidPart::origin, RigidPart::unmoving}),
                Componentwise(AlmostEquals(NonRotatingPileUp::origin +
                                               Displacement<NonRotatingPileUp>(
                                                   {17.0 / 9.0 * Metre,
                                                    1.0 / 3.0 * Metre,
                                                    -4.0 / 9.0 * Metre}),
                                           1122),
                              AlmostEquals(Velocity<NonRotatingPileUp>(
                                               {170.0 / 9.0 * Metre / Second,
                                                10.0 / 3.0 * Metre / Second,
                                                -40.0 / 9.0 * Metre / Second}),
                                           1373)));
    EXPECT_THAT(pile_up.apparent_part_rigid_motion(), IsEmpty());
  }

  MockFunction<void()> deletion_callback_;

  PartId const part_id1_ = 111;
  PartId const part_id2_ = 222;
  Mass const mass1_ = 1 * Kilogram;
  Mass const mass2_ = 2 * Kilogram;
  InertiaTensor<RigidPart> inertia_tensor1_;
  InertiaTensor<RigidPart> inertia_tensor2_;

  // Centre of mass of `p1_` and `p2_` in `Barycentric`, in SI units:
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

TEST_F(PileUpTest, MidStepIntrinsicForce) {
  // An empty ephemeris; the parameters don't matter, since there are no bodies
  // to integrate.
  // NOTE(egg): ... except we have to put a body because `Ephemeris` doesn't
  // want to be empty.  We put a tiny one very far.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  bodies.emplace_back(make_not_null_unique<MassiveBody>(1 * Kilogram));
  std::vector<DegreesOfFreedom<Barycentric>> const initial_state{
      DegreesOfFreedom<Barycentric>{
          Barycentric::origin +
              Displacement<Barycentric>(
                  {std::pow(2, 100) * Metre, 0 * Metre, 0 * Metre}),
          Barycentric::unmoving}};
  Ephemeris<Barycentric> ephemeris{
      std::move(bodies),
      initial_state,
      /*initial_time=*/J2000,
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<Barycentric>::FixedStepParameters{
          SymplecticRungeKuttaNyströmIntegrator<
              BlanesMoan2002SRKN6B,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          1 * Second}};

  Time const fixed_step = 10 * Second;
  Ephemeris<Barycentric>::FixedStepParameters const fixed_parameters{
      SymplecticRungeKuttaNyströmIntegrator<
          BlanesMoan2002SRKN6B,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      fixed_step};
  Ephemeris<Barycentric>::AdaptiveStepParameters const adaptive_parameters{
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*length_integration_tolerance*/ 1 * Micro(Metre),
      /*speed_integration_tolerance=*/1 * Micro(Metre) / Second};

  EXPECT_CALL(deletion_callback_, Call()).Times(1);
  TestablePileUp pile_up({&p1_}, J2000,
                         DefaultPsychohistoryParameters(),
                         DefaultHistoryParameters(),
                         &ephemeris,
                         deletion_callback_.AsStdFunction());
  Velocity<Barycentric> const old_velocity =
      p1_.rigid_motion()({RigidPart::origin, RigidPart::unmoving}).velocity();

  EXPECT_OK(pile_up.AdvanceTime(J2000 + 1.5 * fixed_step));
  pile_up.NudgeParts();
  EXPECT_THAT(
      p1_.rigid_motion()({RigidPart::origin, RigidPart::unmoving}).velocity(),
      AlmostEquals(old_velocity, 2));

  Vector<Acceleration, Barycentric> const a{{1729 * Metre / Pow<2>(Second),
                                             -168 * Metre / Pow<2>(Second),
                                             504 * Metre / Pow<2>(Second)}};
  p1_.apply_intrinsic_force(p1_.mass() * a);
  pile_up.RecomputeFromParts();
  EXPECT_OK(pile_up.AdvanceTime(J2000 + 2 * fixed_step));
  pile_up.NudgeParts();
  EXPECT_THAT(
      p1_.rigid_motion()({RigidPart::origin, RigidPart::unmoving}).velocity(),
      AlmostEquals(old_velocity + 0.5 * fixed_step * a, 1));
}

TEST_F(PileUpTest, Serialization) {
  MockEphemeris<Barycentric> ephemeris;
  p1_.apply_intrinsic_force(
      Vector<Force, Barycentric>({1 * Newton, 2 * Newton, 3 * Newton}));
  p2_.apply_intrinsic_force(
      Vector<Force, Barycentric>({11 * Newton, 21 * Newton, 31 * Newton}));
  EXPECT_CALL(deletion_callback_, Call()).Times(2);
  TestablePileUp const pile_up({&p1_, &p2_},
                               J2000,
                               DefaultPsychohistoryParameters(),
                               DefaultHistoryParameters(),
                               &ephemeris,
                               deletion_callback_.AsStdFunction());

  serialization::PileUp message;
  pile_up.WriteToMessage(&message);

  EXPECT_EQ(2, message.part_id_size());
  EXPECT_EQ(part_id1_, message.part_id(0));
  EXPECT_EQ(part_id2_, message.part_id(1));
  EXPECT_EQ(2, message.history().segment_size());
  EXPECT_EQ(1, message.history().segment(0).zfp().timeline_size());
  EXPECT_EQ(1, message.history().segment(1).zfp().timeline_size());
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
    std::abort();
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
  TestablePileUp const pile_up({&p1_, &p2_},
                               J2000,
                               DefaultPsychohistoryParameters(),
                               DefaultHistoryParameters(),
                               &ephemeris,
                               deletion_callback_.AsStdFunction());

  serialization::PileUp message;
  pile_up.WriteToMessage(&message);

  // Clear the children to simulate pre-Cesàro serialization.
  message.mutable_history()->clear_children();
  EXPECT_EQ(1, message.history().segment(0).zfp().timeline_size());
  EXPECT_EQ(1, message.history().segment(1).zfp().timeline_size());

  auto const part_id_to_part = [this](PartId const part_id) {
    if (part_id == part_id1_) {
      return &p1_;
    }
    if (part_id == part_id2_) {
      return &p2_;
    }
    LOG(FATAL) << "Unexpected part id " << part_id;
    std::abort();
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
          Return(absl::OkStatus())));
  EXPECT_OK(p->DeformAndAdvanceTime(J2000 + 1 * Second));
}

}  // namespace ksp_plugin
}  // namespace principia
