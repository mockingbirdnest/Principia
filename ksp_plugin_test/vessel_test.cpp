
#include "ksp_plugin/vessel.hpp"

#include <limits>
#include <list>
#include <memory>
#include <set>
#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin_test/plugin_io.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rotating_body.hpp"
#include "physics/mock_ephemeris.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace ksp_plugin {

using base::not_null;
using base::make_not_null_unique;
using geometry::Barycentre;
using geometry::Bivector;
using geometry::Displacement;
using geometry::InertiaTensor;
using geometry::InfiniteFuture;
using geometry::Instant;
using geometry::Position;
using geometry::R3x3Matrix;
using geometry::Vector;
using geometry::Velocity;
using interface::ReadPluginFromFile;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::MassiveBody;
using physics::MockEphemeris;
using physics::RigidMotion;
using physics::RotatingBody;
using quantities::Force;
using quantities::Mass;
using quantities::MomentOfInertia;
using quantities::Pow;
using quantities::Torque;
using quantities::si::Degree;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::EqualsProto;
using testing_utilities::AppendTrajectoryTimeline;
using testing_utilities::NewAcceleratedTrajectoryTimeline;
using testing_utilities::NewCircularTrajectoryTimeline;
using testing_utilities::NewLinearTrajectoryTimeline;
using ::testing::AnyNumber;
using ::testing::DoAll;
using ::testing::ElementsAre;
using ::testing::MockFunction;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::_;

class VesselTest : public testing::Test {
 protected:
  VesselTest()
      : body_(MassiveBody::Parameters(1 * Kilogram),
              RotatingBody<Barycentric>::Parameters(
                  /*mean_radius=*/1 * Metre,
                  /*reference_angle=*/0 * Degree,
                  /*reference_instant=*/t0_,
                  /*angular_frequency=*/1 * Radian / Second,
                  /*right_ascension_of_pole=*/0 * Degree,
                  /*declination_of_pole=*/90 * Degree)),
        celestial_(&body_),
        inertia_tensor1_(MakeWaterSphereInertiaTensor(mass1_)),
        inertia_tensor2_(MakeWaterSphereInertiaTensor(mass2_)),
        vessel_("123",
                "vessel",
                &celestial_,
                &ephemeris_,
                DefaultPredictionParameters(),
                DefaultDownsamplingParameters()) {
    auto p1 = make_not_null_unique<Part>(
        part_id1_,
        "p1",
        mass1_,
        EccentricPart::origin,
        inertia_tensor1_,
        RigidMotion<EccentricPart, Barycentric>::MakeNonRotatingMotion(p1_dof_),
        /*deletion_callback=*/nullptr);
    auto p2 = make_not_null_unique<Part>(
        part_id2_,
        "p2",
        mass2_,
        EccentricPart::origin,
        inertia_tensor2_,
        RigidMotion<EccentricPart, Barycentric>::MakeNonRotatingMotion(p2_dof_),
        /*deletion_callback=*/nullptr);
    p1_ = p1.get();
    p2_ = p2.get();
    vessel_.AddPart(std::move(p1));
    vessel_.AddPart(std::move(p2));
  }

  bool IsCollapsible() const {
    return vessel_.IsCollapsible();
  }

  void WriteCheckpointToMessage(
      not_null<serialization::Vessel*> const message) const {
    vessel_.checkpointer_->WriteToMessage(message->mutable_checkpoint());
  }

  MockEphemeris<Barycentric> ephemeris_;
  RotatingBody<Barycentric> const body_;
  Celestial const celestial_;
  PartId const part_id1_ = 111;
  PartId const part_id2_ = 222;
  Mass const mass1_ = 1 * Kilogram;
  Mass const mass2_ = 2 * Kilogram;
  InertiaTensor<RigidPart> inertia_tensor1_;
  InertiaTensor<RigidPart> inertia_tensor2_;

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

  Instant const t0_;
  Part* p1_;
  Part* p2_;
  Vessel vessel_;
};

TEST_F(VesselTest, Parent) {
  Celestial other_celestial(&body_);
  EXPECT_EQ(&celestial_, vessel_.parent());
  vessel_.set_parent(&other_celestial);
  EXPECT_EQ(&other_celestial, vessel_.parent());
}

TEST_F(VesselTest, KeepAndFreeParts) {
  std::set<PartId> remaining_part_ids;
  vessel_.ForAllParts([&remaining_part_ids](Part const& part) {
    remaining_part_ids.insert(part.part_id());
  });
  EXPECT_THAT(remaining_part_ids, ElementsAre(part_id1_, part_id2_));
  EXPECT_EQ(part_id1_, vessel_.part(part_id1_)->part_id());
  EXPECT_EQ(part_id2_, vessel_.part(part_id2_)->part_id());
  remaining_part_ids.clear();

  vessel_.KeepPart(part_id2_);
  vessel_.FreeParts();
  vessel_.ForAllParts([&remaining_part_ids](Part const& part) {
    remaining_part_ids.insert(part.part_id());
  });
  EXPECT_THAT(remaining_part_ids, ElementsAre(part_id2_));
  EXPECT_EQ(part_id2_, vessel_.part(part_id2_)->part_id());
}

TEST_F(VesselTest, PrepareHistory) {
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 2 * Second));
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 2 * Second, _, _))
      .Times(AnyNumber());
  vessel_.CreateHistoryIfNeeded(t0_ + 1 * Second);

  auto const expected_dof = Barycentre<DegreesOfFreedom<Barycentric>, Mass>(
      {p1_dof_, p2_dof_}, {mass1_, mass2_});

  EXPECT_EQ(1, vessel_.psychohistory()->size());
  EXPECT_EQ(t0_ + 1 * Second,
            vessel_.psychohistory()->back().time);
  EXPECT_THAT(
      vessel_.psychohistory()->back().degrees_of_freedom,
      Componentwise(AlmostEquals(expected_dof.position(), 0),
                    AlmostEquals(expected_dof.velocity(), 8)));
}

TEST_F(VesselTest, AdvanceTime) {
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 2 * Second));
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 2 * Second, _, _))
      .Times(AnyNumber());
  vessel_.CreateHistoryIfNeeded(t0_);

  AppendTrajectoryTimeline<Barycentric>(
      NewLinearTrajectoryTimeline<Barycentric>(p1_dof_,
                                               /*Δt=*/0.5 * Second,
                                               /*t0=*/t0_,
                                               /*t1=*/t0_ + 0.5 * Second,
                                               /*t2=*/t0_ + 1.5 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p1_->AppendToHistory(time, degrees_of_freedom);
      });
  AppendTrajectoryTimeline<Barycentric>(
      NewLinearTrajectoryTimeline<Barycentric>(p2_dof_,
                                               /*Δt=*/0.5 * Second,
                                               /*t0=*/t0_,
                                               /*t1=*/t0_ + 0.5 * Second,
                                               /*t2=*/t0_ + 1.5 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p2_->AppendToHistory(time, degrees_of_freedom);
      });

  vessel_.AdvanceTime();

  auto const expected_vessel_psychohistory = NewLinearTrajectoryTimeline(
      Barycentre<DegreesOfFreedom<Barycentric>, Mass>({p1_dof_, p2_dof_},
                                                      {mass1_, mass2_}),
      /*Δt=*/0.5 * Second,
      /*t1=*/t0_,
      /*t2=*/t0_ + 1.1 * Second);

  EXPECT_EQ(3, vessel_.history()->size() + vessel_.psychohistory()->size() - 1);
  auto it1 = vessel_.history()->begin();
  auto it2 = expected_vessel_psychohistory.begin();
  for (;
       it1 != vessel_.psychohistory()->end() &&
       it2 != expected_vessel_psychohistory.end();
       ++it1, ++it2) {
    EXPECT_EQ(it1->time, it2->time);
    EXPECT_THAT(
        it1->degrees_of_freedom,
        Componentwise(AlmostEquals(it2->degrees_of_freedom.position(), 0, 1),
                      AlmostEquals(it2->degrees_of_freedom.velocity(), 0, 8)));
  }
}

TEST_F(VesselTest, Prediction) {
  EXPECT_CALL(ephemeris_, t_min_locked())
      .WillRepeatedly(Return(t0_));
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 2 * Second));

  // The call to fill the prognostication until t_max.
  auto const expected_vessel_prediction = NewLinearTrajectoryTimeline(
      Barycentre<DegreesOfFreedom<Barycentric>, Mass>({p1_dof_, p2_dof_},
                                                      {mass1_, mass2_}),
      /*Δt=*/0.5 * Second,
      /*t1=*/t0_,
      /*t2=*/t0_ + 2 * Second);
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 2 * Second, _, _))
      .WillOnce(DoAll(
          AppendPointsToDiscreteTrajectory(&expected_vessel_prediction),
          Return(absl::OkStatus())))
      .WillRepeatedly(Return(absl::OkStatus()));

  // The call to extend the exphemeris.  Irrelevant since we won't be looking at
  // these points.
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .WillRepeatedly(Return(absl::OkStatus()));

  vessel_.CreateHistoryIfNeeded(t0_);
  // Polling for the integration to happen.
  do {
    vessel_.RefreshPrediction(t0_ + 1 * Second);
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(100ms);
  } while (vessel_.prediction()->back().time == t0_);

  EXPECT_EQ(3, vessel_.prediction()->size());
  auto it1 = vessel_.prediction()->begin();
  auto it2 = expected_vessel_prediction.begin();
  for (;
       it1 != vessel_.prediction()->end() &&
       it2 != expected_vessel_prediction.end();
       ++it1, ++it2) {
    EXPECT_EQ(it1->time, it2->time);
    EXPECT_THAT(
        it1->degrees_of_freedom,
        Componentwise(AlmostEquals(it2->degrees_of_freedom.position(), 0, 0),
                      AlmostEquals(it2->degrees_of_freedom.velocity(), 0, 8)));
  }
}

TEST_F(VesselTest, PredictBeyondTheInfinite) {
  EXPECT_CALL(ephemeris_, t_min_locked())
      .WillRepeatedly(Return(t0_));
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 5 * Second));

  // The call to fill the prognostication until t_max.
  auto const expected_vessel_prediction1 = NewLinearTrajectoryTimeline(
      Barycentre<DegreesOfFreedom<Barycentric>, Mass>({p1_dof_, p2_dof_},
                                                      {mass1_, mass2_}),
      /*Δt=*/0.5 * Second,
      /*t1=*/t0_,
      /*t2=*/t0_ + 5.5 * Second);
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 5 * Second, _, _))
      .WillOnce(DoAll(
          AppendPointsToDiscreteTrajectory(&expected_vessel_prediction1),
          Return(absl::OkStatus())))
      .WillRepeatedly(Return(absl::OkStatus()));

  // The call to extend the exphemeris by many points.
  auto const expected_vessel_prediction2 = NewLinearTrajectoryTimeline(
      Barycentre<DegreesOfFreedom<Barycentric>, Mass>({p1_dof_, p2_dof_},
                                                      {mass1_, mass2_}),
      /*Δt=*/0.5 * Second,
      /*t1=*/t0_ + 5.5 * Second,
      /*t2=*/t0_ + FlightPlan::max_ephemeris_steps_per_frame * Second);
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .WillOnce(DoAll(
          AppendPointsToDiscreteTrajectory(&expected_vessel_prediction2),
          Return(absl::OkStatus())))
      .WillRepeatedly(Return(absl::OkStatus()));

  vessel_.CreateHistoryIfNeeded(t0_);
  // Polling for the integration to happen.
  do {
    vessel_.RefreshPrediction();
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(100ms);
  } while (vessel_.prediction()->size() <
           expected_vessel_prediction1.size() +
               expected_vessel_prediction2.size());

  auto it = expected_vessel_prediction1.begin();
  for (auto const& [time, degrees_of_freedom] : *vessel_.prediction()) {
    EXPECT_EQ(time, it->time);
    EXPECT_THAT(
        degrees_of_freedom,
        Componentwise(AlmostEquals(it->degrees_of_freedom.position(), 0, 0),
                      AlmostEquals(it->degrees_of_freedom.velocity(), 0, 8)));
    if (it->time == t0_ + 5 * Second) {
      it = expected_vessel_prediction2.begin();
    } else {
      ++it;
    }
  }
}

TEST_F(VesselTest, FlightPlan) {
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 2 * Second));
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 2 * Second, _, _))
      .Times(AnyNumber());
  std::vector<not_null<MassiveBody const*>> const bodies;
  ON_CALL(ephemeris_, bodies()).WillByDefault(ReturnRef(bodies));
  vessel_.CreateHistoryIfNeeded(t0_);

  EXPECT_FALSE(vessel_.has_flight_plan());
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, t0_ + 3 * Second, _, _))
      .WillOnce(Return(absl::OkStatus()));
  vessel_.CreateFlightPlan(t0_ + 3.0 * Second,
                           10 * Kilogram,
                           DefaultPredictionParameters(),
                           DefaultBurnParameters());
  EXPECT_TRUE(vessel_.has_flight_plan());
  EXPECT_EQ(0, vessel_.flight_plan().number_of_manœuvres());
  EXPECT_EQ(1, vessel_.flight_plan().number_of_segments());
  vessel_.DeleteFlightPlan();
  EXPECT_FALSE(vessel_.has_flight_plan());
}

TEST_F(VesselTest, IsCollapsible) {
  {
    auto const pile_up =
        std::make_shared<PileUp>(/*parts=*/std::list<not_null<Part*>>{p1_, p2_},
                                 Instant{},
                                 DefaultPsychohistoryParameters(),
                                 DefaultHistoryParameters(),
                                 &ephemeris_,
                                 /*deletion_callback=*/nullptr);
    p1_->set_containing_pile_up(pile_up);
    p2_->set_containing_pile_up(pile_up);

    // Same pile-up.
    EXPECT_TRUE(IsCollapsible());

    // Force.
    p1_->apply_intrinsic_force(
        Vector<Force, Barycentric>({1 * Newton, 0 * Newton, 0 * Newton}));
    EXPECT_FALSE(IsCollapsible());

    // Torque.
    p1_->clear_intrinsic_force();
    p1_->apply_intrinsic_torque(
        Bivector<Torque, Barycentric>({0 * Newton * Metre * Radian,
                                       1 * Newton * Metre * Radian,
                                       0 * Newton * Metre * Radian}));
    EXPECT_TRUE(IsCollapsible());
  }

  {
    auto p3 = make_not_null_unique<Part>(
        333,
        "p3",
        mass2_,
        EccentricPart::origin,
        inertia_tensor2_,
        RigidMotion<EccentricPart, Barycentric>::MakeNonRotatingMotion(p2_dof_),
        /*deletion_callback=*/nullptr);
    auto const pile_up = std::make_shared<PileUp>(
        /*parts=*/std::list<not_null<Part*>>{p1_, p2_, p3.get()},
        Instant{},
        DefaultPsychohistoryParameters(),
        DefaultHistoryParameters(),
        &ephemeris_,
        /*deletion_callback=*/nullptr);
    p1_->reset_containing_pile_up();
    p2_->reset_containing_pile_up();
    p1_->set_containing_pile_up(pile_up);
    p2_->set_containing_pile_up(pile_up);
    p3->set_containing_pile_up(pile_up);

    // A pile-up with an extra part.
    EXPECT_FALSE(IsCollapsible());
  }
}

TEST_F(VesselTest, Checkpointing) {
  MockFunction<int(not_null<PileUp const*>)>
      serialization_index_for_pile_up;
  EXPECT_CALL(serialization_index_for_pile_up, Call(_))
      .Times(2)
      .WillRepeatedly(Return(0));

  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 30 * Second));
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 30 * Second, _, _))
      .Times(AnyNumber());
  vessel_.CreateHistoryIfNeeded(t0_);

  auto const pile_up =
      std::make_shared<PileUp>(/*parts=*/std::list<not_null<Part*>>{p1_, p2_},
                                Instant{},
                                DefaultPsychohistoryParameters(),
                                DefaultHistoryParameters(),
                                &ephemeris_,
                                /*deletion_callback=*/nullptr);
  p1_->set_containing_pile_up(pile_up);
  p2_->set_containing_pile_up(pile_up);

  // Free-fall trajectory.  This creates a checkpoint because the history is not
  // collapsible at the beginning.
  AppendTrajectoryTimeline<Barycentric>(
      NewLinearTrajectoryTimeline<Barycentric>(p1_dof_,
                                               /*Δt=*/1 * Second,
                                               /*t1=*/t0_ + 1 * Second,
                                               /*t2=*/t0_ + 11 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p1_->AppendToHistory(time, degrees_of_freedom);
      });
  AppendTrajectoryTimeline<Barycentric>(
      NewLinearTrajectoryTimeline<Barycentric>(p2_dof_,
                                               /*Δt=*/1 * Second,
                                               /*t1=*/t0_ + 1 * Second,
                                               /*t2=*/t0_ + 11 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p2_->AppendToHistory(time, degrees_of_freedom);
      });

  vessel_.DetectCollapsibilityChange();
  vessel_.AdvanceTime();

  // Apply a force.  This segment is not collapsible.
  auto const p1_force =
      Vector<Force, Barycentric>({1 * Newton, 0 * Newton, 0 * Newton});
  p1_->apply_intrinsic_force(p1_force);
  AppendTrajectoryTimeline<Barycentric>(
      NewAcceleratedTrajectoryTimeline(p1_dof_,
                                       /*acceleration=*/p1_force / mass1_,
                                       /*Δt=*/1 * Second,
                                       /*t1=*/t0_ + 11 * Second,
                                       /*t2=*/t0_ + 26 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p1_->AppendToHistory(time, degrees_of_freedom);
      });
  AppendTrajectoryTimeline<Barycentric>(
      NewLinearTrajectoryTimeline<Barycentric>(p2_dof_,
                                               /*Δt=*/1 * Second,
                                               /*t1=*/t0_ + 11 * Second,
                                               /*t2=*/t0_ + 26 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p2_->AppendToHistory(time, degrees_of_freedom);
      });

  vessel_.DetectCollapsibilityChange();
  vessel_.AdvanceTime();

  // Remove the force.  This creates a checkpoint because we closed a non-
  // collapsible segment.
  p1_->clear_intrinsic_force();
  AppendTrajectoryTimeline<Barycentric>(
      NewLinearTrajectoryTimeline<Barycentric>(p1_dof_,
                                               /*Δt=*/1 * Second,
                                               /*t1=*/t0_ + 26 * Second,
                                               /*t2=*/t0_ + 31 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p1_->AppendToHistory(time, degrees_of_freedom);
      });
  AppendTrajectoryTimeline<Barycentric>(
      NewLinearTrajectoryTimeline<Barycentric>(p2_dof_,
                                               /*Δt=*/1 * Second,
                                               /*t1=*/t0_ + 26 * Second,
                                               /*t2=*/t0_ + 31 * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p2_->AppendToHistory(time, degrees_of_freedom);
      });

  vessel_.DetectCollapsibilityChange();
  vessel_.AdvanceTime();

  serialization::Vessel message;
  vessel_.WriteToMessage(&message,
                         serialization_index_for_pile_up.AsStdFunction());
  EXPECT_TRUE(message.has_is_collapsible());
  EXPECT_TRUE(message.is_collapsible());
  EXPECT_EQ(2, message.checkpoint_size());
  {
    auto const& checkpoint = message.checkpoint(0);
    EXPECT_EQ(0, checkpoint.time().scalar().magnitude());
    EXPECT_EQ(3, checkpoint.non_collapsible_segment().segment_size());
    auto const& segment0 = checkpoint.non_collapsible_segment().segment(0);
    EXPECT_EQ(1, segment0.number_of_dense_points());
    EXPECT_EQ(1, segment0.zfp().timeline_size());
    EXPECT_EQ(0, segment0.exact(0).instant().scalar().magnitude());
    auto const& segment1 = checkpoint.non_collapsible_segment().segment(1);
    // The |number_of_dense_points| is 0 for the psychohistory and the
    // prediction because they are not subject to downsampling.
    EXPECT_FALSE(segment1.has_number_of_dense_points());
    EXPECT_EQ(1, segment1.zfp().timeline_size());
    EXPECT_EQ(0, segment1.exact(0).instant().scalar().magnitude());
    auto const& segment2 = checkpoint.non_collapsible_segment().segment(2);
    EXPECT_FALSE(segment2.has_number_of_dense_points());
    EXPECT_EQ(1, segment2.zfp().timeline_size());
    EXPECT_EQ(0, segment2.exact(0).instant().scalar().magnitude());
    EXPECT_EQ(
        1,
        checkpoint.non_collapsible_segment().segment_by_left_endpoint_size());
    auto const& segment_by_left_endpoint0 =
        checkpoint.non_collapsible_segment().segment_by_left_endpoint(0);
    EXPECT_EQ(
        0, segment_by_left_endpoint0.left_endpoint().scalar().magnitude());
    EXPECT_EQ(2, segment_by_left_endpoint0.segment());
  }
  {
    auto const& checkpoint = message.checkpoint(1);
    EXPECT_EQ(25, checkpoint.time().scalar().magnitude());
    EXPECT_EQ(5, checkpoint.non_collapsible_segment().segment_size());
    auto const& segment0 = checkpoint.non_collapsible_segment().segment(0);
    EXPECT_EQ(0, segment0.number_of_dense_points());
    EXPECT_EQ(0, segment0.zfp().timeline_size());
    auto const& segment1 = checkpoint.non_collapsible_segment().segment(1);
    EXPECT_EQ(1, segment1.number_of_dense_points());
    EXPECT_EQ(1, segment1.zfp().timeline_size());
    EXPECT_EQ(10, segment1.exact(0).instant().scalar().magnitude());
    auto const& segment2 = checkpoint.non_collapsible_segment().segment(2);
    EXPECT_EQ(16, segment2.number_of_dense_points());
    EXPECT_EQ(16, segment2.zfp().timeline_size());
    EXPECT_EQ(10, segment2.exact(0).instant().scalar().magnitude());
    auto const& segment3 = checkpoint.non_collapsible_segment().segment(3);
    EXPECT_FALSE(segment3.has_number_of_dense_points());
    EXPECT_EQ(1, segment3.zfp().timeline_size());
    EXPECT_EQ(25, segment3.exact(0).instant().scalar().magnitude());
    auto const& segment4 = checkpoint.non_collapsible_segment().segment(4);
    EXPECT_FALSE(segment4.has_number_of_dense_points());
    EXPECT_EQ(1, segment4.zfp().timeline_size());
    EXPECT_EQ(25, segment4.exact(0).instant().scalar().magnitude());
    EXPECT_EQ(
        2,
        checkpoint.non_collapsible_segment().segment_by_left_endpoint_size());
    auto const& segment_by_left_endpoint0 =
        checkpoint.non_collapsible_segment().segment_by_left_endpoint(0);
    EXPECT_EQ(
        10, segment_by_left_endpoint0.left_endpoint().scalar().magnitude());
    EXPECT_EQ(2, segment_by_left_endpoint0.segment());
    auto const& segment_by_left_endpoint1 =
        checkpoint.non_collapsible_segment().segment_by_left_endpoint(1);
    EXPECT_EQ(
        25, segment_by_left_endpoint1.left_endpoint().scalar().magnitude());
    EXPECT_EQ(4, segment_by_left_endpoint1.segment());
  }
}

TEST_F(VesselTest, SerializationSuccess) {
  MockFunction<int(not_null<PileUp const*>)>
      serialization_index_for_pile_up;
  EXPECT_CALL(serialization_index_for_pile_up, Call(_)).Times(0);

  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 2 * Second));
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 2 * Second, _, _))
      .Times(AnyNumber());
  vessel_.CreateHistoryIfNeeded(t0_);

  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 3 * Second, _, _))
      .WillRepeatedly(Return(absl::OkStatus()));

  std::vector<not_null<MassiveBody const*>> const bodies;
  ON_CALL(ephemeris_, bodies()).WillByDefault(ReturnRef(bodies));

  vessel_.CreateFlightPlan(t0_ + 3.0 * Second,
                           10 * Kilogram,
                           DefaultPredictionParameters(),
                           DefaultBurnParameters());

  serialization::Vessel message;
  vessel_.WriteToMessage(&message,
                         serialization_index_for_pile_up.AsStdFunction());
  EXPECT_TRUE(message.has_history());
  EXPECT_TRUE(message.has_flight_plan());

  EXPECT_CALL(ephemeris_, Prolong(_)).Times(2);
  auto const v = Vessel::ReadFromMessage(
      message, &celestial_, &ephemeris_, /*deletion_callback=*/nullptr);
  EXPECT_TRUE(v->has_flight_plan());

  serialization::Vessel second_message;
  v->WriteToMessage(&second_message,
                    serialization_index_for_pile_up.AsStdFunction());
  EXPECT_THAT(message, EqualsProto(second_message));
}

#if !defined(_DEBUG)
TEST_F(VesselTest, TailSerialization) {
  // Must be large enough that truncation happens.
  // TODO(phl): Don't hard-wire numbers.
  constexpr std::int64_t number_of_points = 100'000;

  MockFunction<int(not_null<PileUp const*>)>
      serialization_index_for_pile_up;
  EXPECT_CALL(serialization_index_for_pile_up, Call(_))
      .Times(2)
      .WillRepeatedly(Return(0));

  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(t0_ + 30 * Second));
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(ephemeris_,
              FlowWithAdaptiveStep(_, _, t0_ + 30 * Second, _, _))
      .Times(AnyNumber());
  vessel_.CreateHistoryIfNeeded(t0_);

  auto const pile_up =
      std::make_shared<PileUp>(/*parts=*/std::list<not_null<Part*>>{p1_, p2_},
                                Instant{},
                                DefaultPsychohistoryParameters(),
                                DefaultHistoryParameters(),
                                &ephemeris_,
                                /*deletion_callback=*/nullptr);
  p1_->set_containing_pile_up(pile_up);
  p2_->set_containing_pile_up(pile_up);

  // A long trajectory for each part.
  AppendTrajectoryTimeline<Barycentric>(
      NewCircularTrajectoryTimeline<Barycentric>(
          /*period=*/20 * Second,
          /*r=*/101 * Metre,
          /*Δt=*/1 * Second,
          /*t1=*/t0_ + 1 * Second,
          /*t2=*/t0_ + number_of_points * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p1_->AppendToHistory(time, degrees_of_freedom);
      });
  AppendTrajectoryTimeline<Barycentric>(
      NewCircularTrajectoryTimeline<Barycentric>(
          /*period=*/20 * Second,
          /*r=*/102 * Metre,
          /*Δt=*/1 * Second,
          /*t1=*/t0_ + 1 * Second,
          /*t2=*/t0_ + number_of_points * Second),
      [this](Instant const& time,
             DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
        p2_->AppendToHistory(time, degrees_of_freedom);
      });

  vessel_.DetectCollapsibilityChange();
  vessel_.AdvanceTime();
  EXPECT_EQ(12'569,
            std::distance(vessel_.history()->begin(),
                          vessel_.psychohistory()->begin()));

  serialization::Vessel message;
  vessel_.WriteToMessage(&message,
                         serialization_index_for_pile_up.AsStdFunction());

  EXPECT_EQ(4, message.history().segment_size());
  {
    // Non-collapsible segment of the history, entirely excluded.
    auto const& segment0 = message.history().segment(0);
    EXPECT_EQ(0, segment0.number_of_dense_points());
    EXPECT_EQ(0, segment0.zfp().timeline_size());
  }
  {
    // Collapsible segment of the history (backstory), truncated to the left.
    auto const& segment1 = message.history().segment(1);
    EXPECT_EQ(79, segment1.number_of_dense_points());
    EXPECT_EQ(t0_ + 4'553 * Second,
              Instant::ReadFromMessage(segment1.exact(0).instant()));
    EXPECT_EQ(t0_ + (number_of_points - 1) * Second,
              Instant::ReadFromMessage(segment1.exact(1).instant()));
    EXPECT_EQ(12'000, segment1.zfp().timeline_size());
  }
  {
    // Psychohistory, only one point.
    auto const& segment2 = message.history().segment(2);
    EXPECT_EQ(0, segment2.number_of_dense_points());
    EXPECT_EQ(1, segment2.zfp().timeline_size());
  }
  {
    // Prediction, excluded except for its first point.
    auto const& segment3 = message.history().segment(3);
    EXPECT_EQ(0, segment3.number_of_dense_points());
    EXPECT_EQ(1, segment3.zfp().timeline_size());
  }

  auto const v = Vessel::ReadFromMessage(
      message, &celestial_, &ephemeris_, /*deletion_callback=*/nullptr);
  EXPECT_TRUE(v->history()->empty());
  auto const backstory = std::next(v->trajectory().segments().begin());
  EXPECT_EQ(t0_ + 4'553 * Second, backstory->front().time);
  EXPECT_EQ(t0_ + (number_of_points - 1) * Second, backstory->back().time);
}
#endif

TEST_F(VesselTest, Reanimator) {
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      R"(C:\Users\phl.mantegna\Desktop\entwurf.proto.b64)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");

  auto const vessel = plugin->GetVessel("2d8d4b86-c80a-4aca-8886-8b9cff0953aa");
  EXPECT_EQ("Entwurf", vessel->name());
  // TODO(phl): Probably incorrect, probably confused by empty segments.
  EXPECT_EQ(11'958, vessel->trajectory().size());
  EXPECT_EQ(t0_ + 21'047'221.7658279724 * Second,
            vessel->trajectory().front().time);
  EXPECT_EQ(t0_ + 22'848'117.5886942074 * Second,
            vessel->psychohistory()->back().time);
  for (auto const& segment : vessel->trajectory().segments()) {
    LOG(ERROR) << segment.size();
    if (!segment.empty()) {
      LOG(ERROR) << segment.front().time << " "<<segment.back().time;
    }
  }

  // Wait for reanimation to happen.
  LOG(ERROR) << "Waiting until Herbert West is done...";
  vessel->WaitForReanimation(t0_);
  LOG(ERROR) << "Herbert West is finally done.";

  // TODO(phl): Check that everything is good.
}

}  // namespace ksp_plugin
}  // namespace principia
