
#include "ksp_plugin/vessel.hpp"

#include <limits>
#include <set>

#include "astronomy/epoch.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/integrators.hpp"
#include "physics/massive_body.hpp"
#include "physics/rotating_body.hpp"
#include "physics/mock_ephemeris.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::make_not_null_unique;
using base::Status;
using geometry::Displacement;
using geometry::Position;
using geometry::Velocity;
using physics::MassiveBody;
using physics::MockEphemeris;
using physics::RotatingBody;
using quantities::si::Degree;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::EqualsProto;
using ::testing::AnyNumber;
using ::testing::DoAll;
using ::testing::ElementsAre;
using ::testing::MockFunction;
using ::testing::Return;
using ::testing::_;

class VesselTest : public testing::Test {
 protected:
  VesselTest()
      : body_(MassiveBody::Parameters(1 * Kilogram),
              RotatingBody<Barycentric>::Parameters(
                  /*mean_radius=*/1 * Metre,
                  /*reference_angle=*/0 * Degree,
                  /*reference_instant=*/astronomy::J2000,
                  /*angular_frequency=*/1 * Radian / Second,
                  /*right_ascension_of_pole=*/0 * Degree,
                  /*declination_of_pole=*/90 * Degree)),
        celestial_(&body_),
        vessel_("123",
                "vessel",
                &celestial_,
                &ephemeris_,
                DefaultPredictionParameters()) {
    auto p1 = make_not_null_unique<Part>(part_id1_,
                                         "p1",
                                         mass1_,
                                         p1_dof_,
                                         /*deletion_callback=*/nullptr);
    auto p2 = make_not_null_unique<Part>(part_id2_,
                                         "p2",
                                         mass2_,
                                         p2_dof_,
                                         /*deletion_callback=*/nullptr);
    p1_ = p1.get();
    p2_ = p2.get();
    vessel_.AddPart(std::move(p1));
    vessel_.AddPart(std::move(p2));
  }

  MockEphemeris<Barycentric> ephemeris_;
  RotatingBody<Barycentric> const body_;
  Celestial const celestial_;
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
      .WillRepeatedly(Return(astronomy::J2000 + 2 * Second));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 2 * Second, _, _))
      .Times(AnyNumber());
  vessel_.PrepareHistory(astronomy::J2000 + 1 * Second);

  EXPECT_EQ(1, vessel_.psychohistory().Size());
  EXPECT_EQ(astronomy::J2000 + 1 * Second,
            vessel_.psychohistory().last().time());
  EXPECT_THAT(
      vessel_.psychohistory().last().degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                      Displacement<Barycentric>(
                                          {13.0 / 3.0 * Metre,
                                          4.0 * Metre,
                                          11.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                      {130.0 / 3.0 * Metre / Second,
                                      40.0 * Metre / Second,
                                      110.0 / 3.0 * Metre / Second}), 0)));
}

TEST_F(VesselTest, AdvanceTime) {
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(astronomy::J2000 + 2 * Second));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 2 * Second, _, _))
      .Times(AnyNumber());
  vessel_.PrepareHistory(astronomy::J2000);

  p1_->AppendToHistory(
      astronomy::J2000 + 0.5 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin + Displacement<Barycentric>(
                                    {1.1 * Metre, 2.1 * Metre, 3.1 * Metre}),
          Velocity<Barycentric>({10.1 * Metre / Second,
                                 20.1 * Metre / Second,
                                 30.1 * Metre / Second})));
  p1_->AppendToHistory(
      astronomy::J2000 + 1.0 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin + Displacement<Barycentric>(
                                    {1.2 * Metre, 2.2 * Metre, 3.2 * Metre}),
          Velocity<Barycentric>({10.2 * Metre / Second,
                                 20.2 * Metre / Second,
                                 30.2 * Metre / Second})));
  p2_->AppendToHistory(
      astronomy::J2000 + 0.5 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin + Displacement<Barycentric>(
                                    {6.1 * Metre, 5.1 * Metre, 4.1 * Metre}),
          Velocity<Barycentric>({60.1 * Metre / Second,
                                 50.1 * Metre / Second,
                                 40.1 * Metre / Second})));
  p2_->AppendToHistory(
      astronomy::J2000 + 1.0 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin + Displacement<Barycentric>(
                                    {6.2 * Metre, 5.2 * Metre, 4.2 * Metre}),
          Velocity<Barycentric>({60.2 * Metre / Second,
                                 50.2 * Metre / Second,
                                 40.2 * Metre / Second})));

  vessel_.AdvanceTime();

  EXPECT_EQ(3, vessel_.psychohistory().Size());
  auto it = vessel_.psychohistory().Begin();
  ++it;
  EXPECT_EQ(astronomy::J2000 + 0.5 * Second, it.time());
  EXPECT_THAT(it.degrees_of_freedom(),
              Componentwise(AlmostEquals(Barycentric::origin +
                                         Displacement<Barycentric>(
                                             {13.3 / 3.0 * Metre,
                                              4.1 * Metre,
                                              11.3 / 3.0 * Metre}), 2),
                    AlmostEquals(Velocity<Barycentric>(
                                      {130.3 / 3.0 * Metre / Second,
                                       40.1 * Metre / Second,
                                       110.3 / 3.0 * Metre / Second}), 1)));
  ++it;
  EXPECT_EQ(astronomy::J2000 + 1.0 * Second, it.time());
  EXPECT_THAT(it.degrees_of_freedom(),
              Componentwise(AlmostEquals(Barycentric::origin +
                                         Displacement<Barycentric>(
                                             {13.6 / 3.0 * Metre,
                                              4.2 * Metre,
                                              11.6 / 3.0 * Metre}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                                      {130.6 / 3.0 * Metre / Second,
                                       40.2 * Metre / Second,
                                       110.6 / 3.0 * Metre / Second}), 0)));
}

TEST_F(VesselTest, Prediction) {
  EXPECT_CALL(ephemeris_, t_min_locked())
      .WillRepeatedly(Return(astronomy::J2000));
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(astronomy::J2000 + 2 * Second));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::InfiniteFuture, _, _))
      .WillOnce(
          DoAll(AppendToDiscreteTrajectory(
                    astronomy::J2000 + 1.0 * Second,
                    DegreesOfFreedom<Barycentric>(
                        Barycentric::origin +
                            Displacement<Barycentric>(
                                {15.0 / 3.0 * Metre, 6.0 * Metre, 5.0 * Metre}),
                        Velocity<Barycentric>({150.0 / 3.0 * Metre / Second,
                                               60.0 * Metre / Second,
                                               50.0 * Metre / Second}))),
                Return(Status::OK)))
      .WillRepeatedly(Return(Status::OK));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 2 * Second, _, _))
      .WillOnce(
          DoAll(AppendToDiscreteTrajectory(
                    astronomy::J2000 + 1.0 * Second,
                    DegreesOfFreedom<Barycentric>(
                        Barycentric::origin +
                            Displacement<Barycentric>(
                                {14.0 / 3.0 * Metre, 5.0 * Metre, 4.0 * Metre}),
                        Velocity<Barycentric>({140.0 / 3.0 * Metre / Second,
                                               50.0 * Metre / Second,
                                               40.0 * Metre / Second}))),
                Return(Status::OK)))
      .WillRepeatedly(Return(Status::OK));

  vessel_.PrepareHistory(astronomy::J2000);
  // Polling for the integration to happen.
  do {
    vessel_.RefreshPrediction(astronomy::J2000 + 1 * Second);
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(100ms);
  } while (vessel_.prediction().last().time() == astronomy::J2000);

  EXPECT_EQ(2, vessel_.prediction().Size());
  auto it = vessel_.prediction().Begin();
  EXPECT_EQ(astronomy::J2000, it.time());
  EXPECT_THAT(
      it.degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                      Displacement<Barycentric>(
                                          {13.0 / 3.0 * Metre,
                                           4.0 * Metre,
                                           11.0 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                      {130.0 / 3.0 * Metre / Second,
                                       40.0 * Metre / Second,
                                       110.0 / 3.0 * Metre / Second}), 0)));
  ++it;
  EXPECT_EQ(astronomy::J2000 + 1.0 * Second, it.time());
  EXPECT_THAT(
      it.degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                      Displacement<Barycentric>(
                                          {14.0 / 3.0 * Metre,
                                           5.0 * Metre,
                                           4.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                      {140.0 / 3.0 * Metre / Second,
                                       50.0 * Metre / Second,
                                       40.0 * Metre / Second}), 0)));
}

TEST_F(VesselTest, PredictBeyondTheInfinite) {
  EXPECT_CALL(ephemeris_, t_min_locked())
      .WillRepeatedly(Return(astronomy::J2000));
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(astronomy::J2000 + 0.5 * Second));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 0.5 * Second, _, _))
      .WillRepeatedly(
          DoAll(AppendToDiscreteTrajectory(
                    astronomy::J2000 + 0.5 * Second,
                    DegreesOfFreedom<Barycentric>(
                        Barycentric::origin +
                            Displacement<Barycentric>({14.0 / 3.0 * Metre,
                                                       5.0 * Metre,
                                                       4.0 * Metre}),
                        Velocity<Barycentric>({140.0 / 3.0 * Metre / Second,
                                               50.0 * Metre / Second,
                                               40.0 * Metre / Second}))),
                Return(Status::OK)));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::InfiniteFuture, _, _))
      .WillRepeatedly(
          DoAll(AppendToDiscreteTrajectory(
                    astronomy::J2000 + 1.0 * Second,
                    DegreesOfFreedom<Barycentric>(
                        Barycentric::origin +
                            Displacement<Barycentric>({5.0 * Metre,
                                                       6.0 * Metre,
                                                       5.0 * Metre}),
                        Velocity<Barycentric>({50.0 * Metre / Second,
                                               60.0 * Metre / Second,
                                               50.0 * Metre / Second}))),
                Return(Status::OK)));
  vessel_.PrepareHistory(astronomy::J2000);
  // Polling for the integration to happen.
  do {
    vessel_.RefreshPrediction();
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(100ms);
  } while (vessel_.prediction().Size() != 3);

  auto it = vessel_.prediction().Begin();
  ++it;
  EXPECT_EQ(astronomy::J2000 + 0.5 * Second, it.time());
  ++it;
  EXPECT_EQ(astronomy::J2000 + 1.0 * Second, it.time());
  EXPECT_THAT(
      it.degrees_of_freedom(),
      Componentwise(AlmostEquals(Barycentric::origin +
                                      Displacement<Barycentric>(
                                          {5.0 * Metre,
                                           6.0 * Metre,
                                           5.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                      {50.0 * Metre / Second,
                                       60.0 * Metre / Second,
                                       50.0 * Metre / Second}), 0)));
}

TEST_F(VesselTest, FlightPlan) {
  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(astronomy::J2000 + 2 * Second));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 2 * Second, _, _))
      .Times(AnyNumber());
  vessel_.PrepareHistory(astronomy::J2000);

  EXPECT_FALSE(vessel_.has_flight_plan());
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 3 * Second, _, _))
      .WillOnce(Return(Status::OK));
  vessel_.CreateFlightPlan(astronomy::J2000 + 3.0 * Second,
                           10 * Kilogram,
                           DefaultPredictionParameters(),
                           DefaultBurnParameters());
  EXPECT_TRUE(vessel_.has_flight_plan());
  EXPECT_EQ(0, vessel_.flight_plan().number_of_manœuvres());
  EXPECT_EQ(1, vessel_.flight_plan().number_of_segments());
  vessel_.DeleteFlightPlan();
  EXPECT_FALSE(vessel_.has_flight_plan());
}

TEST_F(VesselTest, SerializationSuccess) {
  MockFunction<int(not_null<PileUp const*>)>
      serialization_index_for_pile_up;
  EXPECT_CALL(serialization_index_for_pile_up, Call(_)).Times(0);

  EXPECT_CALL(ephemeris_, t_max())
      .WillRepeatedly(Return(astronomy::J2000 + 2 * Second));
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::InfiniteFuture, _, _))
      .Times(AnyNumber());
  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 2 * Second, _, _))
      .Times(AnyNumber());
  vessel_.PrepareHistory(astronomy::J2000);

  EXPECT_CALL(
      ephemeris_,
      FlowWithAdaptiveStep(_, _, astronomy::J2000 + 3 * Second, _, _))
      .WillRepeatedly(Return(Status::OK));
  vessel_.CreateFlightPlan(astronomy::J2000 + 3.0 * Second,
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

}  // namespace internal_vessel
}  // namespace ksp_plugin
}  // namespace principia
