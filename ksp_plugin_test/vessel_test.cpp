
#include "ksp_plugin/vessel.hpp"

#include <limits>

#include "astronomy/epoch.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "physics/massive_body.hpp"
#include "physics/mock_ephemeris.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::make_not_null_unique;
using geometry::Displacement;
using geometry::Position;
using geometry::Velocity;
using physics::MassiveBody;
using physics::MockEphemeris;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using ::testing::DoAll;
using ::testing::ElementsAre;
using ::testing::Return;
using ::testing::_;

class VesselTest : public testing::Test {
 protected:
  VesselTest()
      : body_(MassiveBody::Parameters(1 * Kilogram)),
        celestial_(&body_),
        vessel_(&celestial_, &ephemeris_, DefaultPredictionParameters()) {
    auto p1 = make_not_null_unique<Part>(part_id1_,
                                         mass1_,
                                         p1_dof_,
                                         /*deletion_callback=*/nullptr);
    auto p2 = make_not_null_unique<Part>(part_id2_,
                                         mass2_,
                                         p2_dof_,
                                         /*deletion_callback=*/nullptr);
    p1_ = p1.get();
    p2_ = p2.get();
    vessel_.AddPart(std::move(p1));
    vessel_.AddPart(std::move(p2));
  }

  MockEphemeris<Barycentric> ephemeris_;
  MassiveBody const body_;
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
  EXPECT_EQ(celestial_, vessel_.parent());
  vessel_.set_parent(&other_celestial);
  EXPECT_EQ(&other_celestial, vessel_.parent());
}

TEST_F(VesselTest, KeepAndFreeParts) {
  std::set<PartId> remaining_part_ids;
  vessel_.ForAllParts([&remaining_part_ids](Part const& part) {
    remaining_part_ids.insert(part.part_id());
  });
  EXPECT_THAT(remaining_part_ids, ElementsAre(part_id1_, part_id2_));
  EXPECT_EQ(part_id1_, vessel_.part(part_id1_));
  EXPECT_EQ(part_id2_, vessel_.part(part_id2_));

  vessel_.KeepPart(part_id2_);
  vessel_.FreeParts();
  vessel_.ForAllParts([&remaining_part_ids](Part const& part) {
    remaining_part_ids.insert(part.part_id());
  });
  EXPECT_THAT(remaining_part_ids, ElementsAre(part_id1_));
  EXPECT_EQ(part_id1_, vessel_.part(part_id1_));
}

TEST_F(VesselTest, PreparePsychohistory) {
  EXPECT_TRUE(vessel_.psychohistory().Empty());
  vessel_.PreparePsychohistory(astronomy::J2000 + 1 * Second);
  EXPECT_EQ(1, vessel_.psychohistory().Size());
  EXPECT_EQ(astronomy::J2000 + 1 * Second,
            vessel_.psychohistory().End().time());
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
  vessel_.PreparePsychohistory(astronomy::J2000);

  p1_->tail().Append(
      astronomy::J2000 + 0.5 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin + Displacement<Barycentric>(
                                    {1.1 * Metre, 2.1 * Metre, 3.1 * Metre}),
          Velocity<Barycentric>({10.1 * Metre / Second,
                                 20.1 * Metre / Second,
                                 30.1 * Metre / Second})));
  p1_->tail().Append(
      astronomy::J2000 + 1.0 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin + Displacement<Barycentric>(
                                    {1.2 * Metre, 2.2 * Metre, 3.2 * Metre}),
          Velocity<Barycentric>({10.2 * Metre / Second,
                                 20.2 * Metre / Second,
                                 30.2 * Metre / Second})));
  p2_->tail().Append(
      astronomy::J2000 + 0.5 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin + Displacement<Barycentric>(
                                    {6.1 * Metre, 5.1 * Metre, 4.1 * Metre}),
          Velocity<Barycentric>({60.1 * Metre / Second,
                                 50.1 * Metre / Second,
                                 40.1 * Metre / Second})));
  p2_->tail().Append(
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
                                          11.3 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                      {130.3 / 3.0 * Metre / Second,
                                      40.1 * Metre / Second,
                                      110.3 / 3.0 * Metre / Second}), 0)));
  ++it;
  EXPECT_EQ(astronomy::J2000 + 1.0 * Second, it.time());
  EXPECT_THAT(it.degrees_of_freedom(),
              Componentwise(AlmostEquals(Barycentric::origin +
                                      Displacement<Barycentric>(
                                          {13.6 / 3.0 * Metre,
                                          4.2 * Metre,
                                          11.6 / 3.0 * Metre}), 0),
                    AlmostEquals(Velocity<Barycentric>(
                                      {130.6 / 3.0 * Metre / Second,
                                      40.2 * Metre / Second,
                                      110.6 / 3.0 * Metre / Second}), 0)));
}

TEST_F(VesselTest, Prediction) {
  vessel_.PreparePsychohistory(astronomy::J2000);

  EXPECT_CALL(ephemeris_, FlowWithAdaptiveStep(_, _, _, _, _, _))
      .WillOnce(
          DoAll(AppendToDiscreteTrajectory(
                    astronomy::J2000 + 1.0 * Second,
                    DegreesOfFreedom<Barycentric>(
                        Barycentric::origin +
                            Displacement<Barycentric>({14.0 / 3.0 * Metre,
                                                       5.0 * Metre,
                                                       4.0 * Metre}),
                        Velocity<Barycentric>({140.0 / 3.0 * Metre / Second,
                                               50.0 * Metre / Second,
                                               40.0 * Metre / Second}))),
                Return(true)));
  vessel_.UpdatePrediction(astronomy::J2000 + 1 * Second);

  EXPECT_EQ(2, vessel_.prediction().Size());
  auto it = vessel_.prediction().Begin();
  ++it;
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

TEST_F(VesselTest, FlightPlan) {
  vessel_.CreateHistoryAndForkProlongation(t1_, d1_);
  vessel_.AdvanceTimeNotInBubble(t2_);
  EXPECT_FALSE(vessel_.has_flight_plan());
  vessel_.CreateFlightPlan(t3_, 10 * Kilogram, adaptive_parameters_);
  EXPECT_TRUE(vessel_.has_flight_plan());
  EXPECT_EQ(0, vessel_.flight_plan().number_of_manœuvres());
  EXPECT_EQ(1, vessel_.flight_plan().number_of_segments());
  vessel_.DeleteFlightPlan();
  EXPECT_FALSE(vessel_.has_flight_plan());
}

TEST_F(VesselDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Vessel message;
    vessel_.WriteToMessage(&message);
  }, "is_initialized");
  EXPECT_DEATH({
    serialization::Vessel message;
    Vessel::ReadFromMessage(message, ephemeris_.get(), earth_.get());
  }, "message.has_history");
}

TEST_F(VesselTest, SerializationSuccess) {
  serialization::Vessel message;
  vessel_.CreateHistoryAndForkProlongation(t2_, d2_);
  vessel_.AdvanceTimeNotInBubble(t2_);
  vessel_.UpdatePrediction(t3_);
  vessel_.CreateFlightPlan(t3_, 10 * Kilogram, adaptive_parameters_);

  vessel_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_history());
  EXPECT_TRUE(message.has_prediction_fork_time());
  EXPECT_TRUE(message.has_prediction_last_time());
  EXPECT_TRUE(message.has_flight_plan());
  vessel_ = Vessel::ReadFromMessage(message, ephemeris_.get(), earth_.get());
  EXPECT_TRUE(vessel_.is_initialized());
  EXPECT_TRUE(vessel_.has_flight_plan());
}

TEST_F(VesselTest, PredictBeyondTheInfinite) {
  vessel_.CreateHistoryAndForkProlongation(t1_, d1_);
  vessel_.AdvanceTimeNotInBubble(t2_);
  vessel_.set_prediction_adaptive_step_parameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
          /*max_steps=*/5,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second));
  Instant previous_t_max = ephemeris_->t_max();
  for (int i = 0; i < 10; ++i) {
    vessel_.UpdatePrediction(astronomy::InfiniteFuture);
  }
  // We stop prolonging when the ephemeris gets long enough (in this case, a
  // single prolongation suffices).
  EXPECT_THAT(ephemeris_->t_max() - previous_t_max,
              Le(FlightPlan::max_ephemeris_steps_per_frame *
                 ephemeris_fixed_parameters_.step()));
  EXPECT_THAT(vessel_.prediction().last().time(), Lt(ephemeris_->t_max()));

  vessel_.set_prediction_adaptive_step_parameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
          /*max_steps=*/1000,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second));
  previous_t_max = ephemeris_->t_max();
  for (int i = 0; i < 10; ++i) {
    vessel_.UpdatePrediction(astronomy::InfiniteFuture);
  }
  // Here the ephemeris isn't long enough yet; we have prolonged every time.
  EXPECT_THAT((ephemeris_->t_max() - previous_t_max) /
                  (FlightPlan::max_ephemeris_steps_per_frame *
                   ephemeris_fixed_parameters_.step()),
              AllOf(Gt(9), Le(10)));
  EXPECT_THAT(vessel_.prediction().last().time(), Eq(ephemeris_->t_max()));
}

}  // namespace internal_vessel
}  // namespace ksp_plugin
}  // namespace principia
