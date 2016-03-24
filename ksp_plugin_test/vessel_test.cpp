
#include "ksp_plugin/vessel.hpp"

#include "astronomy/frames.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using physics::Ephemeris;
using physics::SolarSystem;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Second;

namespace ksp_plugin {

// TODO(egg): We would want to use a real ephemeris to properly exercise the
// limit cases.
class VesselTest : public testing::Test {
 protected:
  VesselTest()
      : adaptive_parameters_(
            DormandElMikkawyPrince1986RKN434FM<Position<ICRFJ2000Equator>>(),
            /*max_steps=*/1,
            /*length_integration_tolerance=*/1 * Metre,
            /*speed_integration_tolerance=*/1 * Metre / Second),
        ephemeris_fixed_parameters_(
            McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
            /*step=*/10 * Second),
        history_fixed_parameters_(
            McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
            /*step=*/1 * Second) {
    solar_system_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    t0_ = solar_system_.epoch();
    ephemeris_ = solar_system_.MakeEphemeris(
        /*fitting_tolerance=*/1 * Milli(Metre), ephemeris_fixed_parameters_);
    vessel_ = std::make_unique<Vessel>(
        solar_system_.massive_body(*ephemeris_, "Earth"),
        &ephemeris_,
        adaptive_parameters_,
        history_fixed_parameters_);
  }

  SolarSystem<ICRFJ2000Equator> solar_system_;
  std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  Ephemeris<ICRFJ2000Equator>::AdaptiveStepParameters const
      adaptive_parameters_;
  Ephemeris<ICRFJ2000Equator>::FixedStepParameters const
      ephemeris_fixed_parameters_;
  Ephemeris<ICRFJ2000Equator>::FixedStepParameters const
      history_fixed_parameters_;
  std::unique_ptr<Vessel> vessel_;
  DegreesOfFreedom<ICRFJ2000Equator> d1_ = {
      ICRFJ2000Equator::origin +
          Displacement<ICRFJ2000Equator>({1 * Metre, 2 * Metre, 3 * Metre}),
      Velocity<ICRFJ2000Equator>(
          {4 * Metre / Second, 5 * Metre / Second, 6 * Metre / Second})};
  DegreesOfFreedom<ICRFJ2000Equator> d2_ = {
      ICRFJ2000Equator::origin +
          Displacement<ICRFJ2000Equator>({11 * Metre, 12 * Metre, 13 * Metre}),
      Velocity<ICRFJ2000Equator>(
          {14 * Metre / Second, 15 * Metre / Second, 16 * Metre / Second})};
  Instant t0_;
  Instant const t1_ = kUniversalTimeEpoch;
  Instant const t2_ = kUniversalTimeEpoch + 42.3 * Second;
};

using VesselDeathTest = VesselTest;

TEST_F(VesselDeathTest, Uninitialized) {
  EXPECT_DEATH({
    vessel_->history();
  }, "is_initialized");
  EXPECT_DEATH({
    vessel_->prolongation();
  }, "is_initialized");
}

TEST_F(VesselTest, Initialization) {
  EXPECT_FALSE(vessel_->is_initialized());
  vessel_->CreateHistoryAndForkProlongation(t2_, d2_);
  EXPECT_TRUE(vessel_->is_initialized());
  auto const& prolongation = vessel_->prolongation();
  EXPECT_EQ(t2_, prolongation.last().time());
  auto const& history = vessel_->history();
  EXPECT_EQ(t2_, history.last().time());
  EXPECT_FALSE(vessel_->has_flight_plan());
  EXPECT_FALSE(vessel_->has_prediction());
}

TEST_F(VesselTest, Dirty) {
  EXPECT_FALSE(vessel_->is_dirty());
  vessel_->set_dirty();
  EXPECT_TRUE(vessel_->is_dirty());
}

TEST_F(VesselTest, Parent) {
  MassiveBody body(2 * Kilogram);
  Celestial celestial(&body);

  EXPECT_EQ(&parent_, vessel_->parent());
  vessel_->set_parent(&celestial);
  EXPECT_EQ(&celestial, vessel_->parent());
}

TEST_F(VesselTest, AdvanceTime) {
  vessel_->CreateHistoryAndForkProlongation(t1_, d1_);
  vessel_->AdvanceTimeNotInBubble(t2_);
  EXPECT_EQ(kUniversalTimeEpoch + 42 * Second,
            vessel_->history().last().time());
  EXPECT_EQ(t2_, vessel_->prolongation().last().time());
}

TEST_F(VesselDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Vessel message;
    vessel_->WriteToMessage(&message);
  }, "is_initialized");
  EXPECT_DEATH({
    serialization::Vessel message;
    Vessel::ReadFromMessage(message, &ephemeris_, &parent_);
  }, "message.has_history");
}

TEST_F(VesselTest, SerializationSuccess) {
  serialization::Vessel message;
  EXPECT_FALSE(message.has_history());
  vessel_->CreateHistoryAndForkProlongation(t2_, d2_);
  vessel_->WriteToMessage(&message);
  EXPECT_TRUE(message.has_history());
  vessel_ = Vessel::ReadFromMessage(message, &ephemeris_, &parent_);
  EXPECT_TRUE(vessel_->is_initialized());
}

}  // namespace ksp_plugin
}  // namespace principia
