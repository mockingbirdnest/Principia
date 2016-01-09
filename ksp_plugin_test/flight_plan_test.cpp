#include "ksp_plugin/flight_plan.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {

using physics::BodyCentredNonRotatingDynamicFrame;
using quantities::si::Kilogram;
using quantities::si::Newton;

namespace ksp_plugin {

class FlightPlanTest : public testing::Test {
 protected:
  using TestNavigationFrame =
      BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>;

  FlightPlanTest() {
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    bodies.emplace_back(
        make_not_null_unique<MassiveBody>(1 * Pow<3>(Metre) / Pow<2>(Second)));
    std::vector<DegreesOfFreedom<Barycentric>> initial_state{
        {Barycentric::origin, Velocity<Barycentric>()}};
    ephemeris_ = std::make_unique<Ephemeris<Barycentric>>(
        std::move(bodies),
        initial_state,
        /*initial_time=*/t0_ - 2 * π * Second,
        integrators::McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
        /*step=*/1 * Second,
        /*fitting_tolerance=*/1 * Milli(Metre));
    navigation_frame_ = std::make_unique<TestNavigationFrame>(
        ephemeris_.get(),
        ephemeris_->bodies().back());
    root_.Append(t0_ - 2 * π * Second,
                 {Barycentric::origin + Displacement<Barycentric>(
                                            {1 * Metre, 0 * Metre, 0 * Metre}),
                  Velocity<Barycentric>({0 * Metre / Second,
                                         1 * Metre / Second,
                                         0 * Metre / Second})});
    root_.Append(t0_ + 2 * π * Second,
                 {Barycentric::origin + Displacement<Barycentric>(
                                            {1 * Metre, 0 * Metre, 0 * Metre}),
                  Velocity<Barycentric>({0 * Metre / Second,
                                         1 * Metre / Second,
                                         0 * Metre / Second})});
    flight_plan_ = std::make_unique<FlightPlan>(
        &root_,
        /*initial_time=*/t0_,
        /*final_time=*/t0_ + 1.5 * Second,
        /*initial_mass=*/1 * Kilogram,
        ephemeris_.get(),
        integrators::DormandElMikkawyPrince1986RKN434FM<
            Position<Barycentric>>(),
        /*length_integration_tolerance=*/ 1 * Milli(Metre),
        /*speed_integration_tolerance=*/ 1 * Milli(Metre) / Second);
  }

  Instant const t0_;
  std::unique_ptr<TestNavigationFrame> navigation_frame_;
  std::unique_ptr<Ephemeris<Barycentric>> ephemeris_;
  DiscreteTrajectory<Barycentric> root_;
  std::unique_ptr<FlightPlan> flight_plan_;
};

TEST_F(FlightPlanTest, Append) {
  auto const first_burn = [this]() -> Burn {
    return {/*thrust=*/1 * Newton,
            /*specific_impulse=*/1 * Newton * Second / Kilogram,
            make_not_null_unique<TestNavigationFrame>(*navigation_frame_),
            /*initial_time=*/t0_ + 1 * Second,
            Velocity<Frenet<Navigation>>(
                {1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second})};
  };
  auto const second_burn = [this, first_burn]() -> Burn {
    auto burn = first_burn();
    burn.initial_time += 1 * Second;
    return burn;
  };
  // Burn ends after final time.
  EXPECT_FALSE(flight_plan_->Append(first_burn()));
  EXPECT_EQ(0, flight_plan_->size());
  flight_plan_->SetFinalTime(t0_ + 42 * Second);
  EXPECT_TRUE(flight_plan_->Append(first_burn()));
  EXPECT_EQ(1, flight_plan_->size());
  EXPECT_FALSE(flight_plan_->Append(first_burn()));
  EXPECT_EQ(1, flight_plan_->size());
  EXPECT_TRUE(flight_plan_->Append(second_burn()));
  EXPECT_EQ(2, flight_plan_->size());
}

TEST_F(FlightPlanTest, Remove) {
  auto const first_burn = [this]() -> Burn {
    return {/*thrust=*/1 * Newton,
            /*specific_impulse=*/1 * Newton * Second / Kilogram,
            make_not_null_unique<TestNavigationFrame>(*navigation_frame_),
            /*initial_time=*/t0_ + 1 * Second,
            Velocity<Frenet<Navigation>>(
                {1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second})};
  };
  auto const second_burn = [this, first_burn]() -> Burn {
    auto burn = first_burn();
    burn.initial_time += 1 * Second;
    return burn;
  };
  flight_plan_->SetFinalTime(t0_ + 42 * Second);
  EXPECT_TRUE(flight_plan_->Append(first_burn()));
  EXPECT_TRUE(flight_plan_->Append(second_burn()));
  EXPECT_EQ(2, flight_plan_->size());
  flight_plan_->RemoveLast();
  EXPECT_EQ(1, flight_plan_->size());
  flight_plan_->RemoveLast();
  EXPECT_EQ(0, flight_plan_->size());
  // Check that appending still works.
  EXPECT_TRUE(flight_plan_->Append(first_burn()));
  EXPECT_EQ(1, flight_plan_->size());
}

TEST_F(FlightPlanTest, Replace) {
  auto const first_burn = [this]() -> Burn {
    return {/*thrust=*/1 * Newton,
            /*specific_impulse=*/1 * Newton * Second / Kilogram,
            make_not_null_unique<TestNavigationFrame>(*navigation_frame_),
            /*initial_time=*/t0_ + 1 * Second,
            Velocity<Frenet<Navigation>>(
                {1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second})};
  };
  auto const second_burn = [this, first_burn]() -> Burn {
    auto burn = first_burn();
    burn.Δv *= 10;
    return burn;
  };
  flight_plan_->SetFinalTime(t0_ + 1.7 * Second);
  EXPECT_TRUE(flight_plan_->Append(first_burn()));
  Mass const old_final_mass =
      flight_plan_->Get(flight_plan_->size() - 1).final_mass();
  EXPECT_EQ(1, flight_plan_->size());
  EXPECT_FALSE(flight_plan_->ReplaceLast(second_burn()));
  EXPECT_EQ(old_final_mass,
            flight_plan_->Get(flight_plan_->size() - 1).final_mass());
  EXPECT_EQ(1, flight_plan_->size());
  flight_plan_->SetFinalTime(t0_ + 42 * Second);
  EXPECT_TRUE(flight_plan_->ReplaceLast(second_burn()));
  EXPECT_GT(old_final_mass,
            flight_plan_->Get(flight_plan_->size() - 1).final_mass());
  EXPECT_EQ(1, flight_plan_->size());
}

}  // namespace ksp_plugin
}  // namespace principia
