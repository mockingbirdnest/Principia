
#include "ksp_plugin/flight_plan.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::MassiveBody;
using quantities::si::Kilogram;
using quantities::si::Milli;
using quantities::si::Newton;

namespace ksp_plugin {

class FlightPlanTest : public testing::Test {
 protected:
  using TestNavigationFrame =
      BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>;

  FlightPlanTest() {
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    bodies.emplace_back(
        make_not_null_unique<MassiveBody>(
            MassiveBody::Parameters(1 * Pow<3>(Metre) / Pow<2>(Second),
                                    1 * Metre)));
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

using FlightPlanDeathTest = FlightPlanTest;

TEST_F(FlightPlanDeathTest, DestroyingFirstSegment) {
  EXPECT_DEATH({
    root_.ForgetAfter(t0_ - 7 * Second);
  }, "Destroying the first segment of flight plan");
}

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
  EXPECT_EQ(0, flight_plan_->number_of_manœuvres());
  flight_plan_->SetFinalTime(t0_ + 42 * Second);
  EXPECT_TRUE(flight_plan_->Append(first_burn()));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_FALSE(flight_plan_->Append(first_burn()));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_TRUE(flight_plan_->Append(second_burn()));
  EXPECT_EQ(2, flight_plan_->number_of_manœuvres());
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
  EXPECT_EQ(2, flight_plan_->number_of_manœuvres());
  flight_plan_->RemoveLast();
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  flight_plan_->RemoveLast();
  EXPECT_EQ(0, flight_plan_->number_of_manœuvres());
  // Check that appending still works.
  EXPECT_TRUE(flight_plan_->Append(first_burn()));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
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
      flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
          final_mass();
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_FALSE(flight_plan_->ReplaceLast(second_burn()));
  EXPECT_EQ(old_final_mass,
            flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
                final_mass());
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  flight_plan_->SetFinalTime(t0_ + 42 * Second);
  EXPECT_TRUE(flight_plan_->ReplaceLast(second_burn()));
  EXPECT_GT(old_final_mass,
            flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
                final_mass());
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
}

TEST_F(FlightPlanTest, Segments) {
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
  EXPECT_EQ(3, flight_plan_->number_of_segments());
  EXPECT_TRUE(flight_plan_->Append(second_burn()));
  EXPECT_EQ(5, flight_plan_->number_of_segments());

  std::vector<Instant> times;
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;

  int last_times_size = times.size();
  Instant last_t = t0_ - 2 * π * Second;
  for (int i = 0; i < flight_plan_->number_of_segments(); ++i) {
    flight_plan_->GetSegment(i, &begin, &end);
    for (auto it = begin; it != end; ++it) {
      Instant const& t = it.time();
      EXPECT_LE(last_t, t);
      EXPECT_LE(t, t0_ + 42 * Second);
      times.push_back(t);
    }
    EXPECT_LT(last_times_size, times.size());
    last_times_size = times.size();
  }
}

TEST_F(FlightPlanTest, Serialization) {
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

  serialization::FlightPlan message;
  flight_plan_->WriteToMessage(&message);
  EXPECT_TRUE(message.has_initial_mass());
  EXPECT_TRUE(message.has_initial_time());
  EXPECT_TRUE(message.has_final_time());
  EXPECT_TRUE(message.has_length_integration_tolerance());
  EXPECT_TRUE(message.has_speed_integration_tolerance());
  EXPECT_EQ(2, message.manoeuvre_size());

  // We need a copy of |root_| otherwise both flight plans have segments that
  // point into |root_| and both want to destroy the forks.  Might as well do
  // the copy using serialization, since it's how it works in real life.
  serialization::Trajectory serialized_trajectory;
  root_.WriteToMessage(&serialized_trajectory);
  auto const root_read =
      DiscreteTrajectory<Barycentric>::ReadFromMessage(serialized_trajectory);

  std::unique_ptr<FlightPlan> flight_plan_read =
      FlightPlan::ReadFromMessage(message, root_read.get(), ephemeris_.get());
  EXPECT_EQ(t0_ - 2 * π * Second, flight_plan_read->initial_time());
  EXPECT_EQ(t0_ + 42 * Second, flight_plan_read->final_time());
  EXPECT_EQ(2, flight_plan_read->number_of_manœuvres());
  EXPECT_EQ(5, flight_plan_read->number_of_segments());
}

}  // namespace ksp_plugin
}  // namespace principia
