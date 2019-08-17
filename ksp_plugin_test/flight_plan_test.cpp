
#include "ksp_plugin/flight_plan.hpp"

#include <limits>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_flight_plan {

using base::Error;
using base::make_not_null_unique;
using geometry::Barycentre;
using geometry::Displacement;
using geometry::Position;
using geometry::Velocity;
using integrators::EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::Fine1987RKNG34;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Frenet;
using physics::MassiveBody;
using quantities::Force;
using quantities::Pow;
using quantities::SpecificImpulse;
using quantities::Sqrt;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Newton;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::StatusIs;
using testing_utilities::operator""_⑴;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::MockFunction;

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
        Ephemeris<Barycentric>::AccuracyParameters(
            /*fitting_tolerance=*/1 * Milli(Metre),
            /*geopotential_tolerance=*/0x1p-24),
        Ephemeris<Barycentric>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<Barycentric>>(),
            /*step=*/10 * Minute));
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
        /*initial_mass=*/1 * Kilogram,
        /*initial_time=*/root_.Begin().time(),
        /*initial_degrees_of_freedom=*/root_.Begin().degrees_of_freedom(),
        /*desired_final_time=*/t0_ + 1.5 * Second,
        ephemeris_.get(),
        Ephemeris<Barycentric>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Position<Barycentric>>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Position<Barycentric>>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second));
  }

  NavigationManœuvre::Burn MakeTangentBurn(
      Force const& thrust,
      SpecificImpulse const& specific_impulse,
      Instant const& initial_time,
      Speed const& Δv) {
    NavigationManœuvre::Intensity intensity;
    intensity.Δv = Velocity<Frenet<Navigation>>({Δv,
                                                 0 * Metre / Second,
                                                 0 * Metre / Second});
    NavigationManœuvre::Timing timing;
    timing.initial_time = initial_time;
    return {intensity,
            timing,
            thrust,
            specific_impulse,
            make_not_null_unique<TestNavigationFrame>(*navigation_frame_),
            /*is_inertially_fixed=*/true};
  }

  NavigationManœuvre::Burn MakeFirstBurn() {
    NavigationManœuvre::Intensity intensity;
    intensity.Δv = Velocity<Frenet<Navigation>>({1 * Metre / Second,
                                                 0 * Metre / Second,
                                                 0 * Metre / Second});
    NavigationManœuvre::Timing timing;
    timing.initial_time = t0_ + 1 * Second;
    return {intensity,
            timing,
            /*thrust=*/1 * Newton,
            /*specific_impulse=*/1 * Newton * Second / Kilogram,
            make_not_null_unique<TestNavigationFrame>(*navigation_frame_),
            /*is_inertially_fixed=*/true};
  }

  NavigationManœuvre::Burn MakeSecondBurn() {
    auto burn = MakeFirstBurn();
    *burn.timing.initial_time += 1 * Second;
    return burn;
  }

  NavigationManœuvre::Burn MakeThirdBurn() {
    auto burn = MakeFirstBurn();
    *burn.intensity.Δv *= 10;
    return burn;
  }

  Instant const t0_;
  std::unique_ptr<TestNavigationFrame> navigation_frame_;
  std::unique_ptr<Ephemeris<Barycentric>> ephemeris_;
  DiscreteTrajectory<Barycentric> root_;
  std::unique_ptr<FlightPlan> flight_plan_;
};

using FlightPlanDeathTest = FlightPlanTest;

TEST_F(FlightPlanTest, Singular) {
  // A test mass falling from x₀ = 1 m at vanishing initial speed onto a body
  // with gravitational parameter μ = 1 m³/s².  A singularity occurs for
  // t² = π² (x₀/2)³ / μ.
  auto const& μ = ephemeris_->bodies().back()->gravitational_parameter();
  Length const x0 = 1 * Metre;
  Instant const singularity = t0_ + π * Sqrt(Pow<3>(x0 / 2) / μ);
  flight_plan_.reset();
  root_.ForgetAfter(root_.Begin().time());
  // NOTE(egg): In order for to avoid singular Frenet frames NaNing everything,
  // we offset our test particle by 100 ε.  The resulting system is still
  // extremely stiff, indeed the integrator detects a singularity at the exact
  // same time.  We could avoid doing this if we had absolute direction
  // specification for manœuvres.
  root_.Append(t0_,
               {Barycentric::origin +
                    Displacement<Barycentric>(
                        {x0,
                         100 * std::numeric_limits<double>::epsilon() * Metre,
                         0 * Metre}),
                Velocity<Barycentric>()});
  flight_plan_ = std::make_unique<FlightPlan>(
      /*initial_mass=*/1 * Kilogram,
      /*initial_time=*/root_.last().time(),
      /*initial_degrees_of_freedom=*/root_.last().degrees_of_freedom(),
      /*final_time=*/singularity + 100 * Second,
      ephemeris_.get(),
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<Barycentric>>(),
          /*max_steps=*/1000,
          /*length_integration_tolerance=*/1 * Milli(Metre),
          /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Position<Barycentric>>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second));
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  flight_plan_->GetSegment(0, begin, end);
  DiscreteTrajectory<Barycentric>::Iterator back = end;
  --back;
  EXPECT_THAT(AbsoluteError(singularity, back.time()), Lt(1e-4 * Second));
  // Attempting to put a burn past the singularity fails.
  EXPECT_THAT(
      flight_plan_->Append(
          MakeTangentBurn(/*thrust=*/1 * Newton,
                          /*specific_impulse=*/1 * Newton * Second / Kilogram,
                          /*initial_time=*/singularity + 1 * Milli(Second),
                          /*Δv=*/1 * Metre / Second)),
      StatusIs(integrators::termination_condition::VanishingStepSize));
  EXPECT_EQ(1, flight_plan_->number_of_anomalous_manœuvres());

  // Add another manœuvre and check the status.
  EXPECT_THAT(
      flight_plan_->Append(
          MakeTangentBurn(/*thrust=*/1 * Newton,
                          /*specific_impulse=*/1 * Newton * Second / Kilogram,
                          /*initial_time=*/singularity + 10 * Second,
                          /*Δv=*/1 * Metre / Second)),
      StatusIs(integrators::termination_condition::VanishingStepSize));
  EXPECT_EQ(2, flight_plan_->number_of_anomalous_manœuvres());

  // Check that RemoveLast returns the proper statuses.
  EXPECT_THAT(flight_plan_->RemoveLast(),
              StatusIs(integrators::termination_condition::VanishingStepSize));
  flight_plan_->RemoveLast();

  // The singularity occurs during the burn: we're boosting towards the
  // singularity, so we reach the singularity in less than π / 2√2 s, before the
  // end of the burn which lasts 10 (1 - 1/e) s.
  // The derivation of analytic expression for the time at which we reach the
  // singularity is left as an exercise to the reader.
  EXPECT_THAT(
      flight_plan_->Append(
          MakeTangentBurn(/*thrust=*/1 * Newton,
                          /*specific_impulse=*/1 * Newton * Second / Kilogram,
                          /*initial_time=*/t0_ + 0.5 * Second,
                          /*Δv=*/1 * Metre / Second)),
      StatusIs(integrators::termination_condition::VanishingStepSize));
  EXPECT_EQ(0, flight_plan_->number_of_anomalous_manœuvres());

  flight_plan_->GetSegment(1, begin, end);
  back = end;
  --back;
  EXPECT_THAT(back.time(), Lt(singularity));
  EXPECT_NE(begin, back);
  flight_plan_->GetSegment(2, begin, end);
  back = end;
  --back;
  EXPECT_EQ(begin, back);

  // The singularity occurs after the burn: we're boosting away from the
  // singularity, so we reach the singularity in more than π / 2√2 s, after the
  // end of the burn which lasts (1 - 1/e)/10 s.
  // The proof of existence of the singularity, as well as the derivation of
  // analytic expression for the time at which we reach the singularity, are
  // left as an exercise to the reader.
  EXPECT_THAT(
      flight_plan_->Replace(
          MakeTangentBurn(/*thrust=*/10 * Newton,
                          /*specific_impulse=*/1 * Newton * Second / Kilogram,
                          /*initial_time=*/t0_ + 0.5 * Second,
                          /*Δv=*/-1 * Metre / Second),
          /*index=*/0),
      StatusIs(integrators::termination_condition::VanishingStepSize));
  EXPECT_EQ(0, flight_plan_->number_of_anomalous_manœuvres());

  flight_plan_->GetSegment(1, begin, end);
  back = end;
  --back;
  EXPECT_THAT(back.time(), Eq(t0_ + 0.5 * Second + (1 - 1 / e) / 10 * Second));
  EXPECT_NE(begin, back);
  flight_plan_->GetSegment(2, begin, end);
  back = end;
  --back;
  EXPECT_THAT(back.time(), AllOf(Gt(singularity), Lt(t0_ + 2 * Second)));
  EXPECT_NE(begin, back);
}

TEST_F(FlightPlanTest, Append) {
  // Burn ends after final time.
  EXPECT_THAT(flight_plan_->Append(MakeFirstBurn()),
              StatusIs(FlightPlan::does_not_fit));
  EXPECT_EQ(0, flight_plan_->number_of_manœuvres());
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_THAT(flight_plan_->Append(MakeFirstBurn()),
              StatusIs(FlightPlan::does_not_fit));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_OK(flight_plan_->Append(MakeSecondBurn()));
  EXPECT_EQ(2, flight_plan_->number_of_manœuvres());
}

TEST_F(FlightPlanTest, ForgetBefore) {
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_OK(flight_plan_->Append(MakeSecondBurn()));
  EXPECT_EQ(2, flight_plan_->number_of_manœuvres());
  EXPECT_EQ(5, flight_plan_->number_of_segments());

  // Find the extremities of each segment.
  std::vector<Instant> begin_times;
  std::vector<Instant> last_times;
  for (int i = 0; i < flight_plan_->number_of_segments(); ++i) {
    DiscreteTrajectory<Barycentric>::Iterator begin;
    DiscreteTrajectory<Barycentric>::Iterator end;
    flight_plan_->GetSegment(i, begin, end);
    --end;
    begin_times.push_back(begin.time());
    last_times.push_back(end.time());
  }

  // Do the forgetting.
  MockFunction<void()> non_empty;
  for (int i = 0; i < flight_plan_->number_of_segments(); ++i) {
    flight_plan_->ForgetBefore(begin_times[i], non_empty.AsStdFunction());
    EXPECT_EQ(begin_times[i], flight_plan_->initial_time());
    if (i % 2 == 0) {
      // A coast.
      EXPECT_EQ(5 - i, flight_plan_->number_of_segments());
      EXPECT_EQ((5 - i) / 2, flight_plan_->number_of_manœuvres());
    } else {
      // A burn.
      EXPECT_EQ(6 - i, flight_plan_->number_of_segments());
      EXPECT_EQ((6 - i) / 2, flight_plan_->number_of_manœuvres());
    }

    Instant const mid_time =
        Barycentre<Instant, double>({begin_times[i], last_times[i]}, {1, 1});
    flight_plan_->ForgetBefore(mid_time, non_empty.AsStdFunction());
    if (i % 2 == 0) {
      // A coast.
      EXPECT_LE(mid_time, flight_plan_->initial_time());
      EXPECT_EQ(5 - i, flight_plan_->number_of_segments());
      EXPECT_EQ((5 - i) / 2, flight_plan_->number_of_manœuvres());
    } else {
      // A burn.
      EXPECT_EQ(begin_times[i + 1], flight_plan_->initial_time());
      EXPECT_EQ(4 - i, flight_plan_->number_of_segments());
      EXPECT_EQ((4 - i) / 2, flight_plan_->number_of_manœuvres());
    }

    flight_plan_->ForgetBefore(last_times[i], non_empty.AsStdFunction());
    if (i % 2 == 0) {
      // A coast.
      EXPECT_EQ(last_times[i], flight_plan_->initial_time());
      EXPECT_EQ(5 - i, flight_plan_->number_of_segments());
      EXPECT_EQ((5 - i) / 2, flight_plan_->number_of_manœuvres());
    } else {
      // A burn.
      EXPECT_EQ(begin_times[i + 1], flight_plan_->initial_time());
      EXPECT_EQ(4 - i, flight_plan_->number_of_segments());
      EXPECT_EQ((4 - i) / 2, flight_plan_->number_of_manœuvres());
    }
  }

  // Forget after the end.
  MockFunction<void()> empty;
  EXPECT_CALL(empty, Call()).Times(1);
  flight_plan_->ForgetBefore(last_times.back() + 1 * Second,
                             empty.AsStdFunction());
}

TEST_F(FlightPlanTest, RemoveLast) {
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  EXPECT_OK(flight_plan_->Append(MakeSecondBurn()));
  EXPECT_EQ(2, flight_plan_->number_of_manœuvres());
  flight_plan_->RemoveLast();
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  flight_plan_->RemoveLast();
  EXPECT_EQ(0, flight_plan_->number_of_manœuvres());
  // Check that appending still works.
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
}

TEST_F(FlightPlanTest, Replace) {
  flight_plan_->SetDesiredFinalTime(t0_ + 1.7 * Second);
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  Mass const old_final_mass =
      flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
          final_mass();
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_THAT(flight_plan_->Replace(MakeThirdBurn(), /*index=*/0),
              StatusIs(FlightPlan::does_not_fit));
  EXPECT_EQ(old_final_mass,
            flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
                final_mass());
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  EXPECT_OK(flight_plan_->Replace(MakeThirdBurn(), /*index=*/0));
  EXPECT_GT(old_final_mass,
            flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
                final_mass());
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
}

TEST_F(FlightPlanTest, Segments) {
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  EXPECT_EQ(3, flight_plan_->number_of_segments());
  EXPECT_OK(flight_plan_->Append(MakeSecondBurn()));
  EXPECT_EQ(5, flight_plan_->number_of_segments());

  std::vector<Instant> times;
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;

  int last_times_size = times.size();
  Instant last_t = t0_ - 2 * π * Second;
  for (int i = 0; i < flight_plan_->number_of_segments(); ++i) {
    flight_plan_->GetSegment(i, begin, end);
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

TEST_F(FlightPlanTest, SetAdaptiveStepParameter) {
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  EXPECT_OK(flight_plan_->Append(MakeSecondBurn()));
  EXPECT_EQ(5, flight_plan_->number_of_segments());
  flight_plan_->GetSegment(4, begin, end);
  --end;
  EXPECT_EQ(t0_ + 42 * Second, end.time());

  auto const adaptive_step_parameters =
      flight_plan_->adaptive_step_parameters();
  auto const generalized_adaptive_step_parameters =
      flight_plan_->generalized_adaptive_step_parameters();

  // Reduce |max_steps|.  This causes many segments to become truncated so the
  // call to |SetAdaptiveStepParameters| returns false and the flight plan is
  // unaffected.
  EXPECT_THAT(flight_plan_->SetAdaptiveStepParameters(
        Ephemeris<Barycentric>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Position<Barycentric>>(),
            /*max_steps=*/1,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Position<Barycentric>>(),
            /*max_steps=*/1,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)),
      StatusIs(integrators::termination_condition::ReachedMaximalStepCount));
  EXPECT_EQ(2, flight_plan_->number_of_anomalous_manœuvres());

  flight_plan_->SetAdaptiveStepParameters(adaptive_step_parameters,
                                          generalized_adaptive_step_parameters);

  EXPECT_EQ(5, flight_plan_->number_of_segments());
  flight_plan_->GetSegment(4, begin, end);
  --end;
  EXPECT_EQ(t0_ + 42 * Second, end.time());

  // Increase |max_steps|.  It works.
  EXPECT_OK(flight_plan_->SetAdaptiveStepParameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<Barycentric>>(),
          /*max_steps=*/10000,
          /*length_integration_tolerance=*/1 * Milli(Metre),
          /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Position<Barycentric>>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Milli(Metre),
          /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)));

  EXPECT_EQ(5, flight_plan_->number_of_segments());
  flight_plan_->GetSegment(4, begin, end);
  --end;
  EXPECT_EQ(t0_ + 42 * Second, end.time());
}

TEST_F(FlightPlanTest, GuidedBurn) {
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  auto unguided_burn = MakeFirstBurn();
  unguided_burn.thrust /= 10;
  EXPECT_OK(flight_plan_->Append(std::move(unguided_burn)));
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  DiscreteTrajectory<Barycentric>::Iterator last;
  flight_plan_->GetAllSegments(begin, end);
  last = --end;
  Speed const unguided_final_speed =
      last.degrees_of_freedom().velocity().Norm();
  auto guided_burn = MakeFirstBurn();
  guided_burn.thrust /= 10;
  guided_burn.is_inertially_fixed = false;
  EXPECT_OK(flight_plan_->Replace(std::move(guided_burn), /*index=*/0));
  flight_plan_->GetAllSegments(begin, end);
  last = --end;
  Speed const guided_final_speed = last.degrees_of_freedom().velocity().Norm();
  EXPECT_THAT(guided_final_speed, IsNear(1.40_⑴ * unguided_final_speed));
}

TEST_F(FlightPlanTest, Serialization) {
  flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second);
  EXPECT_OK(flight_plan_->Append(MakeFirstBurn()));
  EXPECT_OK(flight_plan_->Append(MakeSecondBurn()));

  serialization::FlightPlan message;
  flight_plan_->WriteToMessage(&message);
  EXPECT_TRUE(message.has_initial_mass());
  EXPECT_TRUE(message.has_initial_time());
  EXPECT_TRUE(message.has_desired_final_time());
  EXPECT_TRUE(message.has_adaptive_step_parameters());
  EXPECT_TRUE(message.adaptive_step_parameters().has_integrator());
  EXPECT_TRUE(message.adaptive_step_parameters().has_max_steps());
  EXPECT_TRUE(
      message.adaptive_step_parameters().has_length_integration_tolerance());
  EXPECT_TRUE(
      message.adaptive_step_parameters().has_speed_integration_tolerance());
  EXPECT_EQ(2, message.manoeuvre_size());

  std::unique_ptr<FlightPlan> flight_plan_read =
      FlightPlan::ReadFromMessage(message, ephemeris_.get());
  EXPECT_EQ(t0_ - 2 * π * Second, flight_plan_read->initial_time());
  EXPECT_EQ(t0_ + 42 * Second, flight_plan_read->desired_final_time());
  EXPECT_EQ(2, flight_plan_read->number_of_manœuvres());
  EXPECT_EQ(5, flight_plan_read->number_of_segments());
}

}  // namespace internal_flight_plan
}  // namespace ksp_plugin
}  // namespace principia
