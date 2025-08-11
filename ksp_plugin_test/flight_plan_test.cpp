#include "ksp_plugin/flight_plan.hpp"

#include <chrono>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <thread>

#include "astronomy/epoch.hpp"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/integrators.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rotating_body.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::MockFunction;
using namespace principia::astronomy::_epoch;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_integrators;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rotating_body;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics;
using namespace std::chrono_literals;

class FlightPlanTest : public testing::Test {
 protected:
  using TestNavigationFrame =
      BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation>;

  FlightPlanTest() {
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    bodies.emplace_back(make_not_null_unique<RotatingBody<Barycentric>>(
        1 * Pow<3>(Metre) / Pow<2>(Second),
        RotatingBody<Barycentric>::Parameters(
            /*mean_radius=*/1 * Metre,
            /*reference_angle=*/0 * Radian,
            /*reference_instant=*/J2000,
            /*angular_frequency=*/1 * Radian / Second,
            /*right_ascension_of_pole=*/0 * Radian,
            /*declination_of_pole=*/0 * Radian)));
    std::vector<DegreesOfFreedom<Barycentric>> initial_state{
        {Barycentric::origin, Barycentric::unmoving}};
    ephemeris_ = std::make_unique<Ephemeris<Barycentric>>(
        std::move(bodies),
        initial_state,
        /*initial_time=*/t0_ - 2 * π * Second,
        Ephemeris<Barycentric>::AccuracyParameters(
            /*fitting_tolerance=*/1 * Milli(Metre),
            /*geopotential_tolerance=*/0x1p-24),
        Ephemeris<Barycentric>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*step=*/10 * Minute));
    EXPECT_OK(ephemeris_->Prolong(t0_ - 2 * π * Second));
    navigation_frame_ = std::make_unique<TestNavigationFrame>(
        ephemeris_.get(),
        ephemeris_->bodies().back());
    EXPECT_OK(root_.Append(
        t0_ - 2 * π * Second,
        {Barycentric::origin + Displacement<Barycentric>(
                                  {1 * Metre, 0 * Metre, 0 * Metre}),
         Velocity<Barycentric>({0 * Metre / Second,
                                1 * Metre / Second,
                                0 * Metre / Second})}));
    EXPECT_OK(root_.Append(
        t0_ + 2 * π * Second,
        {Barycentric::origin + Displacement<Barycentric>(
                                  {1 * Metre, 0 * Metre, 0 * Metre}),
        Velocity<Barycentric>({0 * Metre / Second,
                                1 * Metre / Second,
                                0 * Metre / Second})}));
    flight_plan_ = std::make_unique<FlightPlan>(
        /*initial_mass=*/1 * Kilogram,
        /*initial_time=*/root_.front().time,
        /*initial_degrees_of_freedom=*/root_.front().degrees_of_freedom,
        /*desired_final_time=*/t0_ + 1.5 * Second,
        ephemeris_.get(),
        Ephemeris<Barycentric>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
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
  root_.ForgetAfter(root_.front().time);
  // NOTE(egg): In order for to avoid singular Frenet frames NaNing everything,
  // we offset our test particle by 100 ε.  The resulting system is still
  // extremely stiff, indeed the integrator detects a singularity at the exact
  // same time.  We could avoid doing this if we had absolute direction
  // specification for manœuvres.
  EXPECT_OK(root_.Append(
      t0_,
      {Barycentric::origin +
          Displacement<Barycentric>(
              {x0,
                100 * std::numeric_limits<double>::epsilon() * Metre,
                0 * Metre}),
      Barycentric::unmoving}));
  flight_plan_ = std::make_unique<FlightPlan>(
      /*initial_mass=*/1 * Kilogram,
      /*initial_time=*/root_.back().time,
      /*initial_degrees_of_freedom=*/root_.back().degrees_of_freedom,
      /*final_time=*/singularity + 100 * Second,
      ephemeris_.get(),
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*max_steps=*/1000,
          /*length_integration_tolerance=*/1 * Milli(Metre),
          /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second));
  auto const segment0 = flight_plan_->GetSegment(0);
  DiscreteTrajectory<Barycentric>::iterator back = segment0->end();
  --back;
  EXPECT_THAT(AbsoluteError(singularity, back->time), Lt(1e-4 * Second));
  // Attempting to put a burn past the singularity fails.
  EXPECT_THAT(
      flight_plan_->Insert(
          MakeTangentBurn(/*thrust=*/1 * Newton,
                          /*specific_impulse=*/1 * Newton * Second / Kilogram,
                          /*initial_time=*/singularity + 1 * Milli(Second),
                          /*Δv=*/1 * Metre / Second),
          0),
      StatusIs(termination_condition::VanishingStepSize));
  EXPECT_EQ(1, flight_plan_->number_of_anomalous_manœuvres());

  // Add another manœuvre and check the status.
  EXPECT_THAT(
      flight_plan_->Insert(
          MakeTangentBurn(/*thrust=*/1 * Newton,
                          /*specific_impulse=*/1 * Newton * Second / Kilogram,
                          /*initial_time=*/singularity + 10 * Second,
                          /*Δv=*/1 * Metre / Second),
          1),
      StatusIs(termination_condition::VanishingStepSize));
  EXPECT_EQ(2, flight_plan_->number_of_anomalous_manœuvres());

  // Check that Remove returns the proper statuses.
  EXPECT_THAT(flight_plan_->Remove(1),
              StatusIs(termination_condition::VanishingStepSize));
  EXPECT_THAT(flight_plan_->Remove(0),
              StatusIs(termination_condition::VanishingStepSize));

  // The singularity occurs during the burn: we're boosting towards the
  // singularity, so we reach the singularity in less than π / 2√2 s, before the
  // end of the burn which lasts 10 (1 - 1/e) s.
  // The derivation of analytic expression for the time at which we reach the
  // singularity is left as an exercise to the reader.
  EXPECT_THAT(
      flight_plan_->Insert(
          MakeTangentBurn(/*thrust=*/1 * Newton,
                          /*specific_impulse=*/1 * Newton * Second / Kilogram,
                          /*initial_time=*/t0_ + 0.5 * Second,
                          /*Δv=*/1 * Metre / Second),
          0),
      StatusIs(termination_condition::VanishingStepSize));
  EXPECT_EQ(0, flight_plan_->number_of_anomalous_manœuvres());

  auto segment1 = flight_plan_->GetSegment(1);
  back = segment1->end();
  --back;
  EXPECT_THAT(back->time, Lt(singularity));
  EXPECT_NE(segment1->begin(), back);
  auto segment2 = flight_plan_->GetSegment(2);
  back = segment2->end();
  --back;
  EXPECT_EQ(segment2->begin(), back);

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
      StatusIs(termination_condition::VanishingStepSize));
  EXPECT_EQ(0, flight_plan_->number_of_anomalous_manœuvres());

  segment1 = flight_plan_->GetSegment(1);
  back = segment1->end();
  --back;
  EXPECT_THAT(back->time, Eq(t0_ + 0.5 * Second + (1 - 1 / e) / 10 * Second));
  EXPECT_NE(segment1->begin(), back);
  segment2 = flight_plan_->GetSegment(2);
  back = segment2->end();
  --back;
  EXPECT_THAT(back->time, AllOf(Gt(singularity), Lt(t0_ + 2 * Second)));
  EXPECT_NE(segment2->begin(), back);
}

TEST_F(FlightPlanTest, Append) {
  // Burn ends after final time.
  EXPECT_THAT(flight_plan_->Insert(MakeFirstBurn(), 0),
              StatusIs(FlightPlan::does_not_fit));
  EXPECT_EQ(0, flight_plan_->number_of_manœuvres());
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_THAT(flight_plan_->Insert(MakeFirstBurn(), 1),
              StatusIs(FlightPlan::does_not_fit));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 1));
  EXPECT_EQ(2, flight_plan_->number_of_manœuvres());
}

TEST_F(FlightPlanTest, RemoveLast) {
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 1));
  EXPECT_EQ(2, flight_plan_->number_of_manœuvres());
  EXPECT_OK(flight_plan_->Remove(1));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_OK(flight_plan_->Remove(0));
  EXPECT_EQ(0, flight_plan_->number_of_manœuvres());
  // Check that appending still works.
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
}

TEST_F(FlightPlanTest, Replace) {
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 1.7 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  Mass const old_final_mass =
      flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
          final_mass();
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_OK(flight_plan_->Replace(MakeThirdBurn(), /*index=*/0));
  EXPECT_GT(old_final_mass,
            flight_plan_->GetManœuvre(flight_plan_->number_of_manœuvres() - 1).
                final_mass());
  EXPECT_EQ(1, flight_plan_->number_of_manœuvres());
  EXPECT_LT(t0_ + 1.7 * Second, flight_plan_->desired_final_time());
}

TEST_F(FlightPlanTest, Segments) {
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_EQ(3, flight_plan_->number_of_segments());
  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 1));
  EXPECT_EQ(5, flight_plan_->number_of_segments());

  std::vector<Instant> times;
  int last_times_size = times.size();
  Instant last_t = t0_ - 2 * π * Second;
  for (int i = 0; i < flight_plan_->number_of_segments(); ++i) {
    auto const segment_i = flight_plan_->GetSegment(i);
    for (auto const& [t, _] : *segment_i) {
      EXPECT_LE(last_t, t);
      EXPECT_LE(t, t0_ + 42 * Second);
      times.push_back(t);
    }
    EXPECT_LT(last_times_size, times.size());
    last_times_size = times.size();
  }
}

TEST_F(FlightPlanTest, SetAdaptiveStepParameter) {
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 1));
  EXPECT_EQ(5, flight_plan_->number_of_segments());
  auto segment4 = flight_plan_->GetSegment(4);
  EXPECT_EQ(t0_ + 42 * Second, segment4->back().time);

  auto const adaptive_step_parameters =
      flight_plan_->adaptive_step_parameters();
  auto const generalized_adaptive_step_parameters =
      flight_plan_->generalized_adaptive_step_parameters();

  // Reduce `max_steps`.  This causes many segments to become truncated so the
  // call to `SetAdaptiveStepParameters` returns false and the flight plan is
  // unaffected.
  EXPECT_THAT(flight_plan_->SetAdaptiveStepParameters(
        Ephemeris<Barycentric>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*max_steps=*/1,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
            /*max_steps=*/1,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)),
      StatusIs(termination_condition::ReachedMaximalStepCount));
  EXPECT_EQ(2, flight_plan_->number_of_anomalous_manœuvres());

  EXPECT_OK(flight_plan_->SetAdaptiveStepParameters(
      adaptive_step_parameters, generalized_adaptive_step_parameters));

  EXPECT_EQ(5, flight_plan_->number_of_segments());
  segment4 = flight_plan_->GetSegment(4);
  EXPECT_EQ(t0_ + 42 * Second, segment4->back().time);

  // Increase `max_steps`.  It works.
  EXPECT_OK(flight_plan_->SetAdaptiveStepParameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*max_steps=*/10000,
          /*length_integration_tolerance=*/1 * Milli(Metre),
          /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Milli(Metre),
          /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)));

  EXPECT_EQ(5, flight_plan_->number_of_segments());
  segment4 = flight_plan_->GetSegment(4);
  EXPECT_EQ(t0_ + 42 * Second, segment4->back().time);
}

TEST_F(FlightPlanTest, GuidedBurn) {
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  auto unguided_burn = MakeFirstBurn();
  unguided_burn.thrust /= 10;
  EXPECT_OK(flight_plan_->Insert(std::move(unguided_burn), 0));
  auto const& trajectory = flight_plan_->GetAllSegments();
  Speed const unguided_final_speed =
      trajectory.back().degrees_of_freedom.velocity().Norm();
  auto guided_burn = MakeFirstBurn();
  guided_burn.thrust /= 10;
  guided_burn.is_inertially_fixed = false;
  EXPECT_OK(flight_plan_->Replace(std::move(guided_burn), /*index=*/0));
  Speed const guided_final_speed =
      trajectory.back().degrees_of_freedom.velocity().Norm();
  EXPECT_THAT(guided_final_speed, IsNear(1.40_(1) * unguided_final_speed));
}

TEST_F(FlightPlanTest, Issue2331) {
  FlightPlan flight_plan(
      11024.436950683594 * Kilogram,
      J2000 + 1130.8399999993526 * Second,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin +
              Displacement<Barycentric>({-13607881359.190962 * Metre,
                                         -1183506.3832585660 * Metre,
                                         -109366.09411247486 * Metre}),
          Velocity<Barycentric>({-639.01032564318075 * Metre / Second,
                                 -8663.2671328991710 * Metre / Second,
                                 -129.00845706608661 * Metre / Second})),
      J2000 + 4280.4599499295282 * Second,
      ephemeris_.get(),
      DefaultPredictionParameters(),
      DefaultBurnParameters());

  Force const thrust = 250000.10393791363 * Newton;
  SpecificImpulse const specific_impulse = 3432.3274999999999 * Metre / Second;
  auto const frame =
      make_not_null_shared<TestNavigationFrame>(*navigation_frame_);
  bool const inertially_fixed = true;

  NavigationManœuvre::Intensity intensity0;
  intensity0.Δv =
      Velocity<Frenet<NavigationFrame>>({2035.0000000000005 * Metre / Second,
                                         0 * Metre / Second,
                                         0 * Metre / Second});
  NavigationManœuvre::Timing timing0;
  timing0.initial_time = J2000 + 3894.6399999993528 * Second;
  NavigationManœuvre::Burn const burn0{intensity0,
                                       timing0,
                                       thrust,
                                       specific_impulse,
                                       frame,
                                       inertially_fixed};
  EXPECT_OK(flight_plan.Insert(burn0, 0));

  NavigationManœuvre::Intensity intensity1;
  intensity1.Δv =
      Velocity<Frenet<NavigationFrame>>({819.29427681721018 * Metre / Second,
                                         0 * Metre / Second,
                                         0 * Metre / Second});
  NavigationManœuvre::Timing timing1;
  timing1.initial_time = J2000 + 4258.1383894665723 * Second;
  NavigationManœuvre::Burn const burn1{intensity1,
                                       timing1,
                                       thrust,
                                       specific_impulse,
                                       frame,
                                       inertially_fixed};
  EXPECT_OK(flight_plan.Insert(burn1, 1));

  NavigationManœuvre::Intensity intensity2;
  intensity2.Δv = Velocity<Frenet<NavigationFrame>>();
  NavigationManœuvre::Timing timing2;
  timing2.initial_time = J2000 + 3894.6399999993528 * Second;
  NavigationManœuvre::Burn const burn2{intensity2,
                                       timing2,
                                       thrust,
                                       specific_impulse,
                                       frame,
                                       inertially_fixed};

  // This call used to check-fail.
  EXPECT_OK(flight_plan.Replace(burn2, 0));
}

TEST_F(FlightPlanTest, Serialization) {
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 1));

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

  // This sleep causes the analyzer to continue running and produce an
  // anomalistic period of -1.42708553168389347e+01 s.
  std::this_thread::sleep_for(5s);
}

TEST_F(FlightPlanTest, Copy) {
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 1));

  serialization::FlightPlan message1;
  flight_plan_->WriteToMessage(&message1);

  FlightPlan const flight_plan_copy(*flight_plan_);
  serialization::FlightPlan message2;
  flight_plan_copy.WriteToMessage(&message2);

  EXPECT_TRUE(message2.has_initial_mass());
  EXPECT_TRUE(message2.has_initial_time());
  EXPECT_TRUE(message2.has_desired_final_time());
  EXPECT_TRUE(message2.has_adaptive_step_parameters());
  EXPECT_TRUE(message2.adaptive_step_parameters().has_integrator());
  EXPECT_TRUE(message2.adaptive_step_parameters().has_max_steps());
  EXPECT_TRUE(
      message2.adaptive_step_parameters().has_length_integration_tolerance());
  EXPECT_TRUE(
      message2.adaptive_step_parameters().has_speed_integration_tolerance());
  EXPECT_EQ(2, message2.manoeuvre_size());

  EXPECT_THAT(message2, EqualsProto(message1));
}

TEST_F(FlightPlanTest, Insertion) {
  // Check that we get the same flight plan if we add the manœuvres in the
  // opposite order.
  EXPECT_OK(flight_plan_->SetDesiredFinalTime(t0_ + 42 * Second));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 1));
  EXPECT_THAT(flight_plan_->number_of_manœuvres(), Eq(2));
  EXPECT_THAT(flight_plan_->number_of_anomalous_manœuvres(), Eq(0));

  serialization::FlightPlan inserted_in_order;
  flight_plan_->WriteToMessage(&inserted_in_order);

  // This first removal removes a non-last manœuvre.
  EXPECT_OK(flight_plan_->Remove(0));
  EXPECT_OK(flight_plan_->Remove(0));

  EXPECT_OK(flight_plan_->Insert(MakeSecondBurn(), 0));
  EXPECT_OK(flight_plan_->Insert(MakeFirstBurn(), 0));
  EXPECT_THAT(flight_plan_->number_of_manœuvres(), Eq(2));
  EXPECT_THAT(flight_plan_->number_of_anomalous_manœuvres(), Eq(0));

  serialization::FlightPlan inserted_out_of_order;
  flight_plan_->WriteToMessage(&inserted_out_of_order);

  EXPECT_THAT(inserted_out_of_order, EqualsProto(inserted_in_order));
}

}  // namespace ksp_plugin
}  // namespace principia
