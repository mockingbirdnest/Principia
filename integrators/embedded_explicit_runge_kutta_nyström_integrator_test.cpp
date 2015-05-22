
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"

#include "base/macros.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using quantities::Abs;
using quantities::Length;
using si::Centi;
using si::Metre;
using si::Milli;
using si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;

namespace integrators {

namespace {

double HarmonicOscillatorToleranceRatio(
    Time const& h,
    std::vector<Length> const& q_error_estimate,
    std::vector<Speed> const& v_error_estimate,
    Length const& q_tolerance,
    Speed const& v_tolerance) {
  double const r = std::min(q_tolerance / Abs(q_error_estimate[0]),
                            v_tolerance / Abs(v_error_estimate[0]));
  return r;
}

}  // namespace

class EmbeddedExplicitRungeKuttaNyströmIntegratorTest
    : public ::testing::Test {};

TEST_F(EmbeddedExplicitRungeKuttaNyströmIntegratorTest,
       HarmonicOscillatorBackAndForth) {
  EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator =
      DormandElMikkawyPrince1986RKN434FM();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Time const t_initial = 0 * Second;
  Time const t_final = 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;
  int const steps_forward = 132;
  // We integrate backward with double the tolerance.
  int const steps_backward = 112;

  EmbeddedExplicitRungeKuttaNyströmIntegrator::Solution<Length, Speed> solution;
  integrator.Solve<Length>(
      ComputeHarmonicOscillatorAcceleration,
      {{x_initial}, {v_initial}, t_initial},
      t_final,
      t_final - t_initial,
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, _3, length_tolerance, speed_tolerance),
      0.9,
      &solution);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(3E-4 * Metre), Le(4E-4 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().momenta[0].value),
              AllOf(Ge(2E-3 * Metre / Second), Le(3E-3 * Metre / Second)));
  EXPECT_EQ(t_final, solution.back().time.value);
  EXPECT_EQ(steps_forward, solution.size());
  integrator.Solve<Length>(
      ComputeHarmonicOscillatorAcceleration,
      solution.back(),
      t_initial,
      t_initial - t_final,
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, _3, 2 * length_tolerance, 2 * speed_tolerance),
      0.9,
      &solution);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(1E-3 * Metre), Le(2E-3 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().momenta[0].value),
              AllOf(Ge(2E-3 * Metre / Second), Le(3E-3 * Metre / Second)));
  EXPECT_EQ(t_initial, solution.back().time.value);
  EXPECT_EQ(steps_backward, solution.size() - steps_forward);
}

}  // namespace integrators
}  // namespace principia
