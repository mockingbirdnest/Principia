
#include "integrators/explicit_embedded_runge_kutta_nyström_integrator.hpp"

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
  LOG(INFO) << (r > 1.0 ? "Accepting" : "Rejecting") << " step size "
             << h << " with ratio " << r;
  return r;
}

}  // namespace

class ExplicitEmbeddedRungeKuttaNyströmIntegratorTest
    : public ::testing::Test {};

TEST_F(ExplicitEmbeddedRungeKuttaNyströmIntegratorTest,
       HarmonicOscillatorBackAndForth) {
  ExplicitEmbeddedRungeKuttaNyströmIntegrator const& integrator =
      DormandElMikkawyPrince1986RKN434FM();
  Length const x_0 = 1 * Metre;
  Speed const v_0 = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Time const t_0 = 0 * Second;
  Time const t_final = 10 * period;

  ExplicitEmbeddedRungeKuttaNyströmIntegrator::Solution<Length, Speed> solution;
  integrator.Solve<Length>(
      ComputeHarmonicOscillatorAcceleration,
      {{x_0}, {v_0}, t_0},
      t_final,
      t_final - t_0,
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, _3, 1 * Milli(Metre), 1 * Milli(Metre) / Second),
      0.9,
      &solution);
  EXPECT_THAT(AbsoluteError(x_0, solution.back().positions[0].value),
              AllOf(Ge(1 * Milli(Metre)), Le(2 * Milli(Metre))));
  EXPECT_THAT(AbsoluteError(v_0, solution.back().momenta[0].value),
              AllOf(Ge(1 * Centi(Metre) / Second),
                    Le(2 * Centi(Metre) / Second)));
  EXPECT_EQ(t_final, solution.back().time.value);
  integrator.Solve<Length>(
      ComputeHarmonicOscillatorAcceleration,
      solution.back(),
      t_0,
      t_0 - t_final,
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, _3, 1 * Milli(Metre), 1 * Milli(Metre) / Second),
      0.9,
      &solution);
  EXPECT_THAT(AbsoluteError(x_0, solution.back().positions[0].value),
              AllOf(Ge(3 * Milli(Metre)), Le(4 * Milli(Metre))));
  EXPECT_THAT(AbsoluteError(v_0, solution.back().momenta[0].value),
              AllOf(Ge(1E-5 * Metre / Second),
                    Le(1E-4 * Centi(Metre) / Second)));
  EXPECT_EQ(t_0, solution.back().time.value);
}

}  // namespace integrators
}  // namespace principia
