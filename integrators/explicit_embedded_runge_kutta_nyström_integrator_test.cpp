
#include "integrators/explicit_embedded_runge_kutta_nyström_integrator.hpp"

#include "base/macros.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/integration.hpp"

namespace principia {

using quantities::Abs;
using quantities::Length;
using si::Metre;
using si::Milli;
using si::Second;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

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
  LOG(ERROR) << (r > 1.0 ? "Accepting" : "REJECTING") << " step size "
             << h << " with ratio " << r;
  return r;
}

}  // namespace

class ExplicitEmbeddedRungeKuttaNyströmIntegratorTest
    : public ::testing::Test {};

TEST_F(ExplicitEmbeddedRungeKuttaNyströmIntegratorTest, Dummy) {
  ExplicitEmbeddedRungeKuttaNyströmIntegrator const& integrator =
      DormandElMikkawyPrince1986RKN434FM();

  ExplicitEmbeddedRungeKuttaNyströmIntegrator::Solution<Length, Speed> solution;
  integrator.Solve<Length>(
      ComputeHarmonicOscillatorAcceleration,
      {{1 * Metre}, {0 * Metre / Second}, 0 * Second},
      100 * Second,
      3 * Second,
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, _3, 1 * Milli(Metre), 1 * Milli(Metre) / Second),
      0.9,
      &solution);
}

}  // namespace integrators
}  // namespace principia
