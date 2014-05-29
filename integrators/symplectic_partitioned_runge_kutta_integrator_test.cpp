
#define GLOG_NO_ABBREVIATED_SEVERITIES
#define TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/numerics.hpp"

using principia::testing_utilities::AbsoluteError;
using testing::Lt;

namespace principia {
namespace integrators {

namespace {

inline void compute_harmonic_oscillator_force(double const t,
                                              std::vector<double> const& q,
                                              std::vector<double>* result) {
  (*result)[0] = -q[0];
}

inline void compute_harmonice_oscillator_velocity(std::vector<double> const& p,
                                                  std::vector<double>* result) {
  (*result)[0] = p[0];
}

}  // namespace

class SPRKTest : public testing::Test {
 protected:
  void SetUp() override {
    integrator_.reset(new SPRKIntegrator);
    parameters_.reset(new SPRKIntegrator::Parameters);
    solution_.reset(new SPRKIntegrator::Solution);
  }

  std::unique_ptr<SPRKIntegrator> integrator_;
  std::unique_ptr<SPRKIntegrator::Parameters> parameters_;
  std::unique_ptr<SPRKIntegrator::Solution> solution_;
};

TEST_F(SPRKTest, HarmonicOscillator) {
  parameters_->q0 = {1.0};
  parameters_->p0 = {0.0};
  parameters_->t0 = 0.0;
#ifdef _DEBUG
  parameters_->tmax = 100.0;
#else
  parameters_->tmax = 1000.0;
#endif
  parameters_->Δt = 1.0E-4;
  parameters_->coefficients = integrator_->Order5Optimal();
  parameters_->sampling_period = 1;
  integrator_->Solve(&compute_harmonic_oscillator_force,
                         &compute_harmonice_oscillator_velocity,
                         *parameters_,
                         solution_.get());
  double q_error = 0;
  double p_error = 0;
  for (size_t i = 0; i < solution_->time.quantities.size(); ++i) {
    q_error = std::max(q_error,
                       std::abs(solution_->position[0].quantities[i] -
                                std::cos(solution_->time.quantities[i])));
    p_error = std::max(p_error,
                       std::abs(solution_->momentum[0].quantities[i] +
                                std::sin(solution_->time.quantities[i])));
  }
  LOG(ERROR) << "q_error = " << q_error;
  LOG(ERROR) << "p_error = " << p_error;
  EXPECT_THAT(AbsoluteError(0, q_error), Lt(2E-16 * parameters_->tmax));
  EXPECT_THAT(AbsoluteError(0, p_error), Lt(2E-16 * parameters_->tmax));
}

}  // namespace integrators
}  // namespace principia
