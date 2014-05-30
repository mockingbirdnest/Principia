
#define TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/dimensionless.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

using principia::testing_utilities::AbsoluteError;
using testing::AllOf;
using testing::Gt;
using testing::Lt;
using principia::testing_utilities::Slope;
using principia::testing_utilities::PearsonProductMomentCorrelationCoefficient;

namespace principia {
namespace integrators {

namespace {

using quantities::Dimensionless;

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
  void SetUp() override {}

  SPRKIntegrator             integrator_;
  SPRKIntegrator::Parameters parameters_;
  SPRKIntegrator::Solution   solution_;
};

TEST_F(SPRKTest, HarmonicOscillator) {
  parameters_.q0 = {1.0};
  parameters_.p0 = {0.0};
  parameters_.t0 = 0.0;
#ifdef _DEBUG
  parameters_.tmax = 100.0;
#else
  parameters_.tmax = 1000.0;
#endif
  parameters_.Δt = 1.0E-4;
  parameters_.coefficients = integrator_.Order5Optimal();
  parameters_.sampling_period = 1;
  integrator_.Solve(&compute_harmonic_oscillator_force,
                    &compute_harmonice_oscillator_velocity,
                    parameters_, &solution_);
  double q_error = 0;
  double p_error = 0;
  for (size_t i = 0; i < solution_.time.quantities.size(); ++i) {
    q_error = std::max(q_error,
                       std::abs(solution_.position[0].quantities[i] -
                                std::cos(solution_.time.quantities[i])));
    p_error = std::max(p_error,
                       std::abs(solution_.momentum[0].quantities[i] +
                                std::sin(solution_.time.quantities[i])));
  }
  LOG(ERROR) << "q_error = " << q_error;
  LOG(ERROR) << "p_error = " << p_error;
  EXPECT_THAT(AbsoluteError(0, q_error), Lt(2E-16 * parameters_.tmax));
  EXPECT_THAT(AbsoluteError(0, p_error), Lt(2E-16 * parameters_.tmax));
}

TEST_F(SPRKTest, Convergence) {
  parameters_.q0 = {1.0};
  parameters_.p0 = {0.0};
  parameters_.t0 = 0.0;
  parameters_.tmax = 100;
  parameters_.coefficients = integrator_.Order5Optimal();
  parameters_.sampling_period = 0;
  int const step_sizes = 50;
  std::vector<Dimensionless> log_step_sizes(step_sizes);
  std::vector<Dimensionless> log_q_errors(step_sizes);
  std::vector<Dimensionless> log_p_errors(step_sizes);
  parameters_.Δt = 3;
  for(int i = 0; i < step_sizes; ++i, parameters_.Δt /= 1.1) {
    solution_ = SPRKIntegrator::Solution();
    integrator_.Solve(&compute_harmonic_oscillator_force,
                      &compute_harmonice_oscillator_velocity,
                      parameters_, &solution_);
    log_step_sizes[i] = std::log(parameters_.Δt);
    log_q_errors[i] = std::log(
        std::abs(solution_.position[0].quantities[0] -
                 std::cos(solution_.time.quantities[0])));
    log_p_errors[i] = std::log(
        std::abs(solution_.momentum[0].quantities[0] +
                 std::sin(solution_.time.quantities[0])));
  }
  std::string mathematica_out = "{";
  for (int i = 0; i < step_sizes; ++i) {
    mathematica_out += "{";
    mathematica_out += ToString(log_step_sizes[i]);
    mathematica_out += ",";
    mathematica_out += ToString(log_q_errors[i]);
    mathematica_out += "}";
    if (i + 1 < step_sizes) {
      mathematica_out += ",\n ";
    }
  }
  mathematica_out += "}";
  FLAGS_log_prefix = false;
  google::LogToStderr();
  LOG(INFO) << "Convergence data for q:\n" << mathematica_out;
  Dimensionless const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  Dimensionless const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;
  EXPECT_THAT(q_convergence_order, AllOf(Gt(4.5), Lt(6)));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.9), Lt(1.1)));
  mathematica_out = "{";
  for (int i = 0; i < step_sizes; ++i) {
    mathematica_out += "{";
    mathematica_out += ToString(log_step_sizes[i]);
    mathematica_out += ",";
    mathematica_out += ToString(log_p_errors[i]);
    mathematica_out += "}";
    if (i + 1 < step_sizes) {
      mathematica_out += ",\n ";
    }
  }
  mathematica_out += "}";
  LOG(INFO) << "Convergence data for p:\n" << mathematica_out;
  Dimensionless const p_convergence_order = Slope(log_step_sizes, log_p_errors);
  Dimensionless const p_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << p_convergence_order;
  LOG(INFO) << "Correlation            : " << p_correlation;
  EXPECT_THAT(p_convergence_order, AllOf(Gt(4.5), Lt(6)));
  EXPECT_THAT(p_correlation, AllOf(Gt(0.9), Lt(1.1)));
}

}  // namespace integrators
}  // namespace principia
