
#define TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

using principia::quantities::Abs;
using principia::quantities::Energy;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::Momentum;
using principia::quantities::Pow;
using principia::quantities::Power;
using principia::quantities::SIUnit;
using principia::quantities::Stiffness;
using principia::quantities::Time;
using principia::testing_utilities::AbsoluteError;
using principia::testing_utilities::BidimensionalDatasetMathematicaInput;
using principia::testing_utilities::ComputeHarmonicOscillatorForce;
using principia::testing_utilities::ComputeHarmonicOscillatorVelocity;
using principia::testing_utilities::PearsonProductMomentCorrelationCoefficient;
using principia::testing_utilities::Slope;
using testing::AllOf;
using testing::Gt;
using testing::Lt;

namespace principia {
namespace integrators {

class SPRKTest : public testing::Test {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
  }

 protected:
  void SetUp() override {
    integrator_.Initialize(integrator_.Order5Optimal());
  }

  SPRKIntegrator<Length, Momentum>             integrator_;
  SPRKIntegrator<Length, Momentum>::Parameters parameters_;
  SPRKIntegrator<Length, Momentum>::Solution   solution_;
};

TEST_F(SPRKTest, HarmonicOscillator) {
  parameters_.q0 = {SIUnit<Length>()};
  parameters_.p0 = {Momentum()};
  parameters_.t0 = Time();
#ifdef _DEBUG
  parameters_.tmax = 100.0 * SIUnit<Time>();
#else
  parameters_.tmax = 1000.0 * SIUnit<Time>();
#endif
  parameters_.Δt = 1.0E-4 * SIUnit<Time>();
  parameters_.sampling_period = 1;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
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
  LOG(INFO) << "q_error = " << q_error;
  LOG(INFO) << "p_error = " << p_error;
  EXPECT_THAT(q_error, Lt(2E-16 * parameters_.tmax));
  EXPECT_THAT(p_error, Lt(2E-16 * parameters_.tmax));
}

TEST_F(SPRKTest, Convergence) {
  parameters_.q0 = {1.0};
  parameters_.p0 = {0.0};
  parameters_.t0 = 0.0;
  parameters_.tmax = 100;
  parameters_.sampling_period = 0;
  // For 0.2 * 1.1⁻²¹ < |Δt| < 0.2 , the correlation between step size and error
  // is very strong. It the step is small enough to converge and large enough to
  // stay clear of floating point inaccuracy.
  parameters_.Δt = 0.2;
  int const step_sizes = 22;
  double const step_reduction = 1.1;
  std::vector<double> log_step_sizes(step_sizes);
  std::vector<double> log_q_errors(step_sizes);
  std::vector<double> log_p_errors(step_sizes);
  for (int i = 0; i < step_sizes; ++i, parameters_.Δt /= step_reduction) {
    integrator_.Solve(&ComputeHarmonicOscillatorForce,
                      &ComputeHarmonicOscillatorVelocity,
                      parameters_, &solution_);
    log_step_sizes[i] = std::log10(parameters_.Δt);
    log_q_errors[i] = std::log10(
        std::abs(solution_.position[0].quantities[0] -
                 std::cos(solution_.time.quantities[0])));
    log_p_errors[i] = std::log10(
        std::abs(solution_.momentum[0].quantities[0] +
                 std::sin(solution_.time.quantities[0])));
  }
  double const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;
  LOG(INFO) << "Convergence data for q :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
  EXPECT_THAT(q_convergence_order, AllOf(Gt(4.9), Lt(5.1)));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.999), Lt(1.01)));
  double const p_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const p_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << p_convergence_order;
  LOG(INFO) << "Correlation            : " << p_correlation;
  LOG(INFO) << "Convergence data for p :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
  EXPECT_THAT(p_convergence_order, AllOf(Gt(5.9), Lt(6.1)));
  EXPECT_THAT(p_correlation, AllOf(Gt(0.999), Lt(1.01)));
}

TEST_F(SPRKTest, Symplecticity) {
  parameters_.q0 = {1.0};
  parameters_.p0 = {0.0};
  parameters_.t0 = 0.0;
  Stiffness const k = 1 * SIUnit<Stiffness>();
  Mass const m      = 1 * SIUnit<Mass>();
  Length const q0   = parameters_.q0[0] * SIUnit<Length>();
  Momentum const p0 = parameters_.p0[0] * SIUnit<Momentum>();
  Energy const initial_energy = 0.5 * Pow<2>(p0) / m + 0.5 * k * Pow<2>(q0);
  parameters_.tmax = 500.0;
  parameters_.Δt = 1;
  parameters_.sampling_period = 1;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  std::size_t const length = solution_.time.quantities.size();
  std::vector<Energy> energy_error(length);
  std::vector<Time> time_steps(length);
  Energy max_energy_error = 0 * SIUnit<Energy>();
  for (size_t i = 0; i < length; ++i) {
    Length const q_i   = solution_.position[0].quantities[i] * SIUnit<Length>();
    Momentum const p_i = solution_.momentum[0].quantities[i] *
                             SIUnit<Momentum>();
    time_steps[i] = solution_.time.quantities[i] * SIUnit<Time>();
    energy_error[i] = Abs(0.5 * Pow<2>(p_i) / m + 0.5 * k * Pow<2>(q_i) -
                          initial_energy);
    max_energy_error = std::max(energy_error[i], max_energy_error);
  }
  LOG(INFO) << "Energy error as a function of time:\n" <<
      BidimensionalDatasetMathematicaInput(time_steps, energy_error);
  double const correlation =
      PearsonProductMomentCorrelationCoefficient(time_steps, energy_error);
  LOG(INFO) << "Correlation between time and energy error : " << correlation;
  EXPECT_THAT(correlation, Lt(1E-3));
  Power const slope = Slope(time_steps, energy_error);
  LOG(INFO) << "Slope                                     : " << slope;
  EXPECT_THAT(slope, Lt(1E-10 * SIUnit<Power>()));
  LOG(INFO) << "Maximum energy error                      : " <<
      max_energy_error;
  EXPECT_THAT(max_energy_error, AllOf(Gt(1E-4 * SIUnit<Energy>()),
                                      Lt(1E-3 * SIUnit<Energy>())));
}

}  // namespace integrators
}  // namespace principia
