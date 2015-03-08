
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {

using quantities::Abs;
using quantities::AngularFrequency;
using quantities::Energy;
using quantities::Force;
using quantities::Length;
using quantities::Mass;
using quantities::Momentum;
using quantities::Pow;
using quantities::Power;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Stiffness;
using quantities::Time;
using si::Joule;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::BidimensionalDatasetMathematicaInput;
using testing_utilities::ComputeHarmonicOscillatorForce;
using testing_utilities::ComputeHarmonicOscillatorVelocity;
using testing_utilities::PearsonProductMomentCorrelationCoefficient;
using testing_utilities::RelativeError;
using testing_utilities::Slope;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Ne;
using ::testing::ValuesIn;

namespace integrators {

namespace {
using Integrator = SPRKIntegrator<Length, Momentum>;

struct SPRKTestableProperties {
  Integrator::Scheme const& (Integrator::*scheme)() const;
  std::string name;
  int convergence_order;
  // Convergence, in the sense tested below, of the above order occurs for
  // timesteps smaller than |beginning_of_convergence| when integrating the unit
  // harmonic oscillator.
  Time beginning_of_convergence;
  // The expected errors when integrating the unit harmonic oscillator for
  // 1000 s with a 1 ms timestep.
  Length expected_position_error;
  Momentum expected_momentum_error;
  // The expected errors when integrating the unit harmonic oscillator for
  // 500 s with a 1 s timestep.
  Energy expected_energy_error;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream, SPRKTestableProperties param) {
  return stream << param.name;
}

std::vector<SPRKTestableProperties> Instances() {
  return {
      {&Integrator::Leapfrog, "Leapfrog", 2, 0.4 * Second,
       +4.15606749774469300e-05 * Metre,
       +4.16264386218978200e-05 * Kilogram * Metre / Second,
       +1.25000000000000000e-01 * Joule},
      {&Integrator::PseudoLeapfrog, "PseudoLeapfrog", 2, 0.4 * Second,
       +4.15606749774360880e-05 * Metre,
       +4.16261832564750020e-05 * Kilogram * Metre / Second,
       +9.37500000000000000e-02 * Joule},
      {&Integrator::McLachlanAtela1992Order2Optimal,
       "McLachlanAtela1992Order2Optimal", 2, 0.7 * Second,
       +2.01685999379921760e-05 * Metre,
       +2.02003819904818380e-05 * Kilogram * Metre / Second,
       +1.28869320917870400e-02 * Joule},
      {&Integrator::Ruth1983, "Ruth1983", 3, 0.1 * Second,
       +2.77767866216707680e-11 * Metre,
       +7.01570745942348140e-13 * Kilogram * Metre / Second,
       +1.79248948993825370e-02 * Joule},
      {&Integrator::McLachlanAtela1992Order3Optimal,
       "McLachlanAtela1992Order3Optimal", 3, 0.09 * Second,
       +1.21425465168800710e-11 * Metre,
       +3.68977418063742850e-13 * Kilogram * Metre / Second,
       +6.41094324383406630e-03 * Joule},
      {&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta,
       "CandyRozmus1991ForestRuth1990SynchronousMomenta", 4, 0.5 * Second,
       +6.63488482904872610e-11 * Metre,
       +6.64555094981311710e-11 * Kilogram * Metre / Second,
       +6.98808139117831350e-02 * Joule},
      {&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions,
       "CandyRozmus1991ForestRuth1990SynchronousPositions", 4, 0.5 * Second,
       +6.63488188001881700e-11 * Metre,
       +6.64553134743783860e-11 * Kilogram * Metre / Second,
       +8.12345434555920010e-02 * Joule},
      {&Integrator::McLachlanAtela1992Order4Optimal,
       "McLachlanAtela1992Order4Optimal", 4, 1.0 * Second,
       +1.88161985992252310e-13 * Metre,
       +1.88491583452687910e-13 * Kilogram * Metre / Second,
       +5.73140348145262380e-04 * Joule},
      {&Integrator::McLachlanAtela1992Order5Optimal,
       "McLachlanAtela1992Order5Optimal", 5, 1.1 * Second,
       +7.51005160837259210e-14 * Metre,
       +7.50823014872281650e-14 * Kilogram * Metre / Second,
       +1.14708664439744370e-04 * Joule},
      {&Integrator::Yoshida1990Order6A, "Yoshida1990Order6A", 6, 1.5 * Second,
       +8.31001933931929670e-14 * Metre,
       +8.30759072645292920e-14 * Kilogram * Metre / Second,
       +2.54517718121372030e-03 * Joule},
      {&Integrator::Yoshida1990Order6B, "Yoshida1990Order6B", 6, 1.0 * Second,
       +3.32536082003898060e-13 * Metre,
       +3.32810168313102390e-13 * Kilogram * Metre / Second,
       +8.66720827531614060e-02 * Joule},
      {&Integrator::Yoshida1990Order6C, "Yoshida1990Order6C", 6, 1.0 * Second,
       +9.56665302531689580e-14 * Metre,
       +9.57515317034918210e-14 * Kilogram * Metre / Second,
       +9.43263722475094490e-02 * Joule},
      {&Integrator::Yoshida1990Order8A, "Yoshida1990Order8A", 8, 0.043 * Second,
       +8.31001933931929670e-14 * Metre,
       +8.30759072645292920e-14 * Kilogram * Metre / Second,
       +2.54517718121372030e-03 * Joule},
      {&Integrator::Yoshida1990Order8B, "Yoshida1990Order8B", 8, 0.043 * Second,
       +3.32536082003898060e-13 * Metre,
       +3.32810168313102390e-13 * Kilogram * Metre / Second,
       +8.66720827531614060e-02 * Joule},
      {&Integrator::Yoshida1990Order8C, "Yoshida1990Order8C", 8, 0.5 * Second,
       +9.56665302531689580e-14 * Metre,
       +9.57515317034918210e-14 * Kilogram * Metre / Second,
       +9.43263722475094490e-02 * Joule},
      {&Integrator::Yoshida1990Order8B, "Yoshida1990Order8D", 8, 0.043 * Second,
       +3.32536082003898060e-13 * Metre,
       +3.32810168313102390e-13 * Kilogram * Metre / Second,
       +8.66720827531614060e-02 * Joule},
      {&Integrator::Yoshida1990Order8C, "Yoshida1990Order8E", 8, 0.5 * Second,
       +9.56665302531689580e-14 * Metre,
       +9.57515317034918210e-14 * Kilogram * Metre / Second,
       +9.43263722475094490e-02 * Joule}};
}

}  // namespace

class SPRKTest : public testing::TestWithParam<SPRKTestableProperties> {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
  }

 protected:

  SPRKTest() {
    integrator_.Initialize((integrator_.*GetParam().scheme)());
  }

  Integrator                           integrator_;
  Integrator::Parameters               parameters_;
  std::vector<Integrator::SystemState> solution_;
};

INSTANTIATE_TEST_CASE_P(SPRKTests, SPRKTest, ValuesIn(Instances()));

TEST_P(SPRKTest, ConsistentWeights) {
  // Check that the time argument of the force computation is correct by
  // integrating uniform linear motion.
  // We check this for all schemes.
  Speed const v = 1 * Metre / Second;
  Mass const m = 1 * Kilogram;
  auto compute_force = [v](Time const& t,
                           std::vector<Length> const& q,
                           not_null<std::vector<Force>*> const result) {
    EXPECT_THAT(q[0], AlmostEquals(v * t, 0, 8));
    (*result)[0] = 0 * Newton;
  };
  auto compute_velocity = [m, v](std::vector<Momentum> const& p,
                                 not_null<std::vector<Speed>*> const result) {
    EXPECT_EQ(v, p[0] / m);
    (*result)[0] = p[0] / m;
  };
  parameters_.initial.positions.emplace_back(0 * Metre);
  parameters_.initial.momenta.emplace_back(m * v);
  parameters_.initial.time = 0 * Second;
  parameters_.tmax = 16 * Second;
  parameters_.Δt = 1 * Second;
  parameters_.sampling_period = 5;
  integrator_.Solve(compute_force, compute_velocity, parameters_, &solution_);
  EXPECT_THAT(v * parameters_.tmax,
              AlmostEquals(solution_.back().positions.back().value, 0, 2));
}

TEST_P(SPRKTest, HarmonicOscillator) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Momentum());
  parameters_.initial.time = Time();
#ifdef _DEBUG
  parameters_.tmax = 100.0 * SIUnit<Time>();
#else
  parameters_.tmax = 1000.0 * SIUnit<Time>();
#endif
  parameters_.Δt = 1.0E-3 * SIUnit<Time>();
  parameters_.sampling_period = 1;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  Length q_error;
  Momentum p_error;
  for (std::size_t i = 0; i < solution_.size(); ++i) {
    q_error = std::max(q_error,
                       Abs(solution_[i].positions[0].value -
                           SIUnit<Length>() *
                           Cos(solution_[i].time.value *
                               SIUnit<AngularFrequency>())));
    p_error = std::max(p_error,
                       Abs(solution_[i].momenta[0].value +
                           SIUnit<Momentum>() *
                           Sin(solution_[i].time.value *
                               SIUnit<AngularFrequency>())));
  }
  LOG(INFO) << GetParam();
  LOG(INFO) << "q_error = " << q_error;
  LOG(INFO) << "p_error = " << p_error;
  EXPECT_EQ(GetParam().expected_position_error, q_error);
  EXPECT_EQ(GetParam().expected_momentum_error, p_error);
}

TEST_P(SPRKTest, ExactInexactTMax) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Momentum());
  parameters_.initial.time = Time();
  parameters_.tmax = 10.0 * SIUnit<Time>();
  parameters_.sampling_period = 1;
  parameters_.Δt = (1.0 / 3.000001) * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(30, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Lt(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Ne(0.0 * SIUnit<Time>()));

  parameters_.tmax_is_exact = true;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(30, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));

  parameters_.Δt = (1.0 / 2.999999) * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(29, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Lt(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Ne(0.0 * SIUnit<Time>()));

  parameters_.tmax_is_exact = true;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(30, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));

  parameters_.Δt = 11.0 * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(0, solution_.size());

  parameters_.tmax_is_exact = true;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(1, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));

  parameters_.Δt = 100.0 * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(0, solution_.size());

  parameters_.tmax_is_exact = true;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  EXPECT_EQ(1, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));
}

TEST_P(SPRKTest, Convergence) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Momentum());
  parameters_.initial.time = Time();
  parameters_.tmax = 100 * SIUnit<Time>();
  parameters_.sampling_period = 0;
  // For 0.2 * 1.1⁻²¹ < |Δt| < 0.2 , the correlation between step size and error
  // is very strong. It the step is small enough to converge and large enough to
  // stay clear of floating point inaccuracy.
  parameters_.Δt = GetParam().beginning_of_convergence;
  int const step_sizes = 50;
  double const step_reduction = 1.1;
  std::vector<double> log_step_sizes;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_q_errors;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_p_errors;
  log_step_sizes.reserve(step_sizes);
  for (int i = 0; i < step_sizes; ++i, parameters_.Δt /= step_reduction) {
    integrator_.Solve(&ComputeHarmonicOscillatorForce,
                      &ComputeHarmonicOscillatorVelocity,
                      parameters_, &solution_);
    double log_q_error = std::log10(
        std::abs(solution_[0].positions[0].value / SIUnit<Length>() -
                 Cos(solution_[0].time.value *
                     SIUnit<AngularFrequency>())));
    double log_p_error = std::log10(
        std::abs(solution_[0].momenta[0].value / SIUnit<Momentum>() +
                 Sin(solution_[0].time.value *
                     SIUnit<AngularFrequency>())));
    if (log_q_error <= -13 || log_p_error <= -13) {
      // If we keep going the effects of finite precision will drown out
      // convergence.
      break;
    }
    log_step_sizes.push_back(std::log10(parameters_.Δt / SIUnit<Time>()));
    log_q_errors.push_back(log_q_error);
    log_p_errors.push_back(log_p_error);
  }
  double const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << GetParam();
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;
#if 1
  LOG(INFO) << "Convergence data for q :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
#endif
  EXPECT_THAT(RelativeError(GetParam().convergence_order, q_convergence_order),
              Lt(0.02));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.99), Lt(1.01)));
  double const p_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const p_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << p_convergence_order;
  LOG(INFO) << "Correlation            : " << p_correlation;
#if 0
  LOG(INFO) << "Convergence data for p :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
#endif
  // SPRKs with odd convergence order have a higher convergence order in p.
  EXPECT_THAT(
     RelativeError(((GetParam().convergence_order + 1) / 2) * 2,
                   p_convergence_order),
     Lt(0.02));
  EXPECT_THAT(p_correlation, AllOf(Gt(0.99), Lt(1.01)));
}

TEST_P(SPRKTest, Symplecticity) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Momentum());
  parameters_.initial.time = Time();
  Stiffness const k = SIUnit<Stiffness>();
  Mass const m      = SIUnit<Mass>();
  Length const q0   = parameters_.initial.positions[0].value;
  Momentum const p0 = parameters_.initial.momenta[0].value;
  Energy const initial_energy = 0.5 * Pow<2>(p0) / m + 0.5 * k * Pow<2>(q0);
  parameters_.tmax = 500.0 * SIUnit<Time>();
  parameters_.Δt = SIUnit<Time>();
  parameters_.sampling_period = 1;
  integrator_.Solve(&ComputeHarmonicOscillatorForce,
                    &ComputeHarmonicOscillatorVelocity,
                    parameters_, &solution_);
  std::size_t const length = solution_.size();
  std::vector<Energy> energy_error(length);
  std::vector<Time> time_steps(length);
  Energy max_energy_error = 0 * SIUnit<Energy>();
  for (std::size_t i = 0; i < length; ++i) {
    Length const q_i   = solution_[i].positions[0].value;
    Momentum const p_i = solution_[i].momenta[0].value;
    time_steps[i] = solution_[i].time.value;
    energy_error[i] = Abs(0.5 * Pow<2>(p_i) / m + 0.5 * k * Pow<2>(q_i) -
                          initial_energy);
    max_energy_error = std::max(energy_error[i], max_energy_error);
  }
#if 1
  LOG(INFO) << "Energy error as a function of time:\n" <<
      BidimensionalDatasetMathematicaInput(time_steps, energy_error);
#endif
  double const correlation =
      PearsonProductMomentCorrelationCoefficient(time_steps, energy_error);
  LOG(INFO) << GetParam();
  LOG(INFO) << "Correlation between time and energy error : " << correlation;
  EXPECT_THAT(correlation, Lt(2E-3));
  Power const slope = Slope(time_steps, energy_error);
  LOG(INFO) << "Slope                                     : " << slope;
  EXPECT_THAT(Abs(slope), Lt(2E-6 * SIUnit<Power>()));
  LOG(INFO) << "Maximum energy error                      : " <<
      max_energy_error;
  EXPECT_EQ(GetParam().expected_energy_error, max_energy_error);
}

}  // namespace integrators
}  // namespace principia
