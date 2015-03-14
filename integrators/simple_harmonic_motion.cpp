
#include <algorithm>
#include <string>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {

using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Energy;
using quantities::Force;
using quantities::Length;
using quantities::Mass;
using quantities::Speed;
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
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
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

#define INTEGRATOR(name) &name(), #name

namespace integrators {

namespace {

struct SimpleHarmonicMotionTestableProperties {
  not_null<SRKNIntegrator const*> integrator;
  std::string name;
  int convergence_order;
  // Convergence, in the sense tested below, of the above order occurs for
  // timesteps smaller than |beginning_of_convergence| when integrating the unit
  // harmonic oscillator.
  Time beginning_of_convergence;
  // The expected errors when integrating the unit harmonic oscillator for
  // 1000 s with a 1 ms timestep.
  Length expected_position_error;
  Speed expected_velocity_error;
  // The expected errors when integrating the unit harmonic oscillator for
  // 500 s with a 0.2 s timestep.
  Energy expected_energy_error;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream,
                         SimpleHarmonicMotionTestableProperties param) {
  return stream << param.name;
}

std::vector<SimpleHarmonicMotionTestableProperties> Instances() {
  return {
      {INTEGRATOR(Leapfrog), 2,
       0.4 * Second,
       +4.15606749774469300e-05 * Metre,
       +4.16264386218978200e-05 * Metre / Second,
       +5.05049535215751360e-03 * Joule},
      {INTEGRATOR(PseudoLeapfrog), 2,
       0.4 * Second,
       +4.15606749774360880e-05 * Metre,
       +4.16261832564750020e-05 * Metre / Second,
       +4.99999039863535670e-03 * Joule},
      {INTEGRATOR(McLachlanAtela1992Order2Optimal), 2,
       0.7 * Second,
       +2.01685999379921760e-05 * Metre,
       +2.02003819904818380e-05 * Metre / Second,
       +8.63068191495619530e-05 * Joule},
      {INTEGRATOR(Ruth1983), 3, 0.1 * Second,
       +2.77767866216707680e-11 * Metre,
       +7.01570745942348140e-13 * Metre / Second,
       +1.15535032619074050e-04 * Joule},
      {INTEGRATOR(McLachlanAtela1992Order3Optimal), 3,
       0.09 * Second,
       +1.21425465168800710e-11 * Metre,
       +3.68977418063742850e-13 * Metre / Second,
       +4.72513762963533420e-05 * Joule},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousMomenta), 4,
       0.5 * Second,
       +6.63488482904872610e-11 * Metre,
       +6.64555094981311710e-11 * Metre / Second,
       +6.26780491450595890e-05 * Joule},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousPositions), 4,
       0.5 * Second,
       +6.63488188001881700e-11 * Metre,
       +6.64553134743783860e-11 * Metre / Second,
       +6.26859072369034820e-05 * Joule},
      {INTEGRATOR(McLachlanAtela1992Order4Optimal), 4,
       1.0 * Second,
       +1.88161985992252310e-13 * Metre,
       +1.88491583452687910e-13 * Metre / Second,
       +7.52285331973023830e-07 * Joule},
      {INTEGRATOR(McLachlanAtela1992Order5Optimal), 5,
       1.1 * Second,
       +7.51005160837259210e-14 * Metre,
       +7.50823014872281650e-14 * Metre / Second,
       +3.06327349042234690e-08 * Joule},
      {INTEGRATOR(Yoshida1990Order6A), 6,
       1.5 * Second,
       +8.31001933931929670e-14 * Metre,
       +8.30759072645292920e-14 * Metre / Second,
       +1.28253664799515830e-07 * Joule},
      {INTEGRATOR(Yoshida1990Order6B), 6,
       1.0 * Second,
       +3.32536082003898060e-13 * Metre,
       +3.32810168313102390e-13 * Metre / Second,
       +3.39431978740867280e-06 * Joule},
      {INTEGRATOR(Yoshida1990Order6C), 6,
       1.0 * Second,
       +9.56665302531689580e-14 * Metre,
       +9.57515317034918210e-14 * Metre / Second,
       +3.58056353211289040e-06 * Joule},
      {INTEGRATOR(Yoshida1990Order8A), 8,
       0.043 * Second,
       +6.05702987765965870e-13 * Metre,
       +6.06313610429509710e-13 * Metre / Second,
       +1.49030436414898660e-05 * Joule},
      {INTEGRATOR(Yoshida1990Order8B), 8,
       0.5 * Second,
       +4.91471446872893130e-13 * Metre,
       +4.91932883317502960e-13 * Metre / Second,
       +1.33083072562101280e-07 * Joule},
      {INTEGRATOR(Yoshida1990Order8C), 8,
       0.9 * Second,
       +3.12770642718618320e-13 * Metre,
       +3.13027381793062890e-13 * Metre / Second,
       +4.68151011290274250e-08 * Joule},
      {INTEGRATOR(Yoshida1990Order8D), 8,
       1.1 * Second,
       +2.20323759236862320e-13 * Metre,
       +2.20518048266171720e-13 * Metre / Second,
       +1.58094926039353820e-10 * Joule},
      {INTEGRATOR(Yoshida1990Order8E), 8,
       0.3 * Second,
       +1.38892369827559040e-13 * Metre,
       +1.38979106001357880e-13 * Metre / Second,
       +3.42872149006190340e-08 * Joule}};
}

}  // namespace

class SimpleHarmonicMotionTest
    : public testing::TestWithParam<SimpleHarmonicMotionTestableProperties> {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
  }

 protected:
  SimpleHarmonicMotionTest() : integrator_(GetParam().integrator) {}

  not_null<SRKNIntegrator const*> const integrator_;
  SRKNIntegrator::Parameters<Length, Speed> parameters_;
  SRKNIntegrator::Solution<Length, Speed> solution_;
};

INSTANTIATE_TEST_CASE_P(SimpleHarmonicMotionTests, SimpleHarmonicMotionTest,
                        ValuesIn(Instances()));

TEST_P(SimpleHarmonicMotionTest, Error) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Speed());
  parameters_.initial.time = Time();
#ifdef _DEBUG
  parameters_.tmax = 100.0 * SIUnit<Time>();
#else
  parameters_.tmax = 1000.0 * SIUnit<Time>();
#endif
  parameters_.Δt = 1.0E-3 * SIUnit<Time>();
  parameters_.sampling_period = 1;
  integrator_->SolveTrivialKineticEnergyIncrement<Length>(
      ComputeHarmonicOscillatorAcceleration,
      parameters_, &solution_);
  Length q_error;
  Speed v_error;
  for (std::size_t i = 0; i < solution_.size(); ++i) {
    q_error = std::max(q_error,
                       Abs(solution_[i].positions[0].value -
                           SIUnit<Length>() *
                           Cos(solution_[i].time.value *
                               SIUnit<AngularFrequency>())));
    v_error = std::max(v_error,
                       Abs(solution_[i].momenta[0].value +
                           SIUnit<Speed>() *
                           Sin(solution_[i].time.value *
                               SIUnit<AngularFrequency>())));
  }
  LOG(INFO) << GetParam();
  LOG(INFO) << "q_error = " << q_error;
  LOG(INFO) << "v_error = " << v_error;
#ifdef _DEBUG
  EXPECT_GE(GetParam().expected_position_error, q_error);
  EXPECT_GE(GetParam().expected_velocity_error, v_error);
#else
  EXPECT_EQ(GetParam().expected_position_error, q_error);
  EXPECT_EQ(GetParam().expected_velocity_error, v_error);
#endif
  // Check consistency with the more general integration schemes.
  SPRKIntegrator const* const sprk =
      dynamic_cast<SPRKIntegrator const*>(&*integrator_);
  if (sprk != nullptr) {
    SPRKIntegrator::Parameters<Length, Momentum> parameters;
    parameters.initial.momenta.emplace_back(
        parameters_.initial.momenta.back().value * Kilogram);
    parameters.initial.positions = parameters_.initial.positions;
    parameters.initial.time      = parameters.initial.time;
    parameters.sampling_period   = parameters_.sampling_period;
    parameters.tmax              = parameters_.tmax;
    parameters.tmax_is_exact     = parameters_.tmax_is_exact;
    parameters.Δt                = parameters_.Δt;
    SPRKIntegrator::Solution<Length, Momentum> solution;
    sprk->SolveIncrement<Length, Momentum>(
        ComputeHarmonicOscillatorForce,
        ComputeHarmonicOscillatorVelocity,
        parameters, &solution);
    Length q_error;
    Speed v_error;
    for (std::size_t i = 0; i < solution.size(); ++i) {
      EXPECT_EQ(solution_[i].momenta.back().value * Kilogram,
                solution[i].momenta.back().value);
      EXPECT_EQ(solution_[i].positions.back().value,
                solution[i].positions.back().value);
      EXPECT_EQ(solution_[i].time.value, solution[i].time.value);
    }
  }
}

TEST_P(SimpleHarmonicMotionTest, Convergence) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Speed());
  parameters_.initial.time = Time();
  parameters_.tmax = 100 * SIUnit<Time>();
  parameters_.sampling_period = 0;
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
    integrator_->SolveTrivialKineticEnergyIncrement<Length>(
        &ComputeHarmonicOscillatorAcceleration,
        parameters_,
        &solution_);
    double const log_q_error = std::log10(
        std::abs(solution_[0].positions[0].value / SIUnit<Length>() -
                 Cos(solution_[0].time.value *
                     SIUnit<AngularFrequency>())));
    double const log_p_error = std::log10(
        std::abs(solution_[0].momenta[0].value / SIUnit<Speed>() +
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
#if 0
  LOG(INFO) << "Convergence data for q :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
#endif
  EXPECT_THAT(RelativeError(GetParam().convergence_order, q_convergence_order),
              Lt(0.02));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.99), Lt(1.01)));
  double const v_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const v_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << v_convergence_order;
  LOG(INFO) << "Correlation            : " << v_correlation;
#if 0
  LOG(INFO) << "Convergence data for p :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
#endif
  // SPRKs with odd convergence order have a higher convergence order in p.
  EXPECT_THAT(
     RelativeError(((GetParam().convergence_order + 1) / 2) * 2,
                   v_convergence_order),
     Lt(0.02));
  EXPECT_THAT(v_correlation, AllOf(Gt(0.99), Lt(1.01)));
}

TEST_P(SimpleHarmonicMotionTest, Symplecticity) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Speed());
  parameters_.initial.time = Time();
  Stiffness const k = SIUnit<Stiffness>();
  Mass const m      = SIUnit<Mass>();
  Length const q0   = parameters_.initial.positions[0].value;
  Speed const v0 = parameters_.initial.momenta[0].value;
  Energy const initial_energy = 0.5 * m * Pow<2>(v0) + 0.5 * k * Pow<2>(q0);
  parameters_.tmax = 500.0 * SIUnit<Time>();
  parameters_.Δt = 0.2 * Second;
  parameters_.sampling_period = 1;
  integrator_->SolveTrivialKineticEnergyIncrement<Length>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  std::size_t const length = solution_.size();
  std::vector<Energy> energy_error(length);
  std::vector<Time> time_steps(length);
  Energy max_energy_error = 0 * SIUnit<Energy>();
  for (std::size_t i = 0; i < length; ++i) {
    Length const q_i   = solution_[i].positions[0].value;
    Speed const v_i = solution_[i].momenta[0].value;
    time_steps[i] = solution_[i].time.value;
    energy_error[i] = Abs(0.5 * m * Pow<2>(v_i) + 0.5 * k * Pow<2>(q_i) -
                          initial_energy);
    max_energy_error = std::max(energy_error[i], max_energy_error);
  }
#if 0
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
