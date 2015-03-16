
#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <vector>

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/mathematica.hpp"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"

#define INTEGRATOR(name) &name(), #name

namespace principia {

using base::not_null;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Energy;
using quantities::Length;
using quantities::Pow;
using quantities::Sin;
using quantities::Speed;
using si::Joule;
using si::Kilogram;
using si::Metre;
using si::Radian;
using si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
using testing_utilities::MathematicaAssign;
using testing_utilities::MathematicaEscape;
using testing_utilities::MathematicaFunction;
using testing_utilities::MathematicaPlottableDataset;
using testing_utilities::ToMathematica;

namespace integrators {

namespace {

struct SimpleHarmonicMotionPlottedIntegrator {
  not_null<SRKNIntegrator const*> integrator;
  std::string name;
  int stages;
};

std::vector<SimpleHarmonicMotionPlottedIntegrator> Methods() {
  return {
      {INTEGRATOR(Leapfrog), 1},
      {INTEGRATOR(PseudoLeapfrog), 1},
      {INTEGRATOR(McLachlanAtela1992Order2Optimal), 2},
      {INTEGRATOR(Ruth1983), 3},
      {INTEGRATOR(McLachlanAtela1992Order3Optimal), 3},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousMomenta), 3},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousPositions), 3},
      {INTEGRATOR(McLachlanAtela1992Order4Optimal), 4},
      {INTEGRATOR(McLachlanAtela1992Order5Optimal), 6},
      {INTEGRATOR(Yoshida1990Order6A), 7},
      {INTEGRATOR(Yoshida1990Order6B), 7},
      {INTEGRATOR(Yoshida1990Order6C), 7},
      {INTEGRATOR(Yoshida1990Order8A), 15},
      {INTEGRATOR(Yoshida1990Order8B), 15},
      {INTEGRATOR(Yoshida1990Order8C), 15},
      {INTEGRATOR(Yoshida1990Order8D), 15},
      {INTEGRATOR(Yoshida1990Order8E), 15}};
}

}  // namespace

class MathematicaPlots : public testing::Test {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
  }

 protected:
  SRKNIntegrator::Parameters<Length, Speed> parameters_;
  SRKNIntegrator::Solution<Length, Speed> solution_;
};

TEST_F(MathematicaPlots, SimpleHarmonicMotionWorkErrorGraphs) {
  std::ofstream file;
  file.open("simple_harmonic_motion_graphs.m");
  Length const q_amplitude = 1 * Metre;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Stiffness const k = SIUnit<Stiffness>();
  Mass const m = 1 * Kilogram;
  double const step_reduction = 1.014;
  parameters_.initial.positions.emplace_back(q_amplitude);
  parameters_.initial.momenta.emplace_back(0 * Metre / Second);
  parameters_.initial.time = 0 * Second;
  parameters_.tmax = 50 * Second;
  parameters_.sampling_period = 0;
  std::vector<std::string> q_plots;
  std::vector<std::string> v_plots;
  std::vector<std::string> e_plots;
  std::vector<std::string> names;
  for (auto const& method : Methods()) {
    LOG(INFO) << method.name;
    parameters_.Δt = method.stages * 0.5 * Second;
    std::vector<Length> q_errors;
    std::vector<Speed> v_errors;
    std::vector<Energy> e_errors;
    std::vector<double> evaluations;
    for (int i = 0; i < 500; ++i, parameters_.Δt /= step_reduction) {
      LOG_EVERY_N(INFO, 50) << parameters_.Δt;
      method.integrator->SolveTrivialKineticEnergyIncrement<Length>(
          &ComputeHarmonicOscillatorAcceleration,
          parameters_,
          &solution_);
      q_errors.emplace_back(
          AbsoluteError(q_amplitude * Cos(ω * solution_[0].time.value),
                        solution_[0].positions[0].value));
      v_errors.emplace_back(
          AbsoluteError(-v_amplitude * Sin(ω * solution_[0].time.value),
                        solution_[0].momenta[0].value));
      e_errors.emplace_back(
          AbsoluteError(
              0.5 * Joule,
              (m * Pow<2>(solution_[0].momenta[0].value) +
               k * Pow<2>(solution_[0].positions[0].value)) / 2));
      evaluations.emplace_back(
          method.stages * parameters_.tmax / parameters_.Δt);
    }
    q_plots.emplace_back(MathematicaPlottableDataset(evaluations, q_errors));
    v_plots.emplace_back(MathematicaPlottableDataset(evaluations, v_errors));
    e_plots.emplace_back(MathematicaPlottableDataset(evaluations, e_errors));
    names.emplace_back(MathematicaEscape(method.name));
  }
  file << MathematicaAssign("qErrorPlots", q_plots);
  file << MathematicaAssign("vErrorPlots", v_plots);
  file << MathematicaAssign("eErrorPlots", e_plots);
  file << MathematicaAssign("names", names);
  file.close();
}

}  // namespace integrators
}  // namespace principia
