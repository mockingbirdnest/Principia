
#include "mathematica/integrator_plots.hpp"

#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <vector>

#include "glog/logging.h"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "mathematica/mathematica.hpp"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"

#define INTEGRATOR(name) &integrators::name(), #name

namespace principia {

using base::not_null;
using integrators::SRKNIntegrator;
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

namespace mathematica {

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
      {INTEGRATOR(McLachlan1995S2), 2},
      {INTEGRATOR(Ruth1983), 3},
      {INTEGRATOR(McLachlanAtela1992Order3Optimal), 3},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousMomenta), 3},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousPositions), 3},
      {INTEGRATOR(Suzuki1990), 5},
      {INTEGRATOR(McLachlan1995SS5), 5},
      {INTEGRATOR(McLachlan1995S4), 4},
      {INTEGRATOR(McLachlan1995S5), 5},
      {INTEGRATOR(McLachlanAtela1992Order4Optimal), 4},
      {INTEGRATOR(McLachlan1995SB3A4), 4},
      {INTEGRATOR(McLachlan1995SB3A5), 5},
      {INTEGRATOR(McLachlanAtela1992Order5Optimal), 6},
      {INTEGRATOR(Yoshida1990Order6A), 7},
      {INTEGRATOR(Yoshida1990Order6B), 7},
      {INTEGRATOR(Yoshida1990Order6C), 7},
      {INTEGRATOR(McLachlan1995SS9), 9},
      {INTEGRATOR(OkunborSkeel1994Order6Method13), 7},
      {INTEGRATOR(Yoshida1990Order8A), 15},
      {INTEGRATOR(Yoshida1990Order8B), 15},
      {INTEGRATOR(Yoshida1990Order8C), 15},
      {INTEGRATOR(Yoshida1990Order8D), 15},
      {INTEGRATOR(Yoshida1990Order8E), 15},
      {INTEGRATOR(McLachlan1995SS15), 15},
      {INTEGRATOR(McLachlan1995SS17), 17}};
}

std::string ErrorPlot(std::string const& data, std::string const& legend,
                      std::string const& error_kind) {
  std::vector<double> const range = {1e-16, 1};
  std::vector<std::string> const axes_label = {"Evaluations",
                                               Escape(error_kind)};
  return Apply(
      "ListLogLogPlot",
      {data,
       Option("PlotLegends", Apply("SwatchLegend", {legend})),
       Option("PlotRange", ToMathematica(range)),
       Option("AxesLabel", ToMathematica(axes_label)),
       Option("ImageSize", "1200")});
}

}  // namespace

void GenerateSimpleHarmonicMotionWorkErrorGraphs() {
  SRKNIntegrator::Parameters<Length, Speed> parameters;
  SRKNIntegrator::Solution<Length, Speed> solution;
  std::ofstream file;
  file.open("simple_harmonic_motion_graphs.m");
  Length const q_amplitude = 1 * Metre;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Stiffness const k = SIUnit<Stiffness>();
  Mass const m = 1 * Kilogram;
  double const step_reduction = 1.015;
  parameters.initial.positions.emplace_back(q_amplitude);
  parameters.initial.momenta.emplace_back(0 * Metre / Second);
  parameters.initial.time = 0 * Second;
  parameters.tmax = 50 * Second;
  parameters.sampling_period = 0;
  std::vector<std::string> q_plots;
  std::vector<std::string> v_plots;
  std::vector<std::string> e_plots;
  std::vector<std::string> names;
  for (auto const& method : Methods()) {
    LOG(INFO) << method.name;
    parameters.Δt = method.stages * 1 * Second;
    std::vector<Length> q_errors;
    std::vector<Speed> v_errors;
    std::vector<Energy> e_errors;
    std::vector<double> evaluations;
    for (int i = 0; i < 500; ++i, parameters.Δt /= step_reduction) {
      int const number_of_evaluations =
          method.stages * std::ceil(parameters.tmax / parameters.Δt);
      LOG_IF(INFO, (i + 1) % 50 == 0) << number_of_evaluations;
      method.integrator->SolveTrivialKineticEnergyIncrement<Length>(
          &ComputeHarmonicOscillatorAcceleration,
          parameters,
          &solution);
      q_errors.emplace_back(
          AbsoluteError(q_amplitude * Cos(ω * solution[0].time.value),
                        solution[0].positions[0].value));
      v_errors.emplace_back(
          AbsoluteError(-v_amplitude * Sin(ω * solution[0].time.value),
                        solution[0].momenta[0].value));
      e_errors.emplace_back(
          AbsoluteError(
              0.5 * Joule,
              (m * Pow<2>(solution[0].momenta[0].value) +
               k * Pow<2>(solution[0].positions[0].value)) / 2));
      evaluations.emplace_back(number_of_evaluations);
    }
    q_plots.emplace_back(PlottableDataset(evaluations, q_errors));
    v_plots.emplace_back(PlottableDataset(evaluations, v_errors));
    e_plots.emplace_back(PlottableDataset(evaluations, e_errors));
    names.emplace_back(Escape(method.name));
  }
  file << Assign("qErrorPlots", q_plots);
  file << Assign("vErrorPlots", v_plots);
  file << Assign("eErrorPlots", e_plots);
  file << Assign("names", names);
  file << Export("shm_q_error.png",
                 ErrorPlot("qErrorPlots", "names", "Position error (m)"));
  file << Export("shm_v_error.png",
                 ErrorPlot("vErrorPlots", "names", "Velocity error (m/s)"));
  file << Export("shm_e_error.png",
                 ErrorPlot("eErrorPlots", "names", "Energy error (J)"));
  file.close();
}

}  // namespace mathematica
}  // namespace principia
