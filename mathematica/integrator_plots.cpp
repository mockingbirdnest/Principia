
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
      {INTEGRATOR(BlanesMoan2002S6), 6},
      {INTEGRATOR(McLachlanAtela1992Order4Optimal), 4},
      {INTEGRATOR(McLachlan1995SB3A4), 4},
      {INTEGRATOR(McLachlan1995SB3A5), 5},
      {INTEGRATOR(BlanesMoan2002SRKN6B), 6},
      {INTEGRATOR(McLachlanAtela1992Order5Optimal), 6},
      {INTEGRATOR(Yoshida1990Order6A), 7},
      {INTEGRATOR(Yoshida1990Order6B), 7},
      {INTEGRATOR(Yoshida1990Order6C), 7},
      {INTEGRATOR(McLachlan1995SS9), 9},
      {INTEGRATOR(BlanesMoan2002S10), 10},
      {INTEGRATOR(OkunborSkeel1994Order6Method13), 7},
      {INTEGRATOR(BlanesMoan2002SRKN11B), 11},
      {INTEGRATOR(BlanesMoan2002SRKN14A), 14},
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
  std::vector<std::string> const axes_label = {Escape("Evaluations"),
                                               Escape(error_kind)};
  std::vector<std::string> const variables = {
      Set(data, data),
      Set(legend, legend),
      Set("visible", 
          Apply("ConstantArray", {"True", Apply("Length", {legend})}))};
  std::string const plot_style =
      Apply("Map",
            {Apply("Function",
                   {"x", Apply("Opacity", {Apply("Boole", {"x"})})}),
             "visible"});
  std::vector<std::string> const legend_label = {
      Option("True", Apply("Part", {legend, "i"})),
      Option("False", Apply("Part", {legend, "i"}))};
  std::string const plot_legend =
      Apply("Map",
            {Apply("Function",
                   {"i",
                    Apply("Toggler",
                          {Apply("Dynamic", {
                                 Apply("Part", {"visible", "i"})}),
                           ToMathematica(legend_label)})}),
             Apply("Range", {Apply("Length", {legend})})});
  return Apply(
      "DynamicModule",
      {ToMathematica(variables),
       Apply("Dynamic",
             {Apply("ListLogLogPlot",
                    {data,
                     Option("PlotStyle", plot_style),
                     Option("PlotLegends",
                            Apply("SwatchLegend", {plot_legend})),
                     Option("PlotRange", ToMathematica(range)),
                     Option("AxesLabel", ToMathematica(axes_label)),
                     Option("ImageSize", "1200")})})});
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
  // We use dense sampling in order to compute average errors, this leads to
  // more evaluations than reported for FSAL methods.
  parameters.sampling_period = 1;
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
      std::vector<Length> q_error;
      std::vector<Speed> v_error;
      std::vector<Energy> e_error;
      for (auto const& system_state : solution) {
        q_error.emplace_back(
            AbsoluteError(q_amplitude * Cos(ω * system_state.time.value),
                          system_state.positions[0].value));
        v_error.emplace_back(
            AbsoluteError(-v_amplitude * Sin(ω * system_state.time.value),
                          system_state.momenta[0].value));
        e_error.emplace_back(
            AbsoluteError(
                0.5 * Joule,
                (m * Pow<2>(system_state.momenta[0].value) +
                 k * Pow<2>(system_state.positions[0].value)) / 2));
      }
      // We plot the maximum error, i.e., the L∞ norm of the error.
      // Blanes and Moan (2002), or Blanes, Casas and Ros (2001) tend to use
      // the average error (the normalized L¹ norm) instead.
      q_errors.emplace_back(*std::max_element(q_error.begin(), q_error.end()));
      v_errors.emplace_back(*std::max_element(v_error.begin(), v_error.end()));
      e_errors.emplace_back(*std::max_element(e_error.begin(), e_error.end()));
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
  file << Export("shm_q_error.cdf",
                 ErrorPlot("qErrorPlots", "names", "Position error (m)"));
  file << Export("shm_v_error.cdf",
                 ErrorPlot("vErrorPlots", "names", "Velocity error (m/s)"));
  file << Export("shm_e_error.cdf",
                 ErrorPlot("eErrorPlots", "names", "Energy error (J)"));
  file.close();
}

}  // namespace mathematica
}  // namespace principia
