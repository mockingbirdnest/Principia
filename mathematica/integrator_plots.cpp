
#include "mathematica/integrator_plots.hpp"

#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <vector>

#include "glog/logging.h"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "mathematica/mathematica.hpp"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"

#define INTEGRATOR(name) &integrators::name(), #name

namespace principia {

using base::not_null;
using constants::GravitationalConstant;
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
using testing_utilities::ComputeKeplerAcceleration;

namespace mathematica {

namespace {

struct SimpleHarmonicMotionPlottedIntegrator {
  not_null<SRKNIntegrator const*> integrator;
  std::string name;
  int stages;
};

// This list should be sorted by:
// 1. increasing order;
// 2. SPRKs before SRKNs;
// 3. increasing number of effective stages;
// 4. date;
// 5. author names;
// 6. method name.
std::vector<SimpleHarmonicMotionPlottedIntegrator> Methods() {
  return {
      // Order 2
      {INTEGRATOR(Leapfrog), 1},
      {INTEGRATOR(PseudoLeapfrog), 1},
      {INTEGRATOR(McLachlanAtela1992Order2Optimal), 2},
      {INTEGRATOR(McLachlan1995S2), 2},
      // Order 3
      {INTEGRATOR(Ruth1983), 3},
      {INTEGRATOR(McLachlanAtela1992Order3Optimal), 3},
      // Order 4
      //   SPRKs
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousMomenta), 3},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousPositions), 3},
      {INTEGRATOR(Suzuki1990), 5},
      {INTEGRATOR(McLachlan1995SS5), 5},
      {INTEGRATOR(McLachlan1995S4), 4},
      {INTEGRATOR(McLachlan1995S5), 5},
      {INTEGRATOR(BlanesMoan2002S6), 6},
      //   SRKNs
      {INTEGRATOR(McLachlanAtela1992Order4Optimal), 4},
      {INTEGRATOR(McLachlan1995SB3A4), 4},
      {INTEGRATOR(McLachlan1995SB3A5), 5},
      {INTEGRATOR(BlanesMoan2002SRKN6B), 6},
      // Order 5
      {INTEGRATOR(McLachlanAtela1992Order5Optimal), 6},
      // Order 6
      //   SPRKs
      {INTEGRATOR(Yoshida1990Order6A), 7},
      {INTEGRATOR(Yoshida1990Order6B), 7},
      {INTEGRATOR(Yoshida1990Order6C), 7},
      {INTEGRATOR(McLachlan1995SS9), 9},
      {INTEGRATOR(BlanesMoan2002S10), 10},
      //   SRKNs
      {INTEGRATOR(OkunborSkeel1994Order6Method13), 7},
      {INTEGRATOR(BlanesMoan2002SRKN11B), 11},
      {INTEGRATOR(BlanesMoan2002SRKN14A), 14},
      // Order 8
      {INTEGRATOR(Yoshida1990Order8A), 15},
      {INTEGRATOR(Yoshida1990Order8B), 15},
      {INTEGRATOR(Yoshida1990Order8C), 15},
      {INTEGRATOR(Yoshida1990Order8D), 15},
      {INTEGRATOR(Yoshida1990Order8E), 15},
      {INTEGRATOR(McLachlan1995SS15), 15},
      {INTEGRATOR(McLachlan1995SS17), 17}};
}

}  // namespace

void GenerateSimpleHarmonicMotionWorkErrorGraphs() {
  SRKNIntegrator::Parameters<Length, Speed> parameters;
  SRKNIntegrator::Solution<Length, Speed> solution;
  std::ofstream file;
  file.open("simple_harmonic_motion_graphs.generated.wl");
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
  std::vector<std::string> q_error_data;
  std::vector<std::string> v_error_data;
  std::vector<std::string> e_error_data;
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
          method.stages *
              static_cast<int>(std::floor(parameters.tmax / parameters.Δt));
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
    q_error_data.emplace_back(PlottableDataset(evaluations, q_errors));
    v_error_data.emplace_back(PlottableDataset(evaluations, v_errors));
    e_error_data.emplace_back(PlottableDataset(evaluations, e_errors));
    names.emplace_back(Escape(method.name));
  }
  file << Assign("qErrorData", q_error_data);
  file << Assign("vErrorData", v_error_data);
  file << Assign("eErrorData", e_error_data);
  file << Assign("names", names);
  file.close();
}

void GenerateKeplerProblemWorkErrorGraphs() {
  SRKNIntegrator::Parameters<Length, Speed> parameters;
  SRKNIntegrator::Solution<Length, Speed> solution;
  std::ofstream file;
  file.open("kepler_problem_graphs.generated.wl");
  // Semi-major axis.
  Length const a = 0.5 * Metre;
  // Velocity.
  Speed const v = 0.5 * Metre / Second;
  // Gravitational parameter of the system, μ = G(m + m).
  GravitationalParameter const μ = SIUnit<GravitationalParameter>();
  Mass const m = (μ / GravitationalConstant) / 2;
  AngularFrequency const ω = 1 * Radian / Second;
  double const step_reduction = 1.015;
  // Initial conditions for the two bodies orbiting their barycentre in circular
  // orbits with semi-major axis |a|.
  parameters.initial.positions.emplace_back(2 * a);             // q_x
  parameters.initial.positions.emplace_back(0 * Metre);         // q_y
  parameters.initial.momenta.emplace_back(0 * Metre / Second);  // v_x
  parameters.initial.momenta.emplace_back(2 * v);               // v_y
  parameters.initial.time = 0 * Second;
  parameters.tmax = 50 * Second;
  // We use dense sampling in order to compute average errors, this leads to
  // more evaluations than reported for FSAL methods.
  parameters.sampling_period = 1;
  std::vector<std::string> q_error_data;
  std::vector<std::string> v_error_data;
  std::vector<std::string> e_error_data;
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
          method.stages *
              static_cast<int>(std::floor(parameters.tmax / parameters.Δt));
      LOG_IF(INFO, (i + 1) % 50 == 0) << number_of_evaluations;
      method.integrator->SolveTrivialKineticEnergyIncrement<Length>(
          &ComputeKeplerAcceleration,
          parameters,
          &solution);
      std::vector<Length> q_error;
      std::vector<Speed> v_error;
      std::vector<Energy> e_error;
      for (auto const& system_state : solution) {
        q_error.emplace_back(
            Sqrt(Pow<2>(system_state.positions[0].value -
                        2 * a * Cos(ω * system_state.time.value)) +
                 Pow<2>(system_state.positions[1].value -
                        2 * a * Sin(ω * system_state.time.value))));
        v_error.emplace_back(
            Sqrt(Pow<2>(system_state.momenta[0].value -
                        -2 * v * Sin(ω * system_state.time.value)) +
                 Pow<2>(system_state.momenta[1].value -
                        2 * v * Cos(ω * system_state.time.value))));
        Length const r_actual =
            Sqrt(Pow<2>(system_state.positions[0].value) +
                 Pow<2>(system_state.positions[1].value));
        Speed const v_actual =
            Sqrt(Pow<2>(system_state.momenta[0].value) +
                 Pow<2>(system_state.momenta[1].value)) / 2;
        e_error.emplace_back(
            AbsoluteError(
                2 * (m * v * v / 2) - GravitationalConstant * m * m / (2 * a),
                2 * (m * v_actual * v_actual / 2) -
                    GravitationalConstant * m * m / r_actual));
      }
      // We plot the maximum error, i.e., the L∞ norm of the error.
      // Blanes and Moan (2002), or Blanes, Casas and Ros (2001) tend to use
      // the average error (the normalized L¹ norm) instead.
      q_errors.emplace_back(*std::max_element(q_error.begin(), q_error.end()));
      v_errors.emplace_back(*std::max_element(v_error.begin(), v_error.end()));
      e_errors.emplace_back(*std::max_element(e_error.begin(), e_error.end()));
      evaluations.emplace_back(number_of_evaluations);
    }
    q_error_data.emplace_back(PlottableDataset(evaluations, q_errors));
    v_error_data.emplace_back(PlottableDataset(evaluations, v_errors));
    e_error_data.emplace_back(PlottableDataset(evaluations, e_errors));
    names.emplace_back(Escape(method.name));
  }
  file << Assign("qErrorData", q_error_data);
  file << Assign("vErrorData", v_error_data);
  file << Assign("eErrorData", e_error_data);
  file << Assign("names", names);
  file.close();
}

}  // namespace mathematica
}  // namespace principia
