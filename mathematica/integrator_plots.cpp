
#include "mathematica/integrator_plots.hpp"

#include <algorithm>
#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "glog/logging.h"
#include "integrators/sprk_integrator.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "mathematica/mathematica.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

#define INTEGRATOR(name) &integrators::name(), #name

namespace principia {

using astronomy::JulianYear;
using base::not_null;
using constants::GravitationalConstant;
using geometry::InnerProduct;
using geometry::BarycentreCalculator;
using integrators::SRKNIntegrator;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Energy;
using quantities::Length;
using quantities::Pow;
using quantities::Sin;
using quantities::Speed;
using si::Minute;
using si::Joule;
using si::Kilogram;
using si::Metre;
using si::Radian;
using si::Second;
using testing_utilities::ICRFJ2000Ecliptic;
using testing_utilities::AbsoluteError;
using testing_utilities::ComputeGravitationalAcceleration;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
using testing_utilities::ComputeKeplerAcceleration;
using testing_utilities::SolarSystem;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

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

// Those methods which have converged to the limits of double-precision floating
// point error on the circular Kepler problem tested by
// |GenerateKeplerProblemWorkErrorGraphs| with less than 8e4 evaluations.
std::vector<SimpleHarmonicMotionPlottedIntegrator> ReferenceMethods() {
  return {
      // Order 5
      {INTEGRATOR(McLachlanAtela1992Order5Optimal), 6},
      // Order 6
      //   SPRKs
      {INTEGRATOR(McLachlan1995SS9), 9},
      {INTEGRATOR(BlanesMoan2002S10), 10},
      //   SRKNs
      {INTEGRATOR(OkunborSkeel1994Order6Method13), 7},
      {INTEGRATOR(BlanesMoan2002SRKN11B), 11},
      {INTEGRATOR(BlanesMoan2002SRKN14A), 14},
      // Order 8
      {INTEGRATOR(McLachlan1995SS15), 15},
      {INTEGRATOR(McLachlan1995SS17), 17}};
}

}  // namespace

void GenerateSimpleHarmonicMotionWorkErrorGraphs() {
  SRKNIntegrator::Parameters<Length, Speed> parameters;
  SRKNIntegrator::Solution<Length, Speed> solution;
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
  std::ofstream file;
  file.open("simple_harmonic_motion_graphs.generated.wl");
  file << Assign("qErrorData", q_error_data);
  file << Assign("vErrorData", v_error_data);
  file << Assign("eErrorData", e_error_data);
  file << Assign("names", names);
  file.close();
}

void GenerateKeplerProblemWorkErrorGraphs() {
  SRKNIntegrator::Parameters<Length, Speed> parameters;
  SRKNIntegrator::Solution<Length, Speed> solution;
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
  std::ofstream file;
  file.open("kepler_problem_graphs.generated.wl");
  file << Assign("qErrorData", q_error_data);
  file << Assign("vErrorData", v_error_data);
  file << Assign("eErrorData", e_error_data);
  file << Assign("names", names);
  file.close();
}

void GenerateSolarSystemPlanetsWorkErrorGraph() {
  SRKNIntegrator::Parameters<Position<ICRFJ2000Ecliptic>,
                             Velocity<ICRFJ2000Ecliptic>> parameters;
  Energy initial_energy;
  std::vector<MassiveBody> bodies;
  int const last_planet = SolarSystem::kMercury;
  Time const Δt_reference = 10 * Minute;
  {
    not_null<std::unique_ptr<SolarSystem>> const solar_system =
        SolarSystem::AtСпутник1Launch(SolarSystem::Accuracy::kMajorBodiesOnly);
    SolarSystem::Bodies const solar_system_bodies =
        solar_system->massive_bodies();
    for (int i = 0; i <= last_planet; ++i) {
      DiscreteTrajectory<ICRFJ2000Ecliptic> const& trajectory =
          *solar_system->trajectories()[i];
      bodies.emplace_back(*solar_system_bodies[i]);
      MassiveBody const& body = bodies.back();
      parameters.initial.positions.emplace_back(
          trajectory.last().degrees_of_freedom().position());
      Velocity<ICRFJ2000Ecliptic> const& v =
          trajectory.last().degrees_of_freedom().velocity();
      parameters.initial.momenta.emplace_back(v);
      // Kinetic energy.
      initial_energy += 0.5 * body.mass() * InnerProduct(v, v);
    }

    for (int i = 1; i <= last_planet; ++i) {
      for (int j = 0; j < i; ++j) {
        // Potential energy.
        initial_energy -=
            GravitationalConstant * bodies[i].mass() * bodies[j].mass() /
                (parameters.initial.positions[i].value -
                 parameters.initial.positions[j].value).Norm();
      }
    }
  }
  parameters.initial.time = 0 * Second;
  parameters.tmax = 1 * JulianYear;
  // We use dense sampling in order to compute average errors, this leads to
  // more evaluations than reported for FSAL methods.
  parameters.sampling_period = 1;

  SRKNIntegrator::Solution<Position<ICRFJ2000Ecliptic>,
                           Velocity<ICRFJ2000Ecliptic>> reference_solution;
  {
    std::vector<
        SRKNIntegrator::Solution<
            Position<ICRFJ2000Ecliptic>,
            Velocity<ICRFJ2000Ecliptic>>> reference_solutions;
    LOG(INFO) << "Computing reference solutions";
    for (auto const& method : ReferenceMethods()) {
      LOG(INFO) << method.name;
      parameters.Δt = Δt_reference;
      reference_solutions.emplace_back();
      method.integrator->
          SolveTrivialKineticEnergyIncrement<Position<ICRFJ2000Ecliptic>>(
              std::bind(ComputeGravitationalAcceleration<ICRFJ2000Ecliptic>,
                        _1, _2, _3, std::cref(bodies)),
              parameters,
              &reference_solutions.back());
    }
    int const reference_size = reference_solutions.front().size();
    for (auto const& solution : reference_solutions) {
      CHECK_EQ(solution.size(), reference_size);
    }
    for (int i = 0; i < reference_size; ++i) {
      reference_solution.emplace_back();
    }
    for (int b = 0; b <= last_planet; ++b) {
      for (int i = 0; i < reference_size; ++i) {
        Position<ICRFJ2000Ecliptic>::BarycentreCalculator<double> reference_q;
        BarycentreCalculator<Velocity<ICRFJ2000Ecliptic>, double> reference_v;
        for (auto const& solution : reference_solutions) {
          reference_q.Add(solution[i].positions[b].value, 1);
          reference_v.Add(solution[i].momenta[b].value, 1);
        }
        reference_solution[i].positions.emplace_back(reference_q.Get());
        reference_solution[i].momenta.emplace_back(reference_v.Get());
      }
    }
    LOG(INFO) << "Done";
  }
  SRKNIntegrator::Solution<Position<ICRFJ2000Ecliptic>,
                           Velocity<ICRFJ2000Ecliptic>> solution;
  std::vector<std::string> q_error_data;
  std::vector<std::string> v_error_data;
  std::vector<std::string> e_error_data;
  std::vector<std::string> names;
  double const step_reduction = 1.015;
  for (auto const& method : Methods()) {
    LOG(INFO) << method.name;
    parameters.Δt = method.stages * 30 * Day;
    std::vector<Length> q_errors;
    std::vector<Speed> v_errors;
    std::vector<Energy> e_errors;
    std::vector<double> evaluations;
    for (int i = 0;
         parameters.Δt > 0 * Second;
         ++i,
         parameters.Δt = Δt_reference *
             std::floor((parameters.Δt / step_reduction) / Δt_reference)) {
      int const number_of_evaluations =
          method.stages *
              static_cast<int>(std::floor(parameters.tmax / parameters.Δt));
      LOG_IF(INFO, (i + 1) % 50 == 0) << number_of_evaluations;
      method.integrator->
          SolveTrivialKineticEnergyIncrement<Position<ICRFJ2000Ecliptic>>(
              std::bind(ComputeGravitationalAcceleration<ICRFJ2000Ecliptic>,
                        _1, _2, _3, std::cref(bodies)),
              parameters,
              &solution);
      Length q_error;
      Speed v_error;
      Energy e_error;
      for (auto const& system_state : solution) {
        Energy energy;
        for (int body = 0; body <= last_planet; ++body) {
          int t = static_cast<int>(
                      std::round(system_state.time.value / Δt_reference)) - 1;
          q_error =
              std::max(q_error,
                       (system_state.positions[body].value -
                        reference_solution[t].positions[body].value).Norm());
          v_error =
              std::max(v_error,
                       (system_state.momenta[body].value -
                        reference_solution[t].momenta[body].value).Norm());
          // Kinetic energy.
          energy += 0.5 * bodies[body].mass() *
                        InnerProduct(system_state.momenta[body].value,
                                     system_state.momenta[body].value);
        }
        for (int i = 1; i <= last_planet; ++i) {
          for (int j = 0; j < i; ++j) {
            // Potential energy.
            energy -=
                GravitationalConstant * bodies[i].mass() * bodies[j].mass() /
                    (system_state.positions[i].value -
                     system_state.positions[j].value).Norm();
          }
        }
        e_error = std::max(e_error, AbsoluteError(initial_energy, energy));
      }
      // We plot the maximum error, i.e., the L∞ norm of the error.
      // Blanes and Moan (2002), or Blanes, Casas and Ros (2001) tend to use
      // the average error (the normalized L¹ norm) instead.
      q_errors.emplace_back(q_error);
      v_errors.emplace_back(v_error);
      e_errors.emplace_back(e_error);
      evaluations.emplace_back(number_of_evaluations);
    }
    q_error_data.emplace_back(PlottableDataset(evaluations, q_errors));
    v_error_data.emplace_back(PlottableDataset(evaluations, v_errors));
    e_error_data.emplace_back(PlottableDataset(evaluations, e_errors));
    names.emplace_back(Escape(method.name));
  }
  std::ofstream file;
  file.open("planets_graphs.generated.wl");
  file << Assign("qErrorData", q_error_data);
  file << Assign("vErrorData", v_error_data);
  file << Assign("eErrorData", e_error_data);
  file << Assign("names", names);
  file.close();
}

}  // namespace mathematica
}  // namespace principia
