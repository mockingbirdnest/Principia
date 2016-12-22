
#include "mathematica/integrator_plots.hpp"

#include <algorithm>
#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <vector>

#include "glog/logging.h"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "mathematica/mathematica.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/numerics.hpp"

#define SLMS_INTEGRATOR(name) \
  { (integrators::name<Length>()), #name, 1 }
#define SRKN_INTEGRATOR(name)                     \
  {                                               \
    (integrators::name<Length>()), #name,         \
        (integrators::name<Length>()).evaluations \
  }
#define SPRK_INTEGRATOR(name, composition)                 \
  {                                                        \
    (integrators::name<Length, Speed>()                    \
         .AsRungeKuttaNyströmIntegrator<(composition)>()), \
        #name " " #composition,                            \
        (integrators::name<Length, Speed>()).evaluations   \
  }

namespace principia {

using base::not_null;
using geometry::InnerProduct;
using geometry::BarycentreCalculator;
using geometry::Velocity;
using integrators::ABA;
using integrators::BA;
using integrators::FixedStepSizeIntegrator;
using integrators::IntegrationProblem;
using integrators::SpecialSecondOrderDifferentialEquation;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Energy;
using quantities::si::Joule;
using quantities::si::Kilogram;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Minute;
using quantities::Pow;
using quantities::si::Radian;
using quantities::si::Second;
using quantities::Sin;
using quantities::Speed;
using quantities::SpecificEnergy;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using testing_utilities::AbsoluteError;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
using testing_utilities::ComputeKeplerAcceleration;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace mathematica {

// TODO(egg): it would probably be saner to use Position<Whatever> and make the
// simple harmonic oscillator work in 3d.
using ODE = SpecialSecondOrderDifferentialEquation<Length>;
using Problem = IntegrationProblem<ODE>;

namespace {

struct SimpleHarmonicMotionPlottedIntegrator final {
  FixedStepSizeIntegrator<ODE> const& integrator;
  std::string name;
  int evaluations;
};

// This list should be sorted by:
// 1. increasing order;
// 2. SPRKs before SRKNs before symmetric linear multistep integrators.
// 3. increasing number of evaluations;
// 4. date;
// 5. author names;
// 6. method name.
std::vector<SimpleHarmonicMotionPlottedIntegrator> Methods() {
  return {// Order 2
          SPRK_INTEGRATOR(NewtonDelambreStørmerVerletLeapfrog, ABA),
          SPRK_INTEGRATOR(McLachlanAtela1992Order2Optimal, BA),
          SPRK_INTEGRATOR(McLachlan1995S2, ABA),
          // Order 3
          SPRK_INTEGRATOR(Ruth1983, BA),
          SPRK_INTEGRATOR(McLachlanAtela1992Order3Optimal, BA),
          // Order 4
          SPRK_INTEGRATOR(CandyRozmus1991ForestRuth1990, ABA),
          SPRK_INTEGRATOR(Suzuki1990, ABA),
          SPRK_INTEGRATOR(McLachlan1995SS5, ABA),
          SPRK_INTEGRATOR(McLachlan1995S4, ABA),
          SPRK_INTEGRATOR(McLachlan1995S5, ABA),
          SPRK_INTEGRATOR(BlanesMoan2002S6, ABA),
          SRKN_INTEGRATOR(McLachlanAtela1992Order4Optimal),
          SRKN_INTEGRATOR(McLachlan1995SB3A4),
          SRKN_INTEGRATOR(McLachlan1995SB3A5),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN6B),
          // Order 5
          SRKN_INTEGRATOR(McLachlanAtela1992Order5Optimal),
          // Order 6
          SPRK_INTEGRATOR(Yoshida1990Order6A, ABA),
          SPRK_INTEGRATOR(Yoshida1990Order6B, ABA),
          SPRK_INTEGRATOR(Yoshida1990Order6C, ABA),
          SPRK_INTEGRATOR(McLachlan1995SS9, ABA),
          SPRK_INTEGRATOR(BlanesMoan2002S10, ABA),
          SRKN_INTEGRATOR(OkunborSkeel1994Order6Method13),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN11B),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN14A),
          // Order 8
          SPRK_INTEGRATOR(Yoshida1990Order8A, ABA),
          SPRK_INTEGRATOR(Yoshida1990Order8B, ABA),
          SPRK_INTEGRATOR(Yoshida1990Order8C, ABA),
          SPRK_INTEGRATOR(Yoshida1990Order8D, ABA),
          SPRK_INTEGRATOR(Yoshida1990Order8E, ABA),
          SPRK_INTEGRATOR(McLachlan1995SS15, ABA),
          SPRK_INTEGRATOR(McLachlan1995SS17, ABA),
          SLMS_INTEGRATOR(QuinlanTremaine1990Order8),
          SLMS_INTEGRATOR(Quinlan1999Order8A),
          SLMS_INTEGRATOR(Quinlan1999Order8B),
          // Order 10
          SLMS_INTEGRATOR(QuinlanTremaine1990Order10),
          // Order 12
          SLMS_INTEGRATOR(QuinlanTremaine1990Order12),
          // Order 14
          SLMS_INTEGRATOR(QuinlanTremaine1990Order14)};
}

// Those methods which have converged to the limits of double-precision floating
// point error on the circular Kepler problem tested by
// |GenerateKeplerProblemWorkErrorGraphs| with less than 8e4 evaluations.
std::vector<SimpleHarmonicMotionPlottedIntegrator> ReferenceMethods() {
  return {// Order 5
          SRKN_INTEGRATOR(McLachlanAtela1992Order5Optimal),
          // Order 6
          SPRK_INTEGRATOR(McLachlan1995SS9, ABA),
          SPRK_INTEGRATOR(BlanesMoan2002S10, ABA),
          SRKN_INTEGRATOR(OkunborSkeel1994Order6Method13),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN11B),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN14A),
          // Order 8
          SPRK_INTEGRATOR(McLachlan1995SS15, ABA),
          SPRK_INTEGRATOR(McLachlan1995SS17, ABA)};
}

}  // namespace

void GenerateSimpleHarmonicMotionWorkErrorGraphs() {
  ODE::SystemState initial_state;
  Problem problem;
  int number_of_evaluations;
  problem.equation.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration,
                _1, _2, _3, &number_of_evaluations);
  problem.initial_state = &initial_state;

  Instant const t0;
  Length const q_amplitude = 1 * Metre;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Stiffness const k = SIUnit<Stiffness>();
  Mass const m = 1 * Kilogram;

  initial_state.positions.emplace_back(q_amplitude);
  initial_state.velocities.emplace_back(0 * Metre / Second);
  initial_state.time = t0;

  Instant const tmax = t0 + 50 * Second;

  double const step_reduction = 1.015;


  std::vector<std::string> q_error_data;
  std::vector<std::string> v_error_data;
  std::vector<std::string> e_error_data;

  Length max_q_error;
  Speed max_v_error;
  Energy max_e_error;
  auto append_state = [&max_q_error, &max_v_error, &max_e_error,
                       q_amplitude, v_amplitude, ω, m, k, t0](
      ODE::SystemState const& state) {
    max_q_error = std::max(max_q_error,
        AbsoluteError(q_amplitude * Cos(ω * (state.time.value - t0)),
                      state.positions[0].value));
    max_v_error = std::max(max_v_error,
        AbsoluteError(-v_amplitude * Sin(ω * (state.time.value - t0)),
                      state.velocities[0].value));
    max_e_error = std::max(max_e_error,
        AbsoluteError(0.5 * Joule,
                      (m * Pow<2>(state.velocities[0].value) +
                       k * Pow<2>(state.positions[0].value)) / 2));
  };

  std::vector<std::string> names;
  for (auto const& method : Methods()) {
    LOG(INFO) << "Harmonic oscillator: " << method.name;
    Time Δt = method.evaluations * 1 * Second;
    std::vector<Length> q_errors;
    std::vector<Speed> v_errors;
    std::vector<Energy> e_errors;
    std::vector<double> evaluations;
    for (int i = 0; i < 500; ++i, Δt /= step_reduction) {
      max_q_error = Length{};
      max_v_error = Speed{};
      max_e_error = Energy{};
      number_of_evaluations = 0;
      auto instance = method.integrator.NewInstance(problem, append_state, Δt);
      method.integrator.Solve(tmax, *instance);
      // Log both the actual number of evaluations and a theoretical number that
      // ignores any startup costs; that theoretical number is the one used for
      // plotting.
      int const amortized_evaluations =
          method.evaluations * static_cast<int>(std::floor((tmax - t0) / Δt));
      LOG_IF(INFO, (i + 1) % 50 == 0) << number_of_evaluations << "("
                                      << amortized_evaluations << ")";
      // We plot the maximum error, i.e., the L∞ norm of the error.
      // Blanes and Moan (2002), or Blanes, Casas and Ros (2001) tend to use
      // the average error (the normalized L¹ norm) instead.
      q_errors.emplace_back(max_q_error);
      v_errors.emplace_back(max_v_error);
      e_errors.emplace_back(max_e_error);
      evaluations.emplace_back(amortized_evaluations);
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

void GenerateKeplerProblemWorkErrorGraphs(double const eccentricity) {
  ODE::SystemState initial_state;
  Problem problem;
  int number_of_evaluations;
  problem.equation.compute_acceleration =
      std::bind(ComputeKeplerAcceleration, _1, _2, _3, &number_of_evaluations);
  problem.initial_state = &initial_state;

  Instant const t0;

  GravitationalParameter const μ = SIUnit<GravitationalParameter>();
  MassiveBody b1(μ);
  MasslessBody b2;

  using World = geometry::Frame<serialization::Frame::TestTag,
                                serialization::Frame::TEST,
                                /*frame_is_inertial=*/true>;
  KeplerianElements<World> elements;
  elements.semimajor_axis = 1 * Metre;
  elements.eccentricity = eccentricity;
  KeplerOrbit<World> orbit(b1, b2, elements, t0);

  auto const initial_dof = orbit.StateVectors(t0);
  CHECK_EQ(initial_dof.displacement().coordinates().z, 0 * Metre);
  CHECK_EQ(initial_dof.velocity().coordinates().z, 0 * Metre / Second);

  initial_state.positions.emplace_back(
      initial_dof.displacement().coordinates().x);
  initial_state.positions.emplace_back(
      initial_dof.displacement().coordinates().y);
  initial_state.velocities.emplace_back(initial_dof.velocity().coordinates().x);
  initial_state.velocities.emplace_back(initial_dof.velocity().coordinates().y);
  initial_state.time = t0;

  // Integrate over 8 orbits.
  Instant const tmax =
      t0 + 8 * (2 * π * Radian) / *orbit.elements_at_epoch().mean_motion;

  SpecificEnergy const initial_specific_energy =
      InnerProduct(initial_dof.velocity(), initial_dof.velocity()) / 2 -
      μ / initial_dof.displacement().Norm();

  double const step_reduction = 1.015;

  std::vector<std::string> q_error_data;
  std::vector<std::string> v_error_data;
  std::vector<std::string> e_error_data;

  Length max_q_error;
  Speed max_v_error;
  SpecificEnergy max_e_error;
  auto append_state = [&max_q_error, &max_v_error, &max_e_error, &orbit,
                       μ, initial_specific_energy](
      ODE::SystemState const& state) {
    Displacement<World> q(
        {state.positions[0].value, state.positions[1].value, 0 * Metre});
    Velocity<World> v({state.velocities[0].value,
                       state.velocities[1].value,
                       0 * Metre / Second});
    auto const expected_dof = orbit.StateVectors(state.time.value);
    max_q_error = std::max(max_q_error,
        AbsoluteError(expected_dof.displacement(), q));
    max_v_error = std::max(max_v_error,
        AbsoluteError(expected_dof.velocity(), v));
    max_e_error = std::max(max_e_error,
        AbsoluteError(initial_specific_energy,
                      InnerProduct(v, v) / 2 - μ / q.Norm()));
  };

  std::vector<std::string> names;
  for (auto const& method : Methods()) {
    LOG(INFO) << " Kepler problem with e = " << eccentricity << ": "
              << method.name;
    Time Δt = method.evaluations * 1 * Second;
    std::vector<Length> q_errors;
    std::vector<Speed> v_errors;
    std::vector<SpecificEnergy> e_errors;
    std::vector<double> evaluations;
    for (int i = 0; i < 500; ++i, Δt /= step_reduction) {
      max_q_error = Length{};
      max_v_error = Speed{};
      max_e_error = SpecificEnergy{};
      number_of_evaluations = 0;
      auto instance = method.integrator.NewInstance(problem, append_state, Δt);
      method.integrator.Solve(tmax, *instance);
      // Log both the actual number of evaluations and a theoretical number that
      // ignores any startup costs; that theoretical number is the one used for
      // plotting.
      int const amortized_evaluations =
          method.evaluations * static_cast<int>(std::floor((tmax - t0) / Δt));
      LOG_IF(INFO, (i + 1) % 50 == 0) << number_of_evaluations << "("
                                      << amortized_evaluations << ")";
      // We plot the maximum error, i.e., the L∞ norm of the error.
      // Blanes and Moan (2002), or Blanes, Casas and Ros (2001) tend to use
      // the average error (the normalized L¹ norm) instead.
      q_errors.emplace_back(max_q_error);
      v_errors.emplace_back(max_v_error);
      e_errors.emplace_back(max_e_error);
      evaluations.emplace_back(amortized_evaluations);
    }
    q_error_data.emplace_back(PlottableDataset(evaluations, q_errors));
    v_error_data.emplace_back(PlottableDataset(evaluations, v_errors));
    e_error_data.emplace_back(PlottableDataset(evaluations, e_errors));
    names.emplace_back(Escape(method.name));
  }
  std::ofstream file;
  file.open("kepler_problem_graphs_" + std::to_string(eccentricity) +
            ".generated.wl");
  file << Assign("qErrorData", q_error_data);
  file << Assign("vErrorData", v_error_data);
  file << Assign("eErrorData", e_error_data);
  file << Assign("names", names);
  file.close();
}

}  // namespace mathematica
}  // namespace principia
