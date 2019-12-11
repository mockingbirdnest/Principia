
#include "mathematica/integrator_plots.hpp"

#include <algorithm>
#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <vector>

#include "base/bundle.hpp"
#include "base/file.hpp"
#include "base/status.hpp"
#include "glog/logging.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "mathematica/mathematica.hpp"
#include "serialization/integrators.pb.h"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/numerics.hpp"

#define SLMS_INTEGRATOR(name)                         \
  {                                                   \
    (integrators::SymmetricLinearMultistepIntegrator< \
        integrators::methods::name,                   \
        Length>()),                                   \
        u8###name, 1                                  \
  }
#define SRKN_INTEGRATOR(name)                                 \
  {                                                           \
    (integrators::SymplecticRungeKuttaNyströmIntegrator<      \
        integrators::methods::name,                           \
        Length>()),                                           \
        u8###name, (integrators::methods::name::evaluations)  \
  }
#define SPRK_INTEGRATOR(name, composition)                   \
  {                                                          \
    (integrators::SymplecticRungeKuttaNyströmIntegrator<     \
        integrators::methods::name,                          \
        serialization::FixedStepSizeIntegrator::composition, \
        Length>()),                                          \
        u8###name " " u8###composition,                      \
        (integrators::methods::name::evaluations)            \
  }

namespace principia {

using base::Bundle;
using base::not_null;
using base::OFStream;
using base::Status;
using geometry::BarycentreCalculator;
using geometry::Displacement;
using geometry::Inertial;
using geometry::InnerProduct;
using geometry::Instant;
using geometry::Velocity;
using integrators::FixedStepSizeIntegrator;
using integrators::IntegrationProblem;
using integrators::SpecialSecondOrderDifferentialEquation;
using numerics::DoublePrecision;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Energy;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Pow;
using quantities::Sin;
using quantities::SIUnit;
using quantities::Speed;
using quantities::SpecificEnergy;
using quantities::Stiffness;
using quantities::Time;
using quantities::si::Joule;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using testing_utilities::AbsoluteError;
using testing_utilities::ComputeHarmonicOscillatorAcceleration1D;
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
          SPRK_INTEGRATOR(NewtonDelambreStørmerVerletLeapfrog, BAB),
          SPRK_INTEGRATOR(McLachlanAtela1992Order2Optimal, BA),
          SPRK_INTEGRATOR(McLachlan1995S2, BAB),
          // Order 3
          SPRK_INTEGRATOR(Ruth1983, BA),
          SPRK_INTEGRATOR(McLachlanAtela1992Order3Optimal, BA),
          // Order 4
          SPRK_INTEGRATOR(CandyRozmus1991ForestRuth1990, BAB),
          SPRK_INTEGRATOR(鈴木1990, BAB),
          SPRK_INTEGRATOR(McLachlan1995SS5, BAB),
          SPRK_INTEGRATOR(McLachlan1995S4, BAB),
          SPRK_INTEGRATOR(McLachlan1995S5, BAB),
          SPRK_INTEGRATOR(BlanesMoan2002S6, BAB),
          SRKN_INTEGRATOR(McLachlanAtela1992Order4Optimal),
          SRKN_INTEGRATOR(McLachlan1995SB3A4),
          SRKN_INTEGRATOR(McLachlan1995SB3A5),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN6B),
          // Order 5
          SRKN_INTEGRATOR(McLachlanAtela1992Order5Optimal),
          // Order 6
          SPRK_INTEGRATOR(吉田1990Order6A, BAB),
          SPRK_INTEGRATOR(吉田1990Order6B, BAB),
          SPRK_INTEGRATOR(吉田1990Order6C, BAB),
          SPRK_INTEGRATOR(McLachlan1995SS9, BAB),
          SPRK_INTEGRATOR(BlanesMoan2002S10, BAB),
          SRKN_INTEGRATOR(OkunborSkeel1994Order6Method13),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN11B),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN14A),
          // Order 8
          SPRK_INTEGRATOR(吉田1990Order8A, BAB),
          SPRK_INTEGRATOR(吉田1990Order8B, BAB),
          SPRK_INTEGRATOR(吉田1990Order8C, BAB),
          SPRK_INTEGRATOR(吉田1990Order8D, BAB),
          SPRK_INTEGRATOR(吉田1990Order8E, BAB),
          SPRK_INTEGRATOR(McLachlan1995SS15, BAB),
          SPRK_INTEGRATOR(McLachlan1995SS17, BAB),
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
          SPRK_INTEGRATOR(McLachlan1995SS9, BAB),
          SPRK_INTEGRATOR(BlanesMoan2002S10, BAB),
          SRKN_INTEGRATOR(OkunborSkeel1994Order6Method13),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN11B),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN14A),
          // Order 8
          SPRK_INTEGRATOR(McLachlan1995SS15, BAB),
          SPRK_INTEGRATOR(McLachlan1995SS17, BAB)};
}

}  // namespace

// Templatized to allow for problems where specific energy is more convenient
// than energy.
template<typename Energy>
class WorkErrorGraphGenerator {
 public:
  struct Errors {
    Length q_error;
    Speed v_error;
    Energy e_error;
  };

  WorkErrorGraphGenerator(
      std::function<Status(Instant const& t,
                           std::vector<Length> const& q,
                           std::vector<Acceleration>& result,
                           int* evaluations)> compute_accelerations,
      ODE::SystemState initial_state,
      std::function<Errors(ODE::SystemState const&)> compute_errors,
      Instant const& tmax,
      std::string const& problem_name)
      : methods_(Methods()),
        compute_accelerations_(compute_accelerations),
        initial_state_(initial_state),
        compute_errors_(compute_errors),
        tmax_(tmax),
        problem_name_(problem_name) {
    q_errors_.resize(methods_.size());
    v_errors_.resize(methods_.size());
    e_errors_.resize(methods_.size());
    evaluations_.resize(methods_.size());
    for (int i = 0; i < methods_.size(); ++i) {
      q_errors_[i].resize(integrations_per_integrator_);
      v_errors_[i].resize(integrations_per_integrator_);
      e_errors_[i].resize(integrations_per_integrator_);
      evaluations_[i].resize(integrations_per_integrator_);
    }
  }

  std::string GetMathematicaData() {
    LOG(INFO) << "Using " << std::thread::hardware_concurrency()
              << " worker threads";
    Bundle bundle;
    for (int method_index = 0; method_index < methods_.size(); ++method_index) {
      for (int time_step_index = 0;
           time_step_index < integrations_per_integrator_;
           ++time_step_index) {
        bundle.Add(std::bind(&WorkErrorGraphGenerator::Integrate,
                             this,
                             method_index,
                             time_step_index));
      }
    }
    bundle.Join();

    std::vector<std::string> q_error_data;
    std::vector<std::string> v_error_data;
    std::vector<std::string> e_error_data;
    std::vector<std::string> names;
    for (int i = 0; i < methods_.size(); ++i) {
      q_error_data.emplace_back(
          PlottableDataset(evaluations_[i], q_errors_[i]));
      v_error_data.emplace_back(
          PlottableDataset(evaluations_[i], v_errors_[i]));
      e_error_data.emplace_back(
          PlottableDataset(evaluations_[i], e_errors_[i]));
      names.emplace_back(Escape(methods_[i].name));
    }
    std::string result;
    result += Assign("qErrorData", q_error_data);
    result += Assign("vErrorData", v_error_data);
    result += Assign("eErrorData", e_error_data);
    result += Assign("names", names);
    return result;
  }

 private:
  Status Integrate(int const method_index, int const time_step_index) {
    auto const& method = methods_[method_index];
    Problem problem;
    int number_of_evaluations = 0;
    problem.equation.compute_acceleration = std::bind(
        compute_accelerations_, _1, _2, _3, &number_of_evaluations);
    problem.initial_state = initial_state_;
    auto const t0 = problem.initial_state.time.value;
    Length max_q_error;
    Speed max_v_error;
    Energy max_e_error;
    auto append_state = [this, &max_q_error, &max_v_error, &max_e_error](
      ODE::SystemState const& state) {
      auto const errors = compute_errors_(state);
      max_q_error = std::max(max_q_error, errors.q_error);
      max_v_error = std::max(max_v_error, errors.v_error);
      max_e_error = std::max(max_e_error, errors.e_error);
    };
    Time const Δt = method.evaluations * starting_step_size_per_evaluation_ /
                    std::pow(step_reduction_, time_step_index);
    auto const instance =
        method.integrator.NewInstance(problem, append_state, Δt);

    instance->Solve(tmax_);
    // Log both the actual number of evaluations and a theoretical number
    // that ignores any startup costs; that theoretical number is the one
    // used for plotting.
    int const amortized_evaluations =
        method.evaluations * static_cast<int>(std::floor((tmax_ - t0) / Δt));
    LOG_EVERY_N(INFO, 50) << "[" << method_index << "," << time_step_index
                          << "] " << problem_name_ << ": "
                          << number_of_evaluations << " actual evaluations ("
                          << amortized_evaluations << " amortized) with "
                          << method.name;
    // We plot the maximum error, i.e., the L∞ norm of the error.
    // Blanes and Moan (2002), or Blanes, Casas and Ros (2001) tend to use
    // the average error (the normalized L¹ norm) instead.
    q_errors_[method_index][time_step_index] = max_q_error;
    v_errors_[method_index][time_step_index] = max_v_error;
    e_errors_[method_index][time_step_index] = max_e_error;
    evaluations_[method_index][time_step_index] = amortized_evaluations;

    return Status::OK;
  }

  std::vector<SimpleHarmonicMotionPlottedIntegrator> const methods_;
  std::function<Status(Instant const& t,
                       std::vector<Length> const& q,
                       std::vector<Acceleration>& result,
                       int* evaluations)>
      compute_accelerations_;
  ODE::SystemState initial_state_;
  std::function<Errors(ODE::SystemState const&)> compute_errors_;
  std::vector<std::vector<Length>> q_errors_;
  std::vector<std::vector<Speed>> v_errors_;
  std::vector<std::vector<Energy>> e_errors_;
  std::vector<std::vector<double>> evaluations_;
  Instant const tmax_;
  std::string const problem_name_;
  double const step_reduction_ = 1.015;
  Time const starting_step_size_per_evaluation_ = 1 * Second;
  int const integrations_per_integrator_ = 500;
};


void GenerateSimpleHarmonicMotionWorkErrorGraphs() {
  ODE::SystemState initial_state;
  Instant const t0;
  Length const q_amplitude = 1 * Metre;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Stiffness const k = SIUnit<Stiffness>();
  Mass const m = 1 * Kilogram;

  initial_state.positions.emplace_back(q_amplitude);
  initial_state.velocities.emplace_back(0 * Metre / Second);
  initial_state.time = DoublePrecision<Instant>(t0);

  Instant const tmax = t0 + 50 * Second;
  auto const compute_error = [q_amplitude, v_amplitude, ω, m, k, t0](
      ODE::SystemState const& state) {
    return WorkErrorGraphGenerator<Energy>::Errors{
        AbsoluteError(q_amplitude * Cos(ω * (state.time.value - t0)),
                      state.positions[0].value),
        AbsoluteError(-v_amplitude * Sin(ω * (state.time.value - t0)),
                      state.velocities[0].value),
        AbsoluteError(0.5 * Joule,
                      (m * Pow<2>(state.velocities[0].value) +
                       k * Pow<2>(state.positions[0].value)) / 2)};
  };
  WorkErrorGraphGenerator<Energy> generator(
      ComputeHarmonicOscillatorAcceleration1D,
      initial_state,
      compute_error,
      tmax,
      "Harmonic oscillator");

  OFStream file(TEMP_DIR / "simple_harmonic_motion_graphs.generated.wl");
  file << generator.GetMathematicaData();
}

void GenerateKeplerProblemWorkErrorGraphs(double const eccentricity) {
  ODE::SystemState initial_state;
  Instant const t0;
  GravitationalParameter const μ = SIUnit<GravitationalParameter>();
  MassiveBody b1(μ);
  MasslessBody b2;

  using World = geometry::Frame<serialization::Frame::TestTag,
                                serialization::Frame::TEST,
                                Inertial>;
  KeplerianElements<World> elements;
  elements.semimajor_axis = 1 * Metre;
  elements.eccentricity = eccentricity;
  KeplerOrbit<World> const orbit(b1, b2, elements, t0);

  auto const initial_dof = orbit.StateVectors(t0);
  CHECK_EQ(initial_dof.displacement().coordinates().z, 0 * Metre);
  CHECK_EQ(initial_dof.velocity().coordinates().z, 0 * Metre / Second);

  initial_state.positions.emplace_back(
      initial_dof.displacement().coordinates().x);
  initial_state.positions.emplace_back(
      initial_dof.displacement().coordinates().y);
  initial_state.velocities.emplace_back(initial_dof.velocity().coordinates().x);
  initial_state.velocities.emplace_back(initial_dof.velocity().coordinates().y);
  initial_state.time = DoublePrecision<Instant>(t0);

  // Integrate over 8 orbits.
  Instant const tmax =
      t0 + 8 * (2 * π * Radian) / *orbit.elements_at_epoch().mean_motion;

  SpecificEnergy const initial_specific_energy =
      initial_dof.velocity().Norm²() / 2 -
      μ / initial_dof.displacement().Norm();

  auto const compute_error = [&orbit, μ, initial_specific_energy](
      ODE::SystemState const& state) {
    Displacement<World> q(
        {state.positions[0].value, state.positions[1].value, 0 * Metre});
    Velocity<World> v({state.velocities[0].value,
                       state.velocities[1].value,
                       0 * Metre / Second});
    auto const expected_dof = orbit.StateVectors(state.time.value);
    return WorkErrorGraphGenerator<SpecificEnergy>::Errors{
        AbsoluteError(expected_dof.displacement(), q),
        AbsoluteError(expected_dof.velocity(), v),
        AbsoluteError(initial_specific_energy, v.Norm²() / 2 - μ / q.Norm())};
  };

  WorkErrorGraphGenerator<SpecificEnergy> generator(
      ComputeKeplerAcceleration,
      initial_state,
      compute_error,
      tmax,
      " Kepler problem with e = " + std::to_string(eccentricity));

  OFStream file(TEMP_DIR / ("kepler_problem_graphs_" +
                            std::to_string(eccentricity) + ".generated.wl"));
  file << generator.GetMathematicaData();
}

}  // namespace mathematica
}  // namespace principia
