#include "mathematica/integrator_plots.hpp"

#include <algorithm>
#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "base/bundle.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "glog/logging.h"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/explicit_runge_kutta_integrator.hpp"
#include "integrators/integrators.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/double_precision.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/integrators.pb.h"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/numerics.hpp"

#define SLMS_INTEGRATOR(name)                                                 \
  {static_cast<FixedStepSizeIntegrator<SecondOrderODE> const*>(               \
       &SymmetricLinearMultistepIntegrator<methods::name, SecondOrderODE>()), \
   u8## #name,                                                                \
   {methods::name::order},                                                    \
   "SLMS",                                                                    \
   1}
#define SRKN_INTEGRATOR(name)                                     \
  {static_cast<FixedStepSizeIntegrator<SecondOrderODE> const*>(   \
       &SymplecticRungeKuttaNyströmIntegrator<methods::name,      \
                                              SecondOrderODE>()), \
   u8## #name,                                                    \
   {methods::name::order},                                        \
   "SRKN",                                                        \
   (methods::name::evaluations)}
#define ERK_INTEGRATOR(name)                                           \
  {static_cast<FixedStepSizeIntegrator<FirstOrderODE> const*>(         \
       &ExplicitRungeKuttaIntegrator<methods::name, FirstOrderODE>()), \
   u8## #name,                                                         \
   {methods::name::order},                                             \
   "ERK",                                                              \
   methods::name::first_same_as_last ? methods::name::stages - 1       \
                                     : methods::name::stages}
#define EERK_INTEGRATOR(name)                                                  \
  {static_cast<AdaptiveStepSizeIntegrator<FirstOrderODE> const*>(              \
       &EmbeddedExplicitRungeKuttaIntegrator<methods::name, FirstOrderODE>()), \
   u8## #name,                                                                 \
   {methods::name::lower_order, methods::name::higher_order},                  \
   "EERK",                                                                     \
   0}
#define EERKN_INTEGRATOR(name)                                          \
  {static_cast<AdaptiveStepSizeIntegrator<SecondOrderODE> const*>(      \
       &EmbeddedExplicitRungeKuttaNyströmIntegrator<methods::name,      \
                                                    SecondOrderODE>()), \
   u8## #name,                                                          \
   {methods::name::lower_order, methods::name::higher_order},           \
   std::is_base_of_v<EmbeddedExplicitGeneralizedRungeKuttaNyström,      \
                     methods::name>                                     \
       ? "EEGRKN"                                                       \
       : "EERKN",                                                       \
   0}
#define SPRK_INTEGRATOR(name)                                   \
  {static_cast<FixedStepSizeIntegrator<SecondOrderODE> const*>( \
       &SymplecticRungeKuttaNyströmIntegrator<                  \
           methods::name,                                       \
           serialization::FixedStepSizeIntegrator::BA,          \
           SecondOrderODE>()),                                  \
   u8## #name " BA",                                            \
   {methods::name::order},                                      \
   "SPRK",                                                      \
   (methods::name::evaluations)}
#define SPRK_INTEGRATOR_FSAL(name)                               \
  {static_cast<FixedStepSizeIntegrator<SecondOrderODE> const*>(  \
       &SymplecticRungeKuttaNyströmIntegrator<                   \
           methods::name,                                        \
           serialization::FixedStepSizeIntegrator::ABA,          \
           SecondOrderODE>()),                                   \
   u8## #name " ABA",                                            \
   {methods::name::order},                                       \
   "SPRK",                                                       \
   (methods::name::evaluations)},                                \
  {                                                              \
    static_cast<FixedStepSizeIntegrator<SecondOrderODE> const*>( \
        &SymplecticRungeKuttaNyströmIntegrator<                  \
            methods::name,                                       \
            serialization::FixedStepSizeIntegrator::BAB,         \
            SecondOrderODE>()),                                  \
        u8## #name " BAB", {methods::name::order}, "SPRK",       \
        (methods::name::evaluations)                             \
  }

namespace principia {
namespace mathematica {
namespace _integrator_plots {
namespace internal {

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using namespace principia::base::_bundle;
using namespace principia::base::_file;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_explicit_runge_kutta_integrator;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::integrators::_symplectic_partitioned_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::mathematica::_mathematica;
using namespace principia::numerics::_double_precision;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_integration;
using namespace principia::testing_utilities::_numerics;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;
using SecondOrderODE = SpecialSecondOrderDifferentialEquation<Position<World>>;
using FirstOrderODE =
    ExplicitFirstOrderOrdinaryDifferentialEquation<Instant,
                                                   Position<World>,
                                                   Velocity<World>>;
using Problem = InitialValueProblem<SecondOrderODE>;

namespace {

struct PlottedIntegrator final {
  std::variant<FixedStepSizeIntegrator<SecondOrderODE> const*,
               AdaptiveStepSizeIntegrator<SecondOrderODE> const*,
               FixedStepSizeIntegrator<FirstOrderODE> const*,
               AdaptiveStepSizeIntegrator<FirstOrderODE> const*> const integrator;
  std::string name;
  std::vector<int> orders;
  std::string method_type;
  int evaluations;
};

// This list should be sorted by:
// 1. increasing convergence order;
// 2. method types in the following order:
//    a. SPRKs,
//    b. SRKNs,
//    c. symmetric linear multistep integrators,
//    d. Explicit Runge-Kutta methods,
//    e. EERKs,
//    f. EEGRKNs,
//    g. EERKNs;
// 3. increasing number of evaluations;
// 4. date;
// 5. author names;
// 6. method name.
std::vector<PlottedIntegrator> Methods() {
  PlottedIntegrator meow =
      SRKN_INTEGRATOR(BlanesMoan2002SRKN6B);
  return {// Order 2
          SPRK_INTEGRATOR_FSAL(NewtonDelambreStørmerVerletLeapfrog),
          SPRK_INTEGRATOR(McLachlanAtela1992Order2Optimal),
          SPRK_INTEGRATOR_FSAL(McLachlan1995S2),
          // Order 3
          SPRK_INTEGRATOR(Ruth1983),
          SPRK_INTEGRATOR(McLachlanAtela1992Order3Optimal),
          EERKN_INTEGRATOR(Fine1987RKNG34),
          EERKN_INTEGRATOR(DormandالمكاوىPrince1986RKN434FM),
          // Order 4
          SPRK_INTEGRATOR_FSAL(CandyRozmus1991ForestRuth1990),
          SPRK_INTEGRATOR_FSAL(鈴木1990),
          SPRK_INTEGRATOR_FSAL(McLachlan1995SS5),
          SPRK_INTEGRATOR_FSAL(McLachlan1995S4),
          SPRK_INTEGRATOR_FSAL(McLachlan1995S5),
          SPRK_INTEGRATOR_FSAL(BlanesMoan2002S6),
          SRKN_INTEGRATOR(McLachlanAtela1992Order4Optimal),
          SRKN_INTEGRATOR(McLachlan1995SB3A4),
          SRKN_INTEGRATOR(McLachlan1995SB3A5),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN6B),
          ERK_INTEGRATOR(Kutta1901Vσ1),
          EERK_INTEGRATOR(DormandPrince1986RK547FC),
          EERKN_INTEGRATOR(Fine1987RKNG45),
          // Order 5
          SRKN_INTEGRATOR(McLachlanAtela1992Order5Optimal),
          // Order 6
          SPRK_INTEGRATOR_FSAL(吉田1990Order6A),
          SPRK_INTEGRATOR_FSAL(吉田1990Order6B),
          SPRK_INTEGRATOR_FSAL(吉田1990Order6C),
          SPRK_INTEGRATOR_FSAL(McLachlan1995SS9),
          SPRK_INTEGRATOR_FSAL(BlanesMoan2002S10),
          SRKN_INTEGRATOR(OkunborSkeel1994Order6Method13),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN11B),
          SRKN_INTEGRATOR(BlanesMoan2002SRKN14A),
          // Order 8
          SPRK_INTEGRATOR_FSAL(吉田1990Order8A),
          SPRK_INTEGRATOR_FSAL(吉田1990Order8B),
          SPRK_INTEGRATOR_FSAL(吉田1990Order8C),
          SPRK_INTEGRATOR_FSAL(吉田1990Order8D),
          SPRK_INTEGRATOR_FSAL(吉田1990Order8E),
          SPRK_INTEGRATOR_FSAL(McLachlan1995SS15),
          SPRK_INTEGRATOR_FSAL(McLachlan1995SS17),
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
      std::function<
          absl::Status(Instant const& t,
                       std::vector<Position<World>> const& q,
                       std::vector<Vector<Acceleration, World>>& result,
                       int* evaluations)> compute_accelerations,
      SecondOrderODE::State initial_state,
      std::function<Errors(SecondOrderODE::State const&)> compute_errors,
      std::vector<Instant> const& tmax,
      Length const& first_tolerance,
      Time const& min_step_per_evaluation,
      std::string problem_name)
      : methods_(Methods()),
        compute_accelerations_(std::move(compute_accelerations)),
        initial_state_(std::move(initial_state)),
        compute_errors_(std::move(compute_errors)),
        tmax_(tmax),
        first_tolerance_(first_tolerance),
        problem_name_(std::move(problem_name)),
        step_reduction_(std::exp(std::log(starting_step_size_per_evaluation_ /
                                          min_step_per_evaluation) /
                                 integrations_per_integrator_)) {
    q_errors_.resize(methods_.size());
    v_errors_.resize(methods_.size());
    e_errors_.resize(methods_.size());
    evaluations_.resize(methods_.size());
    for (int i = 0; i < methods_.size(); ++i) {
      q_errors_[i].resize(integrations_per_integrator_);
      v_errors_[i].resize(integrations_per_integrator_);
      e_errors_[i].resize(integrations_per_integrator_);
      evaluations_[i].resize(integrations_per_integrator_);
      for (int j = 0; j < integrations_per_integrator_; ++j) {
        q_errors_[i][j].resize(tmax_.size());
        v_errors_[i][j].resize(tmax_.size());
        e_errors_[i][j].resize(tmax_.size());
        evaluations_[i][j].resize(tmax_.size());
      }
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
    CHECK_OK(bundle.Join());

    std::vector<std::vector<std::vector<std::tuple<double, Length>>>>
        q_error_data;
    std::vector<std::vector<std::vector<std::tuple<double, Speed>>>>
        v_error_data;
    std::vector<std::vector<std::vector<std::tuple<double, Energy>>>>
        e_error_data;
    std::vector<std::string> names;
    std::vector<std::vector<int>> orders;
    std::vector<std::string> types;
    for (int i = 0; i < methods_.size(); ++i) {
      q_error_data.emplace_back();
      v_error_data.emplace_back();
      e_error_data.emplace_back();
      for (int j = 0; j < integrations_per_integrator_; ++j) {
        q_error_data.back().emplace_back();
        v_error_data.back().emplace_back();
        e_error_data.back().emplace_back();
        for (int k = 0; k < tmax_.size(); ++k) {
          q_error_data.back().back().emplace_back(evaluations_[i][j][k],
                                                  q_errors_[i][j][k]);
          v_error_data.back().back().emplace_back(evaluations_[i][j][k],
                                                  v_errors_[i][j][k]);
          e_error_data.back().back().emplace_back(evaluations_[i][j][k],
                                                  e_errors_[i][j][k]);
        }
      }
      names.emplace_back(methods_[i].name);
      orders.emplace_back(methods_[i].orders);
      types.emplace_back(methods_[i].method_type);
    }
    std::string result;
    result += Set("qErrorData", q_error_data, ExpressInSIUnits);
    result += Set("vErrorData", v_error_data, ExpressInSIUnits);
    result += Set("eErrorData", e_error_data, ExpressInSIUnits);
    result += Set("names", names);
    result += Set("orders", orders);
    result += Set("types", types);
    return result;
  }

 private:
  absl::Status Integrate(int const method_index, int const time_step_index) {
    auto const& method = methods_[method_index];
    Problem problem;
    InitialValueProblem<FirstOrderODE> first_order_problem;
    int number_of_evaluations = 0;
    problem.equation.compute_acceleration = std::bind(
        compute_accelerations_, _1, _2, _3, &number_of_evaluations);
    problem.initial_state = initial_state_;
    first_order_problem.equation.compute_derivative =
        [&number_of_evaluations, this](
            Instant const& s,
            std::tuple<Position<World>, Velocity<World>> const& y,
            std::tuple<Velocity<World>, Vector<Acceleration, World>>& yʹ) {
          std::vector<Vector<Acceleration, World>> acceleration(1, {});
          auto const& [q, v] = y;
          auto const status = compute_accelerations_(
              s, {q}, acceleration, &number_of_evaluations);
          yʹ = {v, acceleration[0]};
          return status;
        };
    first_order_problem.initial_state.s = initial_state_.time;
    first_order_problem.initial_state.y = {initial_state_.positions[0],
                                           initial_state_.velocities[0]};
    auto const t0 = problem.initial_state.time.value;
    Length max_q_error;
    Speed max_v_error;
    Energy max_e_error;
    auto append_state = [this, &max_q_error, &max_v_error, &max_e_error](
        SecondOrderODE::State const& state) {
      auto const errors = compute_errors_(state);
      max_q_error = std::max(max_q_error, errors.q_error);
      max_v_error = std::max(max_v_error, errors.v_error);
      max_e_error = std::max(max_e_error, errors.e_error);
    };
    auto append_first_order_state =
    [this, &max_q_error, &max_v_error, &max_e_error](
        FirstOrderODE::State const& state) {
      SecondOrderODE::State second_order_state;
      second_order_state.time = state.s;
      auto const& [q, v] = state.y;
      second_order_state.positions = {q};
      second_order_state.velocities = {v};
      auto const errors = compute_errors_(second_order_state);
      max_q_error = std::max(max_q_error, errors.q_error);
      max_v_error = std::max(max_v_error, errors.v_error);
      max_e_error = std::max(max_e_error, errors.e_error);
    };
    switch (method.integrator.index()) {
      case 0: {
        FixedStepSizeIntegrator<SecondOrderODE> const& integrator =
            *std::get<0>(method.integrator);
        Time const Δt = method.evaluations *
                        starting_step_size_per_evaluation_ /
                        std::pow(step_reduction_, time_step_index);
        auto const instance = integrator.NewInstance(problem, append_state, Δt);
        SolveWithFixedStep(method_index,
                           time_step_index,
                           *instance,
                           method,
                           number_of_evaluations,
                           max_q_error,
                           max_v_error,
                           max_e_error,
                           t0,
                           Δt);
        break;
      }
      case 1: {
        AdaptiveStepSizeIntegrator<SecondOrderODE> const& integrator =
            *std::get<1>(method.integrator);
        Length const tolerance =
            std::exp2(-50.0 * time_step_index / integrations_per_integrator_) *
            first_tolerance_;
        auto const instance = integrator.NewInstance(
            problem,
            append_state,
            /*tolerance_to_error_ratio=*/
            [time_step_index, tolerance, this](
                Time const& current_step_size,
                SecondOrderODE::State const& state,
                SecondOrderODE::State::Error const& error) {
              return tolerance / error.position_error[0].Norm();
            },
            AdaptiveStepSizeIntegrator<SecondOrderODE>::Parameters{
                /*first_step=*/tmax_[0] - t0,
                /*safety_factor=*/0.9,
                /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
                /*last_step_is_exact=*/false});
        SolveWithAdaptiveStep(method_index,
                              time_step_index,
                              *instance,
                              method,
                              number_of_evaluations,
                              max_q_error,
                              max_v_error,
                              max_e_error);
        break;
      }
      case 2: {
        FixedStepSizeIntegrator<FirstOrderODE> const& integrator =
            *std::get<2>(method.integrator);
        Time const Δt = method.evaluations *
                        starting_step_size_per_evaluation_ /
                        std::pow(step_reduction_, time_step_index);
        auto const instance = integrator.NewInstance(
            first_order_problem, append_first_order_state, Δt);
        SolveWithFixedStep(method_index,
                           time_step_index,
                           *instance,
                           method,
                           number_of_evaluations,
                           max_q_error,
                           max_v_error,
                           max_e_error,
                           t0,
                           Δt);
        break;
      }
      case 3: {
        AdaptiveStepSizeIntegrator<FirstOrderODE> const& integrator =
            *std::get<3>(method.integrator);
        Length const tolerance =
            std::exp2(-50.0 * time_step_index / integrations_per_integrator_) *
            first_tolerance_;
        auto const instance = integrator.NewInstance(
            first_order_problem,
            append_first_order_state,
            /*tolerance_to_error_ratio=*/
            [time_step_index, tolerance, this](
                Time const& current_step_size,
                FirstOrderODE::State const& state,
                FirstOrderODE::State::Error const& error) {
              return tolerance / std::get<0>(error).Norm();
            },
            AdaptiveStepSizeIntegrator<FirstOrderODE>::Parameters{
                /*first_step=*/tmax_[0] - t0,
                /*safety_factor=*/0.9,
                /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
                /*last_step_is_exact=*/false});
        SolveWithAdaptiveStep(method_index,
                              time_step_index,
                              *instance,
                              method,
                              number_of_evaluations,
                              max_q_error,
                              max_v_error,
                              max_e_error);
        break;
      }
      default:
        LOG(FATAL) << "Unexpected variant " << method.integrator.index();
    }

    return absl::OkStatus();
  }

  template<typename Instance>
  void SolveWithFixedStep(int const method_index,
                          int const time_step_index,
                          Instance& instance,
                          PlottedIntegrator const& method,
                          int const& number_of_evaluations,
                          Length const& max_q_error,
                          Speed const& max_v_error,
                          Energy const& max_e_error,
                          Instant const& t0,
                          Time const& Δt) {
    for (auto const& [i, tmax] : std::ranges::enumerate_view(tmax_)) {
      CHECK_OK(instance.Solve(tmax));
      // Log both the actual number of evaluations and a theoretical number
      // that ignores any startup costs; that theoretical number is the one
      // used for plotting.
      int const amortized_evaluations =
          method.evaluations * static_cast<int>(std::floor((tmax - t0) / Δt));
      LOG_EVERY_N(INFO, 50)
          << "[" << method_index << "," << time_step_index << "] "
          << problem_name_ << ": " << number_of_evaluations
          << " actual evaluations (" << amortized_evaluations
          << " amortized) with " << method.name;
      // We plot the maximum error, i.e., the L∞ norm of the error.
      // [BM02] or [BCR01a] tend to use the average error (the normalized L¹
      // norm) instead.
      q_errors_[method_index][time_step_index][i] = max_q_error;
      v_errors_[method_index][time_step_index][i] = max_v_error;
      e_errors_[method_index][time_step_index][i] = max_e_error;
      evaluations_[method_index][time_step_index][i] = amortized_evaluations;
    }
  }

  template<typename Instance>
  void SolveWithAdaptiveStep(int const method_index,
                             int const time_step_index,
                             Instance& instance,
                             PlottedIntegrator const& method,
                             int const& number_of_evaluations,
                             Length const& max_q_error,
                             Speed const& max_v_error,
                             Energy const& max_e_error) {
    for (auto const& [i, tmax] : std::ranges::enumerate_view(tmax_)) {
      CHECK_OK(instance.Solve(tmax));
      LOG_EVERY_N(INFO, 50)
          << "[" << method_index << "," << time_step_index << "] "
          << problem_name_ << ": " << number_of_evaluations
          << " actual evaluations with " << method.name;
      // We plot the maximum error, i.e., the L∞ norm of the error.
      // [BM02] or [BCR01a] tend to use the average error (the normalized
      // L¹ norm) instead.
      q_errors_[method_index][time_step_index][i] = max_q_error;
      v_errors_[method_index][time_step_index][i] = max_v_error;
      e_errors_[method_index][time_step_index][i] = max_e_error;
      evaluations_[method_index][time_step_index][i] = number_of_evaluations;
    }
  }

  std::vector<PlottedIntegrator> const methods_;
  std::function<absl::Status(Instant const& t,
                             std::vector<Position<World>> const& q,
                             std::vector<Vector<Acceleration, World>>& result,
                             int* evaluations)>
      compute_accelerations_;
  SecondOrderODE::State initial_state_;
  std::function<Errors(SecondOrderODE::State const&)> compute_errors_;
  std::vector<std::vector<std::vector<Length>>> q_errors_;
  std::vector<std::vector<std::vector<Speed>>> v_errors_;
  std::vector<std::vector<std::vector<Energy>>> e_errors_;
  std::vector<std::vector<std::vector<double>>> evaluations_;
  std::vector<Instant> const tmax_;
  Length const first_tolerance_;
  std::string const problem_name_;
  Time const starting_step_size_per_evaluation_ = 1 * Second;
  int const integrations_per_integrator_ = 500;
  double const step_reduction_;
};


void GenerateSimpleHarmonicMotionWorkErrorGraphs() {
  SecondOrderODE::State initial_state;
  Instant const t0;
  Length const q_amplitude = 1 * Metre;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Stiffness const k = si::Unit<Stiffness>;
  Mass const m = 1 * Kilogram;

  initial_state.positions.emplace_back(
      World::origin + Displacement<World>({q_amplitude, 0 * Metre, 0 * Metre}));
  initial_state.velocities.emplace_back(World::unmoving);
  initial_state.time = DoublePrecision<Instant>(t0);

  std::vector<Instant> const tmax = {
      t0 + 50 * Second, t0 + 57.5 * Second, t0 + 75 * Second};
  auto const compute_error = [q_amplitude, v_amplitude, ω, m, k, t0](
      SecondOrderODE::State const& state) {
    return WorkErrorGraphGenerator<Energy>::Errors{
        AbsoluteError(
            q_amplitude * Vector<double, World>(
                              {Cos(ω * (state.time.value - t0)), 0, 0}) +
                World::origin,
            state.positions[0].value),
        AbsoluteError(
            -v_amplitude *
                Vector<double, World>({Sin(ω * (state.time.value - t0)), 0, 0}),
            state.velocities[0].value),
        AbsoluteError(0.5 * Joule,
                      (m * state.velocities[0].value.Norm²() +
                       k * (state.positions[0].value - World::origin).Norm²()) /
                          2)};
  };
  WorkErrorGraphGenerator<Energy> generator(
      ComputeHarmonicOscillatorAcceleration3D<World>,
      initial_state,
      compute_error,
      tmax,
      1 * Metre,
      // 585 μs yields a step reduction factor of 1.015, which is what the old
      // graphs used.
      585 * Micro(Second),
      "Harmonic oscillator");

  OFStream file(TEMP_DIR / "simple_harmonic_motion_graphs.generated.wl");
  file << generator.GetMathematicaData();
}

void GenerateKeplerProblemWorkErrorGraphs(double const eccentricity) {
  SecondOrderODE::State initial_state;
  Instant const t0;
  GravitationalParameter const μ = si::Unit<GravitationalParameter>;
  MassiveBody b1(μ);
  MasslessBody b2;

  KeplerianElements<World> elements;
  elements.semimajor_axis = 1 * Metre;
  elements.eccentricity = eccentricity;
  elements.argument_of_periapsis = 0 * Degree;
  elements.true_anomaly = 0 * Degree;
  KeplerOrbit<World> const orbit(b1, b2, elements, t0);

  auto const initial_dof = orbit.StateVectors(t0);
  CHECK_EQ(initial_dof.displacement().coordinates().z, 0 * Metre);
  CHECK_EQ(initial_dof.velocity().coordinates().z, 0 * Metre / Second);

  initial_state.positions.emplace_back(World::origin +
                                       initial_dof.displacement());
  initial_state.velocities.emplace_back(initial_dof.velocity());
  initial_state.time = DoublePrecision<Instant>(t0);

  std::vector<Instant> const tmax = {
      t0 + 8 * *orbit.elements_at_epoch().period,
      t0 + 10 * *orbit.elements_at_epoch().period,
      t0 + 12 * *orbit.elements_at_epoch().period};

  SpecificEnergy const initial_specific_energy =
      initial_dof.velocity().Norm²() / 2 -
      μ / initial_dof.displacement().Norm();

  auto const compute_error = [&orbit, μ, initial_specific_energy](
      SecondOrderODE::State const& state) {
    Displacement<World> r = state.positions[0].value - World::origin;
    Velocity<World> v = state.velocities[0].value;
    auto const expected_dof = orbit.StateVectors(state.time.value);
    return WorkErrorGraphGenerator<SpecificEnergy>::Errors{
        AbsoluteError(expected_dof.displacement(), r),
        AbsoluteError(expected_dof.velocity(), v),
        AbsoluteError(initial_specific_energy, v.Norm²() / 2 - μ / r.Norm())};
  };

  WorkErrorGraphGenerator<SpecificEnergy> generator(
      ComputeKeplerAcceleration<World>,
      initial_state,
      compute_error,
      tmax,
      1 * Metre,
      585 * Micro(Second) * Sqrt((1 - eccentricity) / (1 + eccentricity)),
      " Kepler problem with e = " + std::to_string(eccentricity));

  OFStream file(TEMP_DIR / ("kepler_problem_graphs_" +
                            std::to_string(eccentricity) + ".generated.wl"));
  file << generator.GetMathematicaData();
}

}  // namespace internal
}  // namespace _integrator_plots
}  // namespace mathematica
}  // namespace principia
