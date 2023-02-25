#pragma once

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <optional>
#include <vector>

#include "base/jthread.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_nyström_integrator {

using geometry::Sign;
using numerics::DoublePrecision;
using quantities::DebugString;
using quantities::Difference;
using quantities::Quotient;
using namespace principia::base::_not_null;

template<typename Method, typename ODE_>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
EmbeddedExplicitRungeKuttaNyströmIntegrator() {
  // The first node is always 0 in an explicit method.
  CHECK_EQ(0.0, c_[0]);
  if (first_same_as_last) {
    // Check that the conditions for the FSAL property are satisfied, see for
    // instance [DEP87a], equation 3.1.
    CHECK_EQ(1.0, c_[stages_ - 1]);
    CHECK_EQ(0.0, b̂_[stages_ - 1]);
    for (int j = 0; j < stages_ - 1; ++j) {
      CHECK_EQ(b̂_[j], a_(stages_ - 1, j));
    }
  }
}

template<typename Method, typename ODE_>
absl::Status EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::Solve(Instant const& t_final) {
  using Position = typename ODE::DependentVariable;
  using Displacement = typename ODE::DependentVariableDifference;
  using Velocity = typename ODE::DependentVariableDerivative;
  using Acceleration = typename ODE::DependentVariableDerivative2;

  auto const& a = integrator_.a_;
  auto const& b̂ = integrator_.b̂_;
  auto const& b̂ʹ = integrator_.b̂ʹ_;
  auto const& b = integrator_.b_;
  auto const& bʹ = integrator_.bʹ_;
  auto const& c = integrator_.c_;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto& first_use = this->first_use_;
  auto& parameters = this->parameters_;
  auto const& equation = this->equation_;

  // |current_state| gets updated as the integration progresses to allow
  // restartability.

  // State before the last, truncated step.
  std::optional<typename ODE::State> final_state;

  // Argument checks.
  int const dimension = current_state.positions.size();
  Sign const integration_direction = Sign(parameters.first_step);
  if (integration_direction.is_positive()) {
    // Integrating forward.
    CHECK_LT(current_state.time.value, t_final);
  } else {
    // Integrating backward.
    CHECK_GT(current_state.time.value, t_final);
  }
  CHECK(first_use || !parameters.last_step_is_exact)
      << "Cannot reuse an instance where the last step is exact";
  first_use = false;

  // Time step.  Updated as the integration progresses to allow restartability.
  Time& h = this->step_;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<Instant>& t = current_state.time;

  // Position increment (high-order).
  std::vector<Displacement> Δq̂(dimension);
  // Velocity increment (high-order).
  std::vector<Velocity> Δv̂(dimension);
  // Current position.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Position>>& q̂ = current_state.positions;
  // Current velocity.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Velocity>>& v̂ = current_state.velocities;

  // Difference between the low- and high-order approximations.
  typename ODE::State::Error error_estimate;
  error_estimate.position_error.resize(dimension);
  error_estimate.velocity_error.resize(dimension);

  // Current Runge-Kutta-Nyström stage.
  std::vector<Position> q_stage(dimension);
  // Accelerations at each stage.
  // TODO(egg): this is a rectangular container, use something more appropriate.
  std::vector<std::vector<Acceleration>> g(stages_);
  for (auto& g_stage : g) {
    g_stage.resize(dimension);
  }

  bool at_end = false;
  double tolerance_to_error_ratio;

  // The first stage of the Runge-Kutta-Nyström iteration.  In the FSAL case,
  // |first_stage == 1| after the first step, since the first RHS evaluation has
  // already occurred in the previous step.  In the non-FSAL case and in the
  // first step of the FSAL case, |first_stage == 0|.
  int first_stage = 0;

  // The number of steps already performed.
  std::int64_t step_count = 0;

  absl::Status status;
  absl::Status step_status;

  // No step size control on the first step.  If this instance is being
  // restarted we already have a value of |h| suitable for the next step, based
  // on the computation of |tolerance_to_error_ratio_| during the last
  // invocation.
  goto runge_kutta_nyström_step;

  while (!at_end) {
    // Compute the next step with decreasing step sizes until the error is
    // tolerable.
    do {
      // Reset the status as any error returned by a force computation for a
      // rejected step is now moot.
      step_status = absl::OkStatus();

      // Adapt step size.
      // TODO(egg): find out whether there's a smarter way to compute that root,
      // especially since we make the order compile-time.
      h *= parameters.safety_factor *
               std::pow(tolerance_to_error_ratio, 1.0 / (lower_order + 1));
      // TODO(egg): should we check whether it vanishes in double precision
      // instead?
      if (t.value + (t.error + h) == t.value) {
        return absl::Status(termination_condition::VanishingStepSize,
                            "At time " + DebugString(t.value) +
                                ", step size is effectively zero.  "
                                "Singularity or stiff system suspected.");
      }

    runge_kutta_nyström_step:
      // Termination condition.
      if (parameters.last_step_is_exact) {
        Time const time_to_end = (t_final - t.value) - t.error;
        at_end = integration_direction * h >=
                 integration_direction * time_to_end;
        if (at_end) {
          // The chosen step size will overshoot.  Clip it to just reach the
          // end, and terminate if the step is accepted.  Note that while this
          // step size is a good approximation, there is no guarantee that it
          // won't over/undershoot, so we still need to special case the very
          // last stage below.
          h = time_to_end;
          final_state = current_state;
        }
      }

      auto const h² = h * h;

      // Runge-Kutta-Nyström iteration; fills |g|.
      for (int i = first_stage; i < stages_; ++i) {
        Instant const t_stage =
            (parameters.last_step_is_exact && at_end && c[i] == 1.0)
                ? t_final
                : t.value + (t.error + c[i] * h);
        for (int k = 0; k < dimension; ++k) {
          Acceleration Σⱼ_aᵢⱼ_gⱼₖ{};
          for (int j = 0; j < i; ++j) {
            Σⱼ_aᵢⱼ_gⱼₖ += a(i, j) * g[j][k];
          }
          q_stage[k] = q̂[k].value + h * c[i] * v̂[k].value + h² * Σⱼ_aᵢⱼ_gⱼₖ;
        }
        termination_condition::UpdateWithAbort(
            equation.compute_acceleration(t_stage, q_stage, g[i]), step_status);
      }

      // Increment computation and step size control.
      for (int k = 0; k < dimension; ++k) {
        Acceleration Σᵢ_b̂ᵢ_gᵢₖ{};
        Acceleration Σᵢ_bᵢ_gᵢₖ{};
        Acceleration Σᵢ_b̂ʹᵢ_gᵢₖ{};
        Acceleration Σᵢ_bʹᵢ_gᵢₖ{};
        // Please keep the eight assigments below aligned, they become illegible
        // otherwise.
        for (int i = 0; i < stages_; ++i) {
          Σᵢ_b̂ᵢ_gᵢₖ  += b̂[i] * g[i][k];
          Σᵢ_bᵢ_gᵢₖ  += b[i] * g[i][k];
          Σᵢ_b̂ʹᵢ_gᵢₖ += b̂ʹ[i] * g[i][k];
          Σᵢ_bʹᵢ_gᵢₖ += bʹ[i] * g[i][k];
        }
        // The hat-less Δq and Δv are the low-order increments.
        Δq̂[k]                  = h * v̂[k].value + h² * Σᵢ_b̂ᵢ_gᵢₖ;
        Displacement const Δqₖ = h * v̂[k].value + h² * Σᵢ_bᵢ_gᵢₖ;
        Δv̂[k]                  = h * Σᵢ_b̂ʹᵢ_gᵢₖ;
        Velocity const Δvₖ     = h * Σᵢ_bʹᵢ_gᵢₖ;

        error_estimate.position_error[k] = Δqₖ - Δq̂[k];
        error_estimate.velocity_error[k] = Δvₖ - Δv̂[k];
      }
      tolerance_to_error_ratio =
          this->tolerance_to_error_ratio_(h, current_state, error_estimate);
    } while (tolerance_to_error_ratio < 1.0);

    status.Update(step_status);

    if (!parameters.last_step_is_exact && t.value + (t.error + h) > t_final) {
      // We did overshoot.  Drop the point that we just computed and exit.
      final_state = current_state;
      break;
    }

    if (first_same_as_last) {
      using std::swap;
      swap(g.front(), g.back());
      first_stage = 1;
    }

    // Increment the solution with the high-order approximation.
    t.Increment(h);
    for (int k = 0; k < dimension; ++k) {
      q̂[k].Increment(Δq̂[k]);
      v̂[k].Increment(Δv̂[k]);
    }
    RETURN_IF_STOPPED;
    append_state(current_state);
    ++step_count;
    if (step_count == parameters.max_steps && !at_end) {
      return absl::Status(termination_condition::ReachedMaximalStepCount,
                          "Reached maximum step count " +
                              std::to_string(parameters.max_steps) +
                              " at time " + DebugString(t.value) +
                              "; requested t_final is " + DebugString(t_final) +
                              ".");
    }
  }
  // The resolution is restartable from the last non-truncated state.
  CHECK(final_state);
  current_state = *final_state;
  return status;
}

template<typename Method, typename ODE_>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_> const&
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
void EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  AdaptiveStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  [[maybe_unused]] auto* const extension =
      message
          ->MutableExtension(
              serialization::AdaptiveStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::
                  EmbeddedExplicitRungeKuttaNystromIntegratorInstance::
                      extension);
}

template<typename Method, typename ODE_>
template<typename, typename>
not_null<std::unique_ptr<
    typename EmbeddedExplicitRungeKuttaNyströmIntegrator<Method,
                                                         ODE_>::Instance>>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::Instance::
ReadFromMessage(
    serialization::
        EmbeddedExplicitRungeKuttaNystromIntegratorInstance const&
            extension,
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters,
    Time const& time_step,
    bool const first_use,
    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator) {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(new Instance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters,
                                                time_step,
                                                first_use,
                                                integrator));
}

template<typename Method, typename ODE_>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::Instance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters,
    Time const& time_step,
    bool const first_use,
    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator)
    : AdaptiveStepSizeIntegrator<ODE>::Instance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters,
                                                time_step,
                                                first_use),
      integrator_(integrator) {}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
NewInstance(InitialValueProblem<ODE> const& problem,
            AppendState const& append_state,
            ToleranceToErrorRatio const& tolerance_to_error_ratio,
            Parameters const& parameters) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem,
                   append_state,
                   tolerance_to_error_ratio,
                   parameters,
                   /*step=*/parameters.first_step,
                   /*first_use=*/true,
                   *this));
}

template<typename Method, typename ODE_>
void EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_>::
WriteToMessage(not_null<serialization::AdaptiveStepSizeIntegrator*> message)
    const {
  message->set_kind(Method::kind);
}

}  // namespace internal_embedded_explicit_runge_kutta_nyström_integrator

template<typename Method, typename ODE_>
internal_embedded_explicit_runge_kutta_nyström_integrator::
    EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_> const&
EmbeddedExplicitRungeKuttaNyströmIntegrator() {
  static_assert(
      std::is_base_of<methods::EmbeddedExplicitRungeKuttaNyström,
                      Method>::value,
      "Method must be derived from EmbeddedExplicitRungeKuttaNyström");
  static internal_embedded_explicit_runge_kutta_nyström_integrator::
      EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE_> const
          integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
