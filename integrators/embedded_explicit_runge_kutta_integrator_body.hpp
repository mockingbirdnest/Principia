
#pragma once

#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <optional>
#include <vector>

#include "base/for_all_of.hpp"
#include "base/jthread.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_integrator {

using base::for_all_of;
using base::make_not_null_unique;
using geometry::Sign;
using numerics::DoublePrecision;
using quantities::DebugString;
using quantities::Difference;
using quantities::Quotient;

template<typename Method, typename... StateElements>
EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
EmbeddedExplicitRungeKuttaIntegrator() {
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

template<typename Method, typename... StateElements>
absl::Status EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
Instance::Solve(Instant const& t_final) {
  using State = typename ODE::State;
  using StateDifference = typename ODE::StateDifference;
  using StateVariation = typename ODE::StateVariation;

  auto const& a = integrator_.a_;
  auto const& b̂ = integrator_.b̂_;
  auto const& b = integrator_.b_;
  auto const& c = integrator_.c_;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto& first_use = this->first_use_;
  auto& parameters = this->parameters_;
  auto const& equation = this->equation_;

  // |current_state| gets updated as the integration progresses to allow
  // restartability.

  // State before the last, truncated step.
  std::optional<typename ODE::SystemState> final_state;

  // Argument checks.
  int const dimension = current_state.y.size();
  Sign const integration_direction = Sign(parameters.first_time_step);
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
  Time& h = this->time_step_;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<Instant>& t = current_state.time;

  // Position increment (high-order).
  std::vector<Displacement> Δŷ(dimension);
  // Current state.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  auto& ŷ = current_state.y;

  // Difference between the low- and high-order approximations.
  typename ODE::SystemStateError error_estimate;
  error_estimate.position_error.resize(dimension);

  // Current Runge-Kutta stage.
  State y_stage(dimension);
  // Variations at each stage.
  std::vector<StateVariation> f(stages_);
  for (auto& f_stage : f) {
    f_stage.resize(dimension);
  }

  bool at_end = false;
  double tolerance_to_error_ratio;

  // The first stage of the Runge-Kutta iteration.  In the FSAL case,
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
  goto runge_kutta_step;

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

    runge_kutta_step:
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

      // Runge-Kutta iteration; fills |f|.
      for (int i = first_stage; i < stages_; ++i) {
        Instant const t_stage =
            (parameters.last_step_is_exact && at_end && c[i] == 1.0)
                ? t_final
                : t.value + (t.error + c[i] * h);
        for_all_of(
          [](){
          for (int k = 0; k < dimension; ++k) {
            Acceleration Σj_a_ij_f_jk{};
            for (int j = 0; j < i; ++j) {
              Σj_a_ij_f_jk += a(i, j) * f[j][k];
            }
            y_stage[k] = q̂[k].value + h * c[i] * v̂[k].value + h² * Σj_a_ij_f_jk;
          }
          step_status.Update(
              equation.compute_acceleration(t_stage, y_stage, f[i]));
          }
      });

      // Increment computation and step size control.
      for (int k = 0; k < dimension; ++k) {
        Acceleration Σi_b̂_i_f_ik{};
        Acceleration Σi_b_i_f_ik{};
        // Please keep the eight assigments below aligned, they become illegible
        // otherwise.
        for (int i = 0; i < stages_; ++i) {
          Σi_b̂_i_f_ik  += b̂[i] * f[i][k];
          Σi_b_i_f_ik  += b[i] * f[i][k];
        }
        // The hat-less Δq is the low-order increments.
        Δq̂[k]                   = h * v̂[k].value + h² * Σi_b̂_i_f_ik;
        Displacement const Δq_k = h * v̂[k].value + h² * Σi_b_i_f_ik;

        error_estimate.position_error[k] = Δq_k - Δq̂[k];
      }
      tolerance_to_error_ratio =
          this->tolerance_to_error_ratio_(h, error_estimate);
    } while (tolerance_to_error_ratio < 1.0);

    status.Update(step_status);

    if (!parameters.last_step_is_exact && t.value + (t.error + h) > t_final) {
      // We did overshoot.  Drop the point that we just computed and exit.
      final_state = current_state;
      break;
    }

    if (first_same_as_last) {
      using std::swap;
      swap(f.front(), f.back());
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

template<typename Method, typename... StateElements>
EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...> const&
EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename... StateElements>
not_null<std::unique_ptr<
    typename Integrator<ExplicitFirstOrderOrdinaryDifferentialEquation<
        StateElements...>>::Instance>>
EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename... StateElements>
void EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
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

#if 0
template<typename Method, typename... StateElements>
template<typename, typename>
not_null<std::unique_ptr<
    typename EmbeddedExplicitRungeKuttaIntegrator<Method,
                                                  StateElements...>::Instance>>
EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::Instance::
ReadFromMessage(
    serialization::
        EmbeddedExplicitRungeKuttaNystromIntegratorInstance const&
            extension,
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters,
    Time const& time_step,
    bool const first_use,
    EmbeddedExplicitRungeKuttaIntegrator const& integrator) {
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
#endif

template<typename Method, typename... StateElements>
EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters,
    Time const& time_step,
    bool const first_use,
    EmbeddedExplicitRungeKuttaIntegrator const& integrator)
    : AdaptiveStepSizeIntegrator<ODE>::Instance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters,
                                                time_step,
                                                first_use),
      integrator_(integrator) {}

template<typename Method, typename... StateElements>
not_null<std::unique_ptr<typename Integrator<
    ExplicitFirstOrderOrdinaryDifferentialEquation<StateElements...>>::Instance>>
EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
NewInstance(IntegrationProblem<ODE> const& problem,
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
                   /*time_step=*/parameters.first_time_step,
                   /*first_use=*/true,
                   *this));
}

template<typename Method, typename... StateElements>
void EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...>::
WriteToMessage(not_null<serialization::AdaptiveStepSizeIntegrator*> message)
    const {
  message->set_kind(Method::kind);
}

}  // namespace internal_embedded_explicit_runge_kutta_integrator

template<typename Method, typename... StateElements>
internal_embedded_explicit_runge_kutta_integrator::
    EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...> const&
EmbeddedExplicitRungeKuttaIntegrator() {
  static_assert(
      std::is_base_of<methods::EmbeddedExplicitRungeKutta,
                      Method>::value,
      "Method must be derived from EmbeddedExplicitRungeKutta");
  static internal_embedded_explicit_runge_kutta_integrator::
      EmbeddedExplicitRungeKuttaIntegrator<Method, StateElements...> const
          integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
