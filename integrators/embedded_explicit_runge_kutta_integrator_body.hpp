
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

  // State increment (high-order).
  StateDifference Δŷ;
  // Current state.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  auto& ŷ = current_state.y;

  // Difference between the low- and high-order approximations.
  typename ODE::SystemStateError error_estimate;

  // Current Runge-Kutta stage.
  State y_stage;

  StateVariation f;
  std::vector<StateDifference> k(stages_);
  for (auto& k_stage : k) {
    for_all_of(ŷ, k_stage).loop([](auto const& ŷ, auto& k_stage) {
      int const dimension = ŷ.size();
      k_stage.resize(dimension);
    });
  }

  for_all_of(ŷ, error_estimate, f, y_stage)
      .loop([](auto const& ŷ, auto& error_estimate, auto& f, auto& y_stage) {
        int const dimension = ŷ.size();
        error_estimate.resize(dimension);
        f.resize(dimension);
        y_stage.resize(dimension);
      });

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

      // Runge-Kutta iteration; fills |k|.
      for (int i = first_stage; i < stages_; ++i) {
        Instant const t_stage =
            (parameters.last_step_is_exact && at_end && c[i] == 1.0)
                ? t_final
                : t.value + (t.error + c[i] * h);

        StateDifference Σj_a_ij_k_j{};
        for (int j = 0; j < i; ++j) {
          for_all_of(k[j], ŷ, y_stage, Σj_a_ij_k_j)
              .loop([](auto const& k_j,
                       auto const& ŷ,
                       auto& y_stage,
                       auto& Σj_a_ij_k_j) {
                int const dimension = ŷ.size();
                Σj_a_ij_k_j.resize(dimension);
                for (int l = 0; l < dimension; ++l) {
                  Σj_a_ij_k_j[l] += a(i, j) * k_j[l];
                  y_stage[l] = ŷ[l].value + Σj_a_ij_k_j[l];
                }
              });
        }

        step_status.Update(equation.compute_derivative(t_stage, y_stage, f));
        for_all_of(f, k[i]).loop([](auto const& f, auto& k_i) {
          int const dimension = f.size();
          for (int l = 0; l < dimension; ++l) {
            k_i[l] = h * f[l];
          }
        });
      }

      // Increment computation and step size control.
      StateDifference Σi_b̂_i_k_i{};
      StateDifference Σi_b_i_k_i{};
      for (int i = 0; i < stages_; ++i) {
        for_all_of(k[i], ŷ, Δŷ, Σi_b̂_i_k_i, Σi_b_i_k_i, error_estimate)
            .loop([](auto const& k_i,
                     auto const& ŷ,
                     auto& Δŷ,
                     auto& Σi_b̂_i_k_i,
                     auto& Σi_b_i_k_i,
                     auto& error_estimate) {
              int const dimension = ŷ.size();
              Σi_b̂_i_k_i.resize(dimension);
              Σi_b_i_k_i.resize(dimension);
              for (int l = 0; l < dimension; ++l) {
                Σi_b̂_i_k_i[l] += b̂[i] + k_i[l];
                Σi_b_i_k_i[l] += b[i] + k_i[l];
                Δŷ[l] = Σi_b̂_i_k_i[l];
                auto const Δy_l = Σi_b_i_k_i[l];
                error_estimate[l] = Δy_l - Δŷ[l];
              }
            });
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
    for_all_of(Δŷ, ŷ).loop([](auto const& Δŷ, auto& ŷ) {
      int const dimension = ŷ.size();
      for (int k = 0; k < dimension; ++k) {
        ŷ[k].Increment(Δŷ[k]);
      }
    });

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
