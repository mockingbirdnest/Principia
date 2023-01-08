#pragma once

#include "integrators/explicit_runge_kutta_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <optional>
#include <vector>

#include "base/for_all_of.hpp"
#include "base/jthread.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_explicit_runge_kutta_integrator {

using base::for_all_of;
using base::make_not_null_unique;
using geometry::Sign;
using numerics::DoublePrecision;
using quantities::Abs;
using quantities::DebugString;
using quantities::Difference;
using quantities::Quotient;

template<typename Method, typename ODE_>
ExplicitRungeKuttaIntegrator<Method, ODE_>::
ExplicitRungeKuttaIntegrator() {
  // The first node is always 0 in an explicit method.
  CHECK_EQ(0.0, c_[0]);
  if (first_same_as_last) {
    // Check that the conditions for the FSAL property are satisfied, see for
    // instance [DEP87a], equation 3.1.
    CHECK_EQ(1.0, c_[stages_ - 1]);
    CHECK_EQ(0.0, b_[stages_ - 1]);
    for (int j = 0; j < stages_ - 1; ++j) {
      CHECK_EQ(b_[j], a_(stages_ - 1, j));
    }
  }
}

template<typename Method, typename ODE_>
absl::Status
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Solve(typename ODE::IndependentVariable const& s_final) {
  using IndependentVariable = typename ODE::IndependentVariable;
  using IndependentVariableDifference =
      typename ODE::IndependentVariableDifference;
  using DependentVariables = typename ODE::DependentVariables;
  using DependentVariableDifferences =
      typename ODE::DependentVariableDifferences;
  using DependentVariableDerivatives =
      typename ODE::DependentVariableDerivatives;

  auto const& a = integrator_.a_;
  auto const& b = integrator_.b_;
  auto const& c = integrator_.c_;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto const& equation = this->equation_;
  auto const& step = this->step_;

  // Argument checks.
  Sign const integration_direction = Sign(step);
  if (integration_direction.is_positive()) {
    // Integrating forward.
    CHECK_LT(current_state.s.value, s_final);
  } else {
    // Integrating backward.
    CHECK_GT(current_state.s.value, s_final);
  }

  // Step.
  IndependentVariableDifference const& h = step;
  IndependentVariableDifference const abs_h = integration_direction * h;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<IndependentVariable>& s = current_state.s;

  // DependentVariables increment (high-order).
  DependentVariableDifferences Δy;
  // Current state.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  auto& y = current_state.y;

  // Current Runge-Kutta stage.
  DependentVariables y_stage;

  DependentVariableDerivatives f;
  DependentVariableDerivatives last_f;
  std::vector<DependentVariableDifferences> k(stages_);
  for (auto& k_stage : k) {
    for_all_of(y, Δy, k_stage).loop([](auto const& y, auto& Δy, auto& k_stage) {
      int const dimension = y.size();
      Δy.resize(dimension);
      k_stage.resize(dimension);
    });
  }

  for_all_of(y, f, last_f, y_stage)
      .loop([](auto const& y,
               auto& f,
               auto& last_f,
               auto& y_stage) {
        int const dimension = y.size();
        f.resize(dimension);
        last_f.resize(dimension);
        for (auto const& yₗ : y) {
          y_stage.push_back(yₗ.value);
        }
      });

  absl::Status status;

  if (first_same_as_last) {
    status = equation.compute_derivative(s.value, y_stage, last_f);
  }

  while (abs_h <= Abs((s_final - s.value) - s.error)) {
      // Runge-Kutta iteration; fills |k|.
      for (int i = 0; i < stages_; ++i) {
        if (i == 0 && first_same_as_last) {
          // TODO(phl): Use pointers to avoid copying big objects.
          f = last_f;
        } else {
          DependentVariableDifferences Σⱼ_aᵢⱼ_kⱼ{};
          for (int j = 0; j < i; ++j) {
            for_all_of(k[j], y, y_stage, Σⱼ_aᵢⱼ_kⱼ)
                .loop([&a, i, j](auto const& kⱼ,
                                 auto const& y,
                                 auto& y_stage,
                                 auto& Σⱼ_aᵢⱼ_kⱼ) {
                  int const dimension = y.size();
                  Σⱼ_aᵢⱼ_kⱼ.resize(dimension);
                  for (int l = 0; l < dimension; ++l) {
                    Σⱼ_aᵢⱼ_kⱼ[l] += a(i, j) * kⱼ[l];
                  }
                });
          }
          for_all_of(y, Σⱼ_aᵢⱼ_kⱼ, y_stage)
              .loop([](auto const& y, auto const& Σⱼ_aᵢⱼ_kⱼ, auto& y_stage) {
                int const dimension = y.size();
                for (int l = 0; l < dimension; ++l) {
                  y_stage[l] = y[l].value + Σⱼ_aᵢⱼ_kⱼ[l];
                }
              });

          termination_condition::UpdateWithAbort(
              equation.compute_derivative(
                  s.value + (s.error + c[i] * h), y_stage, f),
              status);
        }
        for_all_of(f, k[i]).loop([h](auto const& f, auto& kᵢ) {
          int const dimension = f.size();
          for (int l = 0; l < dimension; ++l) {
            kᵢ[l] = h * f[l];
          }
        });
      }

      // Increment computation.
      DependentVariableDifferences Σᵢ_bᵢ_kᵢ{};
      for (int i = 0; i < stages_; ++i) {
        for_all_of(k[i], y, Δy, Σᵢ_bᵢ_kᵢ)
            .loop([&a, &b, i](auto const& kᵢ,
                              auto const& y,
                              auto& Δy,
                              auto& Σᵢ_bᵢ_kᵢ) {
              int const dimension = y.size();
              Σᵢ_bᵢ_kᵢ.resize(dimension);
              for (int l = 0; l < dimension; ++l) {
                Σᵢ_bᵢ_kᵢ[l] += b[i] * kᵢ[l];
                Δy[l] = Σᵢ_bᵢ_kᵢ[l];
              }
            });
      }

    if (first_same_as_last) {
      last_f = f;
    }

    // Increment the solution.
    s.Increment(h);
    for_all_of(Δy, y).loop([](auto const& Δy, auto& y) {
      int const dimension = y.size();
      for (int l = 0; l < dimension; ++l) {
        y[l].Increment(Δy[l]);
      }
    });

    RETURN_IF_STOPPED;
    append_state(current_state);
    if (absl::IsAborted(status)) {
      return status;
    }
  }

  return status;
}

template<typename Method, typename ODE_>
ExplicitRungeKuttaIntegrator<Method, ODE_> const&
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
integrator() const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
void ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
WriteToMessage(not_null<serialization::IntegratorInstance*> message) const {
  FixedStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  [[maybe_unused]] auto* const extension =
      message
          ->MutableExtension(
              serialization::FixedStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::
                  ExplicitRungeKuttaNystromIntegratorInstance::
                      extension);
}

#if 0
template<typename Method, typename ODE_>
template<typename, typename>
not_null<std::unique_ptr<
    typename ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance>>
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
ReadFromMessage(serialization::
                    ExplicitRungeKuttaNystromIntegratorInstance const&
                        extension,
                InitialValueProblem<ODE> const& problem,
                AppendState const& append_state,
                ToleranceToErrorRatio const& tolerance_to_error_ratio,
                Parameters const& parameters,
                Time const& time_step,
                bool const first_use,
                ExplicitRungeKuttaIntegrator const& integrator) {
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

template<typename Method, typename ODE_>
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Instance(InitialValueProblem<ODE> const& problem,
         AppendState const& append_state,
         typename ODE::IndependentVariableDifference const& step,
         ExplicitRungeKuttaIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem,
                                             append_state,
                                             step),
      integrator_(integrator) {}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ExplicitRungeKuttaIntegrator<Method, ODE_>::
NewInstance(InitialValueProblem<ODE> const& problem,
            AppendState const& append_state,
            typename ODE::IndependentVariableDifference const& step) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename ODE_>
void ExplicitRungeKuttaIntegrator<Method, ODE_>::
WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal_explicit_runge_kutta_integrator

template<typename Method, typename ODE_>
internal_explicit_runge_kutta_integrator::
    ExplicitRungeKuttaIntegrator<Method, ODE_> const&
ExplicitRungeKuttaIntegrator() {
  static_assert(
      std::is_base_of<methods::ExplicitRungeKutta,
                      Method>::value,
      "Method must be derived from ExplicitRungeKutta");
  static internal_explicit_runge_kutta_integrator::
      ExplicitRungeKuttaIntegrator<Method, ODE_> const integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
