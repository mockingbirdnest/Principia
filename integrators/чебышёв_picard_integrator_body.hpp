#pragma once

#include "integrators/чебышёв_picard_integrator.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>

#include "base/for_all_of.hpp"
#include "base/status_utilities.hpp"  // 🧙 For RETURN_IF_ERROR.
#include "base/tags.hpp"
#include "geometry/sign.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/matrix_views.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace integrators {
namespace _чебышёв_picard_integrator {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::base::_tags;
using namespace principia::geometry::_sign;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_matrix_views;
using namespace principia::quantities::_si;

// Strip DoublePrecision from a tuple.
template<typename... T>
DirectSum<T...> StripDoublePrecision(
    DirectSum<DoublePrecision<T>...> const& in) {
  DirectSum<T...> out;
  for_all_of(in, out).loop(
      [](auto const& inᵢ, auto& outᵢ) { outᵢ = inᵢ.value; });
  return out;
}

// Wrap a tuple in DoublePrecision.
template<typename... T>
DirectSum<DoublePrecision<T>...> WrapInDoublePrecision(
    DirectSum<T...> const& in) {
  DirectSum<DoublePrecision<T>...> out;
  for_all_of(in, out).loop([](auto const& inᵢ, auto& outᵢ) { outᵢ = inᵢ; });
  return out;
}

// Constructs a FixedVector of Чебышёв nodes of the second kind.
template<std::int64_t M>
FixedVector<double, M + 1> Nodes() {
  FixedVector<double, M + 1> nodes(uninitialized);
  for (std::int64_t i = 0; i <= M; ++i) {
    nodes[i] = -Cos(π / M * i * Radian);
  }
  return nodes;
}

// Constructs an N×N diagonal matrix defined by a function.
template<std::int64_t N>
FixedMatrix<double, N, N, /*use_heap=*/true> DiagonalMatrix(
    std::function<double(std::int64_t)> const f) {
  FixedMatrix<double, N, N, /*use_heap=*/true> A;
  for (std::int64_t i = 0; i < N; ++i) {
    A(i, i) = f(i);
  }
  return A;
}

// Constructs an S matrix of size (N + 1)×N.
// See equation 1.26 in [Mac15].
template<std::int64_t N>
FixedMatrix<double, N + 1, N, /*use_heap=*/true> S() {
  FixedMatrix<double, N + 1, N, /*use_heap=*/true> S;
  S(0, 0) = 1;
  S(0, 1) = -0.5;
  for (std::int64_t k = 2; k < N; ++k) {
    S(0, k) = (k % 2 == 1 ? 1 : -1) * (1.0 / (k - 1) - 1.0 / (k + 1));
  }
  for (std::int64_t i = 0; i < N; ++i) {
    S(i + 1, i) = 1;
  }
  for (std::int64_t i = 1; i + 2 < N; ++i) {
    S(i, i + 1) = -1;
  }

  return S;
}

// V is is a diagonal (M + 1)×(M + 1) matrix with diagonal [1/M, 2/M, 2/M,
// ..., 1/M].
template<std::int64_t M>
FixedMatrix<double, M + 1, M + 1, /*use_heap=*/true> V() {
  constexpr double one_over_M = 1.0 / M;
  return DiagonalMatrix<M + 1>([](std::int64_t const i) {
    return (i == 0 || i == M) ? one_over_M : 2 * one_over_M;
  });
}

template<ЧебышёвPicardMethod Method>
ЧебышёвPicardMatrices<Method, 1>::ЧебышёвPicardMatrices(
    FixedVector<double, Method::M + 1> const& nodes)
    : CₓCα(uninitialized) {
  // We use the notation from [Mac15], section 1.4.3.

  // ᵝT is a (M + 1)×(N + 1) matrix of Чебышёв polynomials evaluated at nodes.
  // See [Mac15], equation (1.20).
  FixedMatrix<double, M + 1, N + 1, /*use_heap=*/true> ᵝT(uninitialized);

  for (std::int64_t i = 0; i <= M; ++i) {
    auto const τᵢ = nodes[i];
    // The 0-degree polynomial is uniformly 1.
    ᵝT(i, 0) = 1;
    // The 0-degree polynomial is the identity.
    ᵝT(i, 1) = τᵢ;

    // We populate the rest of ᵝT using the recurrence relation.
    for (std::int64_t j = 2; j <= N; ++j) {
      ᵝT(i, j) = 2 * τᵢ * ᵝT(i, j - 1) - ᵝT(i, j - 2);
    }
  }

  // ᵝW is a diagonal (N + 1)×(N + 1) matrix with diagonal [½, 1, 1, ..., ½].
  // See [Mac15], equation (1.20).
  FixedMatrix<double, N + 1, N + 1, /*use_heap=*/true> ᵝW =
      DiagonalMatrix<N + 1>(
          [](std::int64_t const i) { return (i == 0 || i == N) ? 0.5 : 1; });

  FixedMatrix<double, M + 1, N + 1, /*use_heap=*/true> Cₓ = ᵝT * ᵝW;

  // R is a diagonal (N + 1)×(N + 1) matrix.
  // See [Mac15], equation (1.25).
  FixedMatrix<double, N + 1, N + 1, /*use_heap=*/true> R =
      DiagonalMatrix<N + 1>([](std::int64_t const i) {
        if (i == 0) {
          return 1.0;
        } else if (i == N) {
          return 1.0 / N;
        } else {
          return 1.0 / (2 * i);
        }
      });

  // ᶠTᵀ is ᵝTᵀ with the last row removed.
  // See [Mac15], equation (1.22).
  FixedMatrix<double, N, M + 1, /*use_heap=*/true> ᶠTᵀ(uninitialized);
  for (std::int64_t i = 0; i < N; ++i) {
    for (std::int64_t j = 0; j <= M; ++j) {
      ᶠTᵀ(i, j) = ᵝT(j, i);
    }
  }

  // Cα is R * R * ᶠTᵀ * V (we do not assign it to a variable).

  CₓCα = Cₓ * R * S<N>() * ᶠTᵀ * V<M>();
}

template<ЧебышёвPicardMethod Method>
ЧебышёвPicardMatrices<Method, 2>::ЧебышёвPicardMatrices(
    FixedVector<double, Method::M + 1> const& nodes)
    : vCₓᵝCα(uninitialized), xCₓᵅCᵧᵝCα(uninitialized) {
  // We use the notation from [Mac15], section 1.4.4.

  // ᵅT is a (M + 1)×(N + 1) matrix of Чебышёв polynomials evaluated at nodes.
  // See [Mac15], equation (1.48a).
  FixedMatrix<double, M + 1, N + 1, /*use_heap=*/true> ᵅT(uninitialized);
  for (std::int64_t i = 0; i <= M; ++i) {
    auto const τᵢ = nodes[i];
    // The 0-degree polynomial is uniformly 1.
    ᵅT(i, 0) = 1;
    // The 0-degree polynomial is the identity.
    ᵅT(i, 1) = τᵢ;

    // We populate the rest of ᵅT using the recurrence relation.
    for (std::int64_t j = 2; j <= N; ++j) {
      ᵅT(i, j) = 2 * τᵢ * ᵅT(i, j - 1) - ᵅT(i, j - 2);
    }
  }

  // ᵝT is a (M + 1)×N matrix of Чебышёв polynomials evaluated at nodes.
  // See [Mac15], equation (1.48b).
  FixedMatrix<double, M + 1, N, /*use_heap=*/true> ᵝT(uninitialized);
  for (std::int64_t i = 0; i <= M; ++i) {
    for (std::int64_t j = 0; j < N; ++j) {
      ᵝT(i, j) = ᵅT(i, j);
    }
  }

  // ᵅW is a diagonal (N + 1)×(N + 1) matrix with diagonal [½, 1, 1, ..., ½].
  auto const ᵅW = DiagonalMatrix<N + 1>(
      [](std::int64_t const i) { return (i == 0 || i == N) ? 0.5 : 1; });

  // ᵝW is a diagonal N×N matrix with diagonal [½, 1, 1, ..., 1].
  auto const ᵝW =
      DiagonalMatrix<N>([](std::int64_t const i) { return i == 0 ? 0.5 : 1; });

  FixedMatrix<double, M + 1, N, /*use_heap=*/true> const vCₓ = ᵝT * ᵝW;
  FixedMatrix<double, M + 1, N + 1, /*use_heap=*/true> const xCₓ = ᵅT * ᵅW;

  // ᵅR and ᵝR are diagonal matrices of size (N + 1)×(N + 1) and N×N,
  // respectively.
  auto const ᵅR = DiagonalMatrix<N + 1>([](std::int64_t const i) {
    if (i == 0) {
      return 1.0;
    } else if (i == N) {
      return 1.0 / N;
    } else {
      return 1.0 / (2 * i);
    }
  });

  auto const ᵝR = DiagonalMatrix<N>(
      [](std::int64_t const i) { return i == 0 ? 1.0 : 1.0 / (2 * i); });

  // The S matrices are the same as in the first-order case.
  auto const ᵝS = S<N - 1>();
  auto const ᵅS = S<N>();

  // ᶠTᵀ is ᵝTᵀ with the last row removed.
  // See [Mac15], equation (1.22).
  FixedMatrix<double, N - 1, M + 1, /*use_heap=*/true> ᶠTᵀ(uninitialized);
  for (std::int64_t i = 0; i < N - 1; ++i) {
    for (std::int64_t j = 0; j <= M; ++j) {
      ᶠTᵀ(i, j) = ᵝT(j, i);
    }
  }

  FixedMatrix<double, N, M + 1, /*use_heap=*/true> const ᵝCα =
      ᵝR * ᵝS * ᶠTᵀ * V<M>();

  FixedMatrix<double, N + 1, N, /*use_heap=*/true> const ᵅCᵧ = ᵅR * ᵅS;

  vCₓᵝCα = vCₓ * ᵝCα;
  xCₓᵅCᵧᵝCα = xCₓ * ᵅCᵧ * ᵝCα;
}

template<std::int64_t M,
         typename IndependentVariable,
         typename... DependentVariable>
FixedVector<DirectSum<DependentVariable...>, M + 1, true>
ЧебышёвPicardIterationState<
    M,
    ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                   DependentVariable...>>::
    UninitializedDependentVariableMatrix(
        InitialValueProblem<ODE> const& problem) {
  return FixedVector<DirectSum<DependentVariable...>, M + 1, true>(
      uninitialized);
}

template<std::int64_t M,
         typename IndependentVariable,
         typename... DependentVariable>
FixedVector<DirectSum<Derivative<DependentVariable, IndependentVariable>...>,
            M + 1,
            true>
ЧебышёвPicardIterationState<
    M,
    ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                   DependentVariable...>>::
    UninitializedRightHandSideMatrix(InitialValueProblem<ODE> const& problem) {
  return FixedVector<
      DirectSum<Derivative<DependentVariable, IndependentVariable>...>,
      M + 1,
      true>(uninitialized);
}

template<std::int64_t M,
         typename IndependentVariable,
         typename... DependentVariable>
IndependentVariable ЧебышёвPicardIterationState<
    M,
    ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                   DependentVariable...>>::
    GetIndependentVariable(ODE::State const& state) {
  return state.s.value;
}

template<std::int64_t M, typename DependentVariable>
std::pair<UnboundedMatrix<DependentVariable>,
          UnboundedMatrix<Derivative<DependentVariable, Instant>>>
ЧебышёвPicardIterationState<
    M,
    ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>>::
    UninitializedDependentVariableMatrix(
        InitialValueProblem<ODE> const& problem) {
  return {UnboundedMatrix<DependentVariable>(
              M + 1, problem.initial_state.positions.size(), uninitialized),
          UnboundedMatrix<Derivative<DependentVariable, Instant>>(
              M + 1, problem.initial_state.velocities.size(), uninitialized)};
}

template<std::int64_t M, typename DependentVariable>
UnboundedMatrix<Derivative<DependentVariable, Instant, 2>>
ЧебышёвPicardIterationState<
    M,
    ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>>::
    UninitializedRightHandSideMatrix(InitialValueProblem<ODE> const& problem) {
  return UnboundedMatrix<Derivative<DependentVariable, Instant, 2>>(
      M + 1, problem.initial_state.positions.size(), uninitialized);
}

template<std::int64_t M, typename DependentVariable>
Instant ЧебышёвPicardIterationState<
    M,
    ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>>::
    GetIndependentVariable(ODE::State const& state) {
  return state.time.value;
}

template<ЧебышёвPicardMethod Method, typename ODE_>
absl::Status ЧебышёвPicardIntegrator<Method, ODE_>::Instance::Solve(
    ODE::IndependentVariable const& t_final) {
  using DependentVariables = typename ODE::DependentVariables;
  using DependentVariableDerivatives =
      typename ODE::DependentVariableDerivatives;
  using State = typename ODE::State;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto const& equation = this->equation_;
  auto const& step = this->step_;

  // Argument checks.
  Sign const integration_direction = Sign(step);
  if (integration_direction.is_positive()) {
    // Integrating forward.
    CHECK_LT((ЧебышёвPicardIterationState<M, ODE>::GetIndependentVariable(
                 current_state)),
             t_final);
  } else {
    // Integrating backward.
    CHECK_GT((ЧебышёвPicardIterationState<M, ODE>::GetIndependentVariable(
                 current_state)),
             t_final);
  }

  while (integration_direction.is_positive()
             ? ЧебышёвPicardIterationState<M, ODE>::GetIndependentVariable(
                   current_state) < t_final
             : ЧебышёвPicardIterationState<M, ODE>::GetIndependentVariable(
                   current_state) > t_final) {
    auto const t_initial =
        ЧебышёвPicardIterationState<M, ODE>::GetIndependentVariable(
            current_state);

    // Rescale the nodes for feeding into the compute_derivative function.
    t_.clear();
    for (double const node : integrator_.nodes_) {
      t_.push_back(t_initial + (0.5 * node + 0.5) * step);
    }

    // Set the boundary condition and store it in boundary_.
    if constexpr (ODE::order == 1) {
      boundary_[0] = StripDoublePrecision(current_state.y);
      for (std::int64_t i = 1; i <= M; ++i) {
        boundary_[i] = boundary_[0];
      }
    } else {
      for (std::int64_t i = 0; i <= M; ++i) {
        for (std::int64_t j = 0; j < boundary_.first.columns(); ++j) {
          boundary_.second(i, j) = current_state.velocities[j].value;
          boundary_.first(i, j) =
              current_state.positions[j].value +
              current_state.velocities[j].value * (t_[i] - t_[0]);
        }
      }
    }

    // A good starting guess for X⁰ is uniform position (for first-order) or
    // velocity (for second-order); as it happens that's what we just set
    // boundary_ to.
    Xⁱ_ = boundary_;

    bool previous_should_stop = false;
    bool converged = false;
    for (int64_t iteration = 0; iteration < params_.max_iterations;
         ++iteration) {
      // Evaluate the right hand side of the equation.
      for (std::int64_t i = 0; i <= M; ++i) {
        if constexpr (ODE::order == 1) {
          auto const& y = Xⁱ_[i];
          DependentVariableDerivatives yʹᵢ;
          RETURN_IF_ERROR(equation.compute_derivative(t_[i], y, yʹᵢ));

          // Store it in yʹ.
          f_[i] = std::move(yʹᵢ);
        } else {
          DependentVariables y;
          DependentVariableDerivatives yʹ;
          for (std::int64_t j = 0; j < Xⁱ_.first.columns(); ++j) {
            y.push_back(Xⁱ_.first(i, j));
            yʹ.push_back(Xⁱ_.second(i, j));
          }
          typename ODE::DependentVariableDerivatives2 yʺᵢ(Xⁱ_.first.columns());
          RETURN_IF_ERROR(equation.compute_acceleration(t_[i], y, yʹ, yʺᵢ));
          for (std::int64_t j = 0; j < Xⁱ_.first.columns(); ++j) {
            f_(i, j) = yʺᵢ[j];
          }
        }
      }

      // Compute new x.
      if constexpr (ODE::order == 1) {
        Xⁱ⁺¹_ = integrator_.matrices_.CₓCα * (0.5 * step * f_) + boundary_;
      } else {
        Xⁱ⁺¹_.second = integrator_.matrices_.vCₓᵝCα * (0.5 * step * f_) + boundary_.second;
        Xⁱ⁺¹_.first = integrator_.matrices_.xCₓᵅCᵧᵝCα * (0.25 * step * step * f_) + boundary_.first;
      }

      // Check for convergence by applying the stopping criterion.
      bool should_stop = true;
      for (std::int64_t i = 0; i <= M; ++i) {
        if constexpr (ODE::order == 1) {
          should_stop =
              should_stop && params_.stopping_criterion(Xⁱ⁺¹_[i] - Xⁱ_[i]);
        } else {
          typename ODE::State::Error Δ;
          for (std::int64_t j = 0; j < Xⁱ_.first.columns(); ++j) {
            Δ.position_error.push_back((Xⁱ⁺¹_.first)(i, j) - (Xⁱ_.first)(i, j));
            Δ.velocity_error.push_back((Xⁱ⁺¹_.second)(i, j) -
                                       (Xⁱ_.second)(i, j));
          }
          should_stop = should_stop && params_.stopping_criterion(Δ);
        }
      }
      Xⁱ_ = std::move(Xⁱ⁺¹_);

      // We require that ‖Xⁱ⁺¹ - Xⁱ‖ and ‖Xⁱ - Xⁱ⁻¹‖ are _both_ less than
      // the given tolerance to account for nonlinearity issues (as suggested in
      // [BJ12]).
      if (should_stop && previous_should_stop) {
        converged = true;
        break;
      }

      previous_should_stop = should_stop;
      RETURN_IF_STOPPED;
    }

    if (converged) {
      // We have successfully converged!
      for (std::int64_t i = 0; i <= M; ++i) {
        if constexpr (ODE::order == 1) {
          append_state(State(t_[i], Xⁱ_[i]));
        } else {
          DependentVariables q;
          DependentVariableDerivatives v;
          for (std::int64_t j = 0; j < Xⁱ_.first.columns(); ++j) {
            q.push_back((Xⁱ_.first)(i, j));
            v.push_back((Xⁱ_.second)(i, j));
          }
          append_state(State(t_[i], q, v));
        }
      }

      // Set the current state to the final state we appended.
      if constexpr (ODE::order == 1) {
        current_state = State(t_[M], Xⁱ_[M]);
      } else {
        DependentVariables q;
        DependentVariableDerivatives v;
        for (std::int64_t j = 0; j < Xⁱ_.first.columns(); ++j) {
          q.push_back((Xⁱ_.first)(M, j));
          v.push_back((Xⁱ_.second)(M, j));
        }
        current_state = State(t_[M], q, v);
      }
      RETURN_IF_STOPPED;
    } else {
      // We failed to converge.
      return absl::Status(absl::StatusCode::kFailedPrecondition,
                          "Чебышёв-Picard iteration failed to converge.");
    }
  }

  return absl::OkStatus();
}

template<ЧебышёвPicardMethod Method, typename ODE_>
ЧебышёвPicardIntegrator<Method, ODE_> const&
ЧебышёвPicardIntegrator<Method, ODE_>::Instance::integrator() const {
  return integrator_;
}

template<ЧебышёвPicardMethod Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIntegrator<Method, ODE_>::Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<ЧебышёвPicardMethod Method, typename ODE_>
ЧебышёвPicardIntegrator<Method, ODE_>::Instance::Instance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    ЧебышёвPicardIntegrator const& integrator,
    ЧебышёвPicardIterationParams<ODE> const& params)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      integrator_(integrator),
      boundary_(ЧебышёвPicardIterationState<Method::M, ODE_>::
                    UninitializedDependentVariableMatrix(problem)),
      Xⁱ_(ЧебышёвPicardIterationState<Method::M, ODE_>::
              UninitializedDependentVariableMatrix(problem)),
      Xⁱ⁺¹_(ЧебышёвPicardIterationState<Method::M, ODE_>::
                UninitializedDependentVariableMatrix(problem)),
      f_(ЧебышёвPicardIterationState<Method::M, ODE_>::
             UninitializedRightHandSideMatrix(problem)),
      params_(params) {
  t_.reserve(M + 1);
}

template<ЧебышёвPicardMethod Method, typename ODE_>
ЧебышёвPicardIntegrator<Method, ODE_>::ЧебышёвPicardIntegrator()
    : nodes_(Nodes<M>()), matrices_(nodes_) {}

template<ЧебышёвPicardMethod Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIntegrator<Method, ODE_>::NewInstance(
    InitialValueProblem<ODE_> const& problem,
    AppendState const& append_state,
    Time const& step) const {
  return NewInstance(
      problem, append_state, step, ЧебышёвPicardIterationParams<ODE>());
}

template<ЧебышёвPicardMethod Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIntegrator<Method, ODE_>::NewInstance(
    InitialValueProblem<ODE_> const& problem,
    AppendState const& append_state,
    Time const& step,
    ЧебышёвPicardIterationParams<ODE> const& params) const {
  // Cannot use `make_not_null_unique` because the constructor of `Instance` is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this, params));
}

template<ЧебышёвPicardMethod Method, typename ODE_>
void ЧебышёвPicardIntegrator<Method, ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  LOG(FATAL) << "Serialization of ЧебышёвPicardIntegrator is not yet supported";
  std::abort();
}

}  // namespace internal
}  // namespace _чебышёв_picard_integrator
}  // namespace integrators
}  // namespace principia
