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
  return DiagonalMatrix<M + 1>([one_over_M](std::int64_t const i) {
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
    CHECK_LT(current_state.s.value, t_final);
  } else {
    // Integrating backward.
    CHECK_GT(current_state.s.value, t_final);
  }

  while (integration_direction.is_positive()
             ? current_state.s.value < t_final
             : current_state.s.value > t_final) {
    auto const t_initial = current_state.s.value;

    // Rescale the nodes for feeding into the compute_derivative function.
    t_.clear();
    for (double const node : integrator_.nodes_) {
      t_.push_back(t_initial + (0.5 * node + 0.5) * step);
    }

    // Set the boundary condition and store it in CₓX₀_.
    CₓX₀_[0] = StripDoublePrecision(current_state.y);
    for (std::int64_t i = 1; i <= M; ++i) {
      CₓX₀_[i] = CₓX₀_[0];
    }

    // A good starting guess for X⁰ is uniform current_state.y; as it happens
    // that's what we just set CₓX₀_ to.
    Xⁱ_ = CₓX₀_;

    bool previous_should_stop = false;
    bool converged = false;
    for (int64_t iteration = 0; iteration < params_.max_iterations;
         ++iteration) {
      // Evaluate the right hand side of the equation.
      for (std::int64_t i = 0; i <= M; ++i) {
        auto const& y = Xⁱ_[i];
        DependentVariableDerivatives yʹᵢ;
        RETURN_IF_ERROR(equation.compute_derivative(t_[i], y, yʹᵢ));

        // Store it in yʹ.
        yʹ_[i] = std::move(yʹᵢ);
      }

      // Compute new x.
      Xⁱ⁺¹_ = integrator_.matrices_.CₓCα * (0.5 * step * yʹ_) + CₓX₀_;

      // Check for convergence by applying the stopping criterion.
      bool should_stop = true;
      for (std::int64_t i = 0; i <= M; ++i) {
        should_stop = should_stop &&
                      params_.stopping_criterion((Xⁱ⁺¹_[i] - Xⁱ_[i]));
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
        append_state(State(t_[i], Xⁱ_[i]));
      }

      // Set the current state to the final state we appended.
      current_state = State(t_[M], Xⁱ_[M]);
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
      params_(params),
      CₓX₀_(uninitialized),
      Xⁱ_(uninitialized),
      Xⁱ⁺¹_(uninitialized),
      yʹ_(uninitialized) {
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
