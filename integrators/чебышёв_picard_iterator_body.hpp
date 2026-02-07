#pragma once

#include "integrators/чебышёв_picard_iterator.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <utility>

#include "base/for_all_of.hpp"
#include "base/status_utilities.hpp"  // 🧙 For RETURN_IF_ERROR.
#include "base/tags.hpp"
#include "geometry/sign.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/matrix_views.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace integrators {
namespace _чебышёв_picard_iterator {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::base::_tags;
using namespace principia::geometry::_sign;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_matrix_views;
using namespace principia::quantities::_si;

template<typename ODE, std::int64_t M>
ODE::DependentVariables DependentVariablesFromMatrixRow(
    FixedMatrix<double,
                M,
                std::tuple_size_v<typename ODE::DependentVariables>> const&
        matrix,
    std::int64_t const row) {
  std::int64_t j = 0;
  typename ODE::DependentVariables y;
  for_all_of(y).loop([&j, &matrix, row](auto& yⱼ) {
    yⱼ = matrix(row, j++) * si::Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
  return y;
}

template<typename ODE, std::int64_t M>
void DependentVariablesToMatrixRow(
    typename ODE::DependentVariables const& y,
    std::int64_t const row,
    FixedMatrix<double, M, std::tuple_size_v<typename ODE::DependentVariables>>&
        matrix) {
  std::int64_t j = 0;
  for_all_of(y).loop([row, &matrix, &j](auto const& yⱼ) {
    matrix(row, j++) = yⱼ / si::Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
}

template<typename ODE, std::int64_t M>
void DependentVariableDerivativesToMatrixRow(
    typename ODE::DependentVariableDerivatives const& y,
    std::int64_t const row,
    FixedMatrix<double, M, std::tuple_size_v<typename ODE::DependentVariables>>&
        matrix) {
  std::int64_t j = 0;
  for_all_of(y).loop([row, &matrix, &j](auto const& yⱼ) {
    matrix(row, j++) = yⱼ / si::Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
}

// Returns max|aᵢⱼ|.
template<std::int64_t M, std::int64_t N>
double LInfinityNorm(FixedMatrix<double, M, N> const& A) {
  double norm = 0.0;
  for (std::int64_t i = 0; i < M; ++i) {
    for (std::int64_t j = 0; j < N; ++j) {
      norm = std::max(norm, std::abs(A(i, j)));
    }
  }
  return norm;
}

template<typename Method, typename ODE_>
absl::Status ЧебышёвPicardIterator<Method, ODE_>::Instance::Solve(
    ODE::IndependentVariable const& t_final) {
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
    for (const double node : integrator_.nodes_) {
      t_.push_back(t_initial + (0.5 * node + 0.5) * step);
    }

    // Set the boundary condition and store it in CₓX₀_.
    std::int64_t j = 0;
    for_all_of(current_state.y).loop([this, &j](auto const& yⱼ) {
      CₓX₀_(0, j++) = yⱼ.value / si::Unit<decltype(yⱼ.value)>;
    });
    for (std::int64_t i = 1; i <= M; ++i) {
      for (std::int64_t j = 0; j < n; ++j) {
        CₓX₀_(i, j) = CₓX₀_(0, j);
      }
    }

    // A good starting guess for X⁰ is uniform current_state.y; as it happens
    // that's what we just set CₓX₀_ to.
    Xⁱ_ = CₓX₀_;

    double previous_norm = std::numeric_limits<float>::infinity();
    bool converged = false;
    for (int64_t iteration = 0; iteration < params_.max_iterations;
         ++iteration) {
      // Evaluate the right hand side of the equation.
      for (int64_t i = 0; i < Xⁱ_.rows(); ++i) {
        auto const y = DependentVariablesFromMatrixRow<ODE, M + 1>(Xⁱ_, i);
        DependentVariableDerivatives yʹᵢ;
        RETURN_IF_ERROR(equation.compute_derivative(t_[i], y, yʹᵢ));

        // Store it in yʹ.
        DependentVariableDerivativesToMatrixRow<ODE, M + 1>(yʹᵢ, i, yʹ_);
      }

      // Compute new x.
      Xⁱ⁺¹_ = integrator_.CₓCα_ * (0.5 * step / Second * yʹ_) + CₓX₀_;

      // Check for convergence by computing the ∞-norm.
      const double norm = LInfinityNorm(Xⁱ⁺¹_ - Xⁱ_);
      Xⁱ_ = std::move(Xⁱ⁺¹_);

      // We require that ‖Xⁱ⁺¹ - Xⁱ‖ and ‖Xⁱ - Xⁱ⁻¹‖ are _both_ less than
      // the given tolerance to account for nonlinearity issues (as suggested in
      // [BJ12]).
      if (std::max(norm, previous_norm) < params_.stopping_criterion) {
        converged = true;
        break;
      }

      previous_norm = norm;
      RETURN_IF_STOPPED;
    }

    if (converged) {
      // We have successfully converged!
      for (std::int64_t i = 0; i < Xⁱ_.rows(); ++i) {
        append_state(
            State(t_[i], DependentVariablesFromMatrixRow<ODE>(Xⁱ_, i)));
      }

      // Set the current state to the final state we appended.
      current_state =
          State(t_[Xⁱ_.rows() - 1],
                DependentVariablesFromMatrixRow<ODE>(Xⁱ_, Xⁱ_.rows() - 1));
      RETURN_IF_STOPPED;
    } else {
      // We failed to converge.
      return absl::Status(absl::StatusCode::kFailedPrecondition,
                          "Чебышёв-Picard iteration failed to converge.");
    }
  }

  return absl::OkStatus();
}

template<typename Method, typename ODE_>
ЧебышёвPicardIterator<Method, ODE_> const&
ЧебышёвPicardIterator<Method, ODE_>::Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIterator<Method, ODE_>::Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
ЧебышёвPicardIterator<Method, ODE_>::Instance::Instance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    ЧебышёвPicardIterator const& integrator,
    ЧебышёвPicardIterationParams const& params)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      params_(params),
      integrator_(integrator),
      CₓX₀_(uninitialized),
      Xⁱ_(uninitialized),
      Xⁱ⁺¹_(uninitialized),
      yʹ_(uninitialized) {
  t_.reserve(M + 1);
}

template<typename Method, typename ODE_>
ЧебышёвPicardIterator<Method, ODE_>::ЧебышёвPicardIterator()
    : nodes_(uninitialized), CₓCα_(uninitialized) {
  // We use the notation from [Mac15], section 1.4.3.

  // Populate nodes.
  for (std::int64_t i = 0; i <= M; i++) {
    nodes_[i] = -Cos(π / M * i * Radian);
  }

  // ᵝT is a (M + 1)×(N + 1) matrix of Чебышёв polynomials evaluated at nodes.
  // See [Mac15], equation (1.20).
  FixedMatrix<double, M + 1, N + 1> ᵝT(uninitialized);

  for (std::int64_t i = 0; i <= M; i++) {
    const auto τᵢ = nodes_[i];
    // The 0-degree polynomial is uniformly 1.
    ᵝT(i, 0) = 1;
    // The 0-degree polynomial is the identity.
    ᵝT(i, 1) = τᵢ;

    // We populate the rest of ᵝT using the recurrence relation.
    for (std::int64_t j = 2; j <= N; j++) {
      ᵝT(i, j) = 2 * τᵢ * ᵝT(i, j - 1) - ᵝT(i, j - 2);
    }
  }

  // ᵝW is a diagonal (N + 1)×(N + 1) matrix with diagonal [½, 1, 1, ..., ½].
  // See [Mac15], equation (1.20).
  FixedMatrix<double, N + 1, N + 1> ᵝW;
  ᵝW(0, 0) = 0.5;
  ᵝW(N, N) = 0.5;
  for (std::int64_t i = 1; i < N; ++i) {
    ᵝW(i, i) = 1;
  }

  FixedMatrix<double, M + 1, N + 1> Cₓ = ᵝT * ᵝW;

  // R is a diagonal (N + 1)×(N + 1) matrix.
  // See [Mac15], equation (1.25).
  FixedMatrix<double, N + 1, N + 1> R;
  R(0, 0) = 1;
  R(N, N) = 1.0 / N;
  for (std::int64_t i = 1; i < N; i++) {
    R(i, i) = 1.0 / (2 * i);
  }

  // S is an (N + 1)×N matrix.
  // See equation 1.26 in [Mac15].
  FixedMatrix<double, N + 1, N> S;
  S(0, 0) = 1;
  S(0, 1) = -0.5;
  for (std::int64_t k = 2; k < N; k++) {
    S(0, k) = (k % 2 == 1 ? 1 : -1) * (1.0 / (k - 1) - 1.0 / (k + 1));
  }
  for (std::int64_t i = 0; i < N; i++) {
    S(i + 1, i) = 1;
  }
  for (std::int64_t i = 1; i + 2 < N; i++) {
    S(i, i + 1) = -1;
  }

  // ᶠTᵀ is ᵝTᵀ with the last row removed.
  // See [Mac15], equation (1.22).
  FixedMatrix<double, N, M + 1> ᶠTᵀ(uninitialized);
  for (std::int64_t i = 0; i < N; i++) {
    for (std::int64_t j = 0; j <= M; j++) {
      ᶠTᵀ(i, j) = ᵝT(j, i);
    }
  }

  // V is is a diagonal (M + 1)×(M + 1) matrix with diagonal [1/M, 2/M, 2/M,
  // ..., 1/M].
  FixedMatrix<double, M + 1, M + 1> V;
  constexpr double one_over_M = 1.0 / M;
  V(0, 0) = one_over_M;
  V(M, M) = one_over_M;
  for (std::int64_t i = 1; i < M; i++) {
    V(i, i) = 2.0 * one_over_M;
  }

  // Cα is R * R * ᶠTᵀ * V (we do not assign it to a variable).

  CₓCα_ = Cₓ * R * S * ᶠTᵀ * V;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIterator<Method, ODE_>::NewInstance(
    InitialValueProblem<ODE_> const& problem,
    AppendState const& append_state,
    Time const& step) const {
  return NewInstance(
      problem, append_state, step, ЧебышёвPicardIterationParams());
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIterator<Method, ODE_>::NewInstance(
    InitialValueProblem<ODE_> const& problem,
    AppendState const& append_state,
    Time const& step,
    ЧебышёвPicardIterationParams const& params) const {
  // Cannot use `make_not_null_unique` because the constructor of `Instance` is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this, params));
}

template<typename Method, typename ODE_>
void ЧебышёвPicardIterator<Method, ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  LOG(FATAL) << "Serialization of ЧебышёвPicardIntegrator is not yet supported";
  std::abort();
}

}  // namespace internal
}  // namespace _чебышёв_picard_iterator
}  // namespace integrators
}  // namespace principia
