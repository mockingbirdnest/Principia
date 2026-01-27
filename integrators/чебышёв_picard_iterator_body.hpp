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

template<typename ODE_>
ODE_::DependentVariables DependentVariablesFromMatrixRow(
    UnboundedMatrix<double> const& matrix,
    std::int64_t const row) {
  std::int64_t j = 0;
  typename ODE_::DependentVariables y;
  for_all_of(y).loop([&matrix, &j, row](auto& yⱼ) {
    yⱼ = matrix(row, j++) * Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
  return y;
}

template<typename ODE_>
void DependentVariablesToMatrixRow(typename ODE_::DependentVariables const& y,
                                   std::int64_t const row,
                                   UnboundedMatrix<double>& matrix) {
  std::int64_t j = 0;
  for_all_of(y).loop([row, &matrix, &j](auto const& yⱼ) {
    matrix(row, j++) = yⱼ / Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
}

template<typename ODE_>
void DependentVariableDerivativesToMatrixRow(
    typename ODE_::DependentVariableDerivatives const& y,
    std::int64_t const row,
    UnboundedMatrix<double>& matrix) {
  std::int64_t j = 0;
  for_all_of(y).loop([row, &matrix, &j](auto const& yⱼ) {
    matrix(row, j++) = yⱼ / Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
}

// Returns max|aᵢⱼ|.
inline double MaxNorm(UnboundedMatrix<double> const& A) {
  double norm = 0.0;
  for (std::int64_t i = 0; i < A.rows(); i++) {
    for (std::int64_t j = 0; j < A.columns(); j++) {
      norm = std::max(norm, std::abs(A(i, j)));
    }
  }
  return norm;
}

template<typename ODE_>
absl::Status ЧебышёвPicardIterator<ODE_>::Instance::Solve(
    ODE::IndependentVariable const& t_final) {
  using DependentVariables = typename ODE::DependentVariables;
  using DependentVariableDerivatives =
      typename ODE::DependentVariableDerivatives;
  using State = typename ODE::State;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto const& equation = this->equation_;
  auto const& step = this->step_;
  auto const& params = integrator_.params();
  auto const n = std::tuple_size<DependentVariables>::value;

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
      CₓX₀_(0, j++) = yⱼ.value / Unit<decltype(yⱼ.value)>;
    });
    for (std::int64_t i = 1; i <= params.M; i++) {
      for (std::int64_t j = 0; j < n; j++) {
        CₓX₀_(i, j) = CₓX₀_(0, j);
      }
    }

    // A good starting guess for X⁰ is uniform current_state.y; as it happens
    // that's what we just set CₓX₀_ to.
    Xⁱ_ = CₓX₀_;

    double prev_norm = std::numeric_limits<float>::infinity();
    bool converged = false;
    for (int64_t iteration = 0; iteration < params.max_iterations;
         iteration++) {
      // Evaluate the right hand side of the equation.
      for (int64_t i = 0; i < Xⁱ_.rows(); i++) {
        auto const y = DependentVariablesFromMatrixRow<ODE>(Xⁱ_, i);
        DependentVariableDerivatives yʹᵢ;
        RETURN_IF_ERROR(equation.compute_derivative(t_[i], y, yʹᵢ));

        // Store it in yʹ.
        DependentVariableDerivativesToMatrixRow<ODE>(yʹᵢ, i, yʹ_);
      }

      // Compute new x.
      Xⁱ⁺¹_ = integrator_.CₓCα_ * ((step / Second) * yʹ_) + CₓX₀_;

      // Check for convergence by computing the ∞-norm.
      const double norm = MaxNorm(Xⁱ⁺¹_ - Xⁱ_);
      Xⁱ_ = std::move(Xⁱ⁺¹_);

      if (std::max(norm, prev_norm) < params.stopping_criterion) {
        converged = true;
        break;
      }

      prev_norm = norm;
    }

    if (converged) {
      // We have successfully converged!
      for (std::int64_t i = 0; i < Xⁱ_.rows(); i++) {
        append_state(
            State(t_[i], DependentVariablesFromMatrixRow<ODE>(Xⁱ_, i)));
      }

      // Set the current state to the final state we appended.
      current_state =
          State(t_[Xⁱ_.rows() - 1],
                DependentVariablesFromMatrixRow<ODE>(Xⁱ_, Xⁱ_.rows() - 1));
    } else {
      // We failed to converge.
      return absl::Status(absl::StatusCode::kFailedPrecondition,
                          "Чебышёв-Picard iteration failed to converge.");
    }
  }

  return absl::OkStatus();
}

template<typename ODE_>
ЧебышёвPicardIterator<ODE_> const&
ЧебышёвPicardIterator<ODE_>::Instance::integrator() const {
  return integrator_;
}

template<typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIterator<ODE_>::Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename ODE_>
ЧебышёвPicardIterator<ODE_>::Instance::Instance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    ЧебышёвPicardIterator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      integrator_(integrator),
      t_(),
      CₓX₀_(integrator.params_.M + 1,
            std::tuple_size<typename ODE::DependentVariables>::value,
            uninitialized),
      Xⁱ_(CₓX₀_.rows(), CₓX₀_.columns(), uninitialized),
      Xⁱ⁺¹_(Xⁱ_.rows(), Xⁱ_.columns(), uninitialized),
      yʹ_(Xⁱ_.rows(), Xⁱ_.columns(), uninitialized) {
  t_.reserve(integrator.nodes_.size());
}

template<typename ODE_>
ЧебышёвPicardIterator<ODE_>::ЧебышёвPicardIterator(
    ЧебышёвPicardIterationParams const& params)
    : params_(params),
      nodes_(params.M + 1, uninitialized),
      CₓCα_(params.M + 1, params.N + 1, uninitialized) {
  // We use the notation from [Mac15], section 1.4.3.
  const std::int64_t M = params_.M;
  const std::int64_t N = params_.N;
  CHECK_GE(M, 1);
  CHECK_GE(N, 1);

  // Populate nodes.
  for (std::int64_t i = 0; i <= M; i++) {
    nodes_[i] = -Cos(π / M * i * Radian);
  }

  // ᵝT is a (M + 1)×(N + 1) matrix of Чебышёв polynomials evaluated at nodes.
  // See [Mac15], equation (1.20).
  UnboundedMatrix<double> ᵝT(M + 1, N + 1, uninitialized);

  for (std::int64_t i = 0; i <= M; i++) {
    const auto τᵢ = nodes_[i];
    // The 0-degree polynomial is uniformly 1.
    ᵝT(i, 0) = 1;
    // The 0-degree polynomial is the identity.
    ᵝT(i, 1) = τᵢ;

    // We populate the rest of ᵝT using the recurrence relation.
    for (std::int64_t j = 2; j <= N; j++) {
      ᵝT(i, j) = 2 * τᵢ * ᵝT(i, j - 1) - ᵝT(i, j - 2);

      // Make sure the zeroes are actually zero.
      if (std::abs(ᵝT(i, j)) < 1e-14)
        ᵝT(i, j) = 0;
    }
  }

  // ᵝW is a diagonal (N + 1)×(N + 1) matrix with diagonal [½, 1, 1, ..., ½].
  // See [Mac15], equation (1.20).
  UnboundedMatrix<double> ᵝW(N + 1, N + 1);
  ᵝW(0, 0) = 0.5;
  ᵝW(N, N) = 0.5;
  for (std::int64_t i = 1; i < N; i++) {
    ᵝW(i, i) = 1;
  }

  UnboundedMatrix<double> Cₓ = ᵝT * ᵝW;

  // R is a diagonal (N + 1)×(N + 1) matrix.
  // See [Mac15], equation (1.25).
  UnboundedMatrix<double> R(N + 1, N + 1);
  R(0, 0) = 1;
  R(N, N) = 1.0 / N;
  for (std::int64_t i = 1; i < N; i++) {
    R(i, i) = 1.0 / (2 * i);
  }

  // S is an (N + 1)×N matrix.
  // See equation 1.26 in [Mac15].
  UnboundedMatrix<double> S(N + 1, N);
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

  // ᶠT is ᵝTᵀ with the last row removed.
  // See [Mac15], equation (1.22).
  UnboundedMatrix<double> ᶠT(N, M + 1, uninitialized);
  for (std::int64_t i = 0; i < N; i++) {
    for (std::int64_t j = 0; j <= M; j++) {
      ᶠT(i, j) = ᵝT(j, i);
    }
  }

  // ᵀV is 1/M * ᵝW (we do not assign it to a variable).
  // Cα is r * s * ᶠT * ᵀV (we do not assign it to a variable).

  CₓCα_ = 1.0 / M * Cₓ * R * S * ᶠT * ᵝW;
}

template<typename ODE_>
ЧебышёвPicardIterationParams const& ЧебышёвPicardIterator<ODE_>::params()
    const {
  return params_;
}

template<typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ЧебышёвPicardIterator<ODE_>::NewInstance(
    InitialValueProblem<ODE_> const& problem,
    AppendState const& append_state,
    Time const& step) const {
  // Cannot use `make_not_null_unique` because the constructor of `Instance` is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename ODE_>
void ЧебышёвPicardIterator<ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(serialization::FixedStepSizeIntegrator::CHEBYSHEV_PICARD);
  auto& message_params = *message->mutable_chebyshev_picard_params();
  message_params.set_m(params_.M);
  message_params.set_n(params_.N);
  message_params.set_max_iterations(params_.max_iterations);
  message_params.set_stopping_criterion(params_.stopping_criterion);
}

}  // namespace internal
}  // namespace _чебышёв_picard_iterator
}  // namespace integrators
}  // namespace principia
