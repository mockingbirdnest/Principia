#pragma once

#include "base/for_all_of.hpp"
#include "base/status_utilities.hpp"  // 🧙 For RETURN_IF_ERROR.
#include "base/tags.hpp"
#include "integrators/chebyshev_picard_iterator.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/matrix_computations.hpp"  // For eigenvalues.
#include "numerics/matrix_views.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace integrators {

using namespace principia::base::_for_all_of;
using namespace principia::base::_tags;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_matrix_views;
using namespace principia::quantities::_si;

namespace _chebyshev_picard_iterator {
namespace internal {

template <typename ODE_>
ODE_::DependentVariables DependentVariablesFromMatrixRow(
    UnboundedMatrix<double> const& matrix, int row) {
  int j = 0;
  typename ODE_::DependentVariables y;
  for_all_of(y).loop([&matrix, &j, row](auto& yⱼ) {
    yⱼ = matrix(row, j++) * Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
  return y;
}

template <typename ODE_>
void DependentVariablesToMatrixRow(typename ODE_::DependentVariables const& y,
                                   int row, UnboundedMatrix<double>& matrix) {
  int j = 0;
  for_all_of(y).loop([row, &matrix, &j](auto const& yⱼ) {
    matrix(row, j++) = yⱼ / Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
}

template <typename ODE_>
void DependentVariableDerivativesToMatrixRow(
    typename ODE_::DependentVariableDerivatives const& y, int row,
    UnboundedMatrix<double>& matrix) {
  int j = 0;
  for_all_of(y).loop([row, &matrix, &j](auto const& yⱼ) {
    matrix(row, j++) = yⱼ / Unit<std::remove_reference_t<decltype(yⱼ)>>;
  });
}

template <typename ODE_>
absl::Status ChebyshevPicardIterator<ODE_>::Instance::Solve(
    ODE::IndependentVariable const& t_final) {
  using IndependentVariable = typename ODE::IndependentVariable;
  using DependentVariables = typename ODE::DependentVariables;
  using DependentVariableDerivatives =
      typename ODE::DependentVariableDerivatives;
  using State = typename ODE::State;

  auto& append_state = this->append_state_;
  auto const& equation = this->equation_;
  auto const& step = this->step_;
  auto const& params = integrator_.params();
  auto& current_state = this->current_state_;
  while (current_state.s.value < t_final) {
    auto const t_initial = current_state.s.value;
    CHECK_LT(t_initial, t_final);

    // Rescale the nodes for feeding into the compute_derivative function.
    std::vector<IndependentVariable> t;
    t.reserve(integrator_.nodes_.size());
    for (const double node : integrator_.nodes_) {
      t.push_back(t_initial + (0.5 * node + 0.5) * step);
    }

    // x is an (N + 1)×n matrix, where n is the dimension of the ODE's
    // dependent variable.
    UnboundedMatrix<double> x(integrator_.cx_.columns(),
                              std::tuple_size<DependentVariables>::value);

    // Set the boundary condition and store it in cₓx₀.
    int j = 0;
    for_all_of(current_state.y).loop([&x, &j](auto const& yⱼ) {
      x(0, j++) = yⱼ.value / Unit<decltype(yⱼ.value)>;
    });

    const UnboundedMatrix<double> cₓx₀ = integrator_.cx_ * 2 * x;

    // Set the initial value of x (this is x⁰, with superscript 0) to the
    // current state.
    for (int i = 1; i < x.rows(); i++) {
      for (int j = 0; j < x.columns(); j++) {
        x(i, j) = x(0, j);
      }
    }

    // The computed derivative. We will multiply this by 0.5 * step to get g.
    UnboundedMatrix<double> yʹ(x.rows(), x.columns(), uninitialized);

    double prev_norm = std::numeric_limits<float>::infinity();
    bool converged = false;
    for (int iteration = 0; iteration < params.max_iterations; iteration++) {
      // Evaluate the right hand side of the equation.
      for (int i = 0; i < x.rows(); i++) {
        auto const y = DependentVariablesFromMatrixRow<ODE>(x, i);
        DependentVariableDerivatives yʹᵢ;
        RETURN_IF_ERROR(equation.compute_derivative(t[i], y, yʹᵢ));

        // Store it in yʹ.
        DependentVariableDerivativesToMatrixRow<ODE>(yʹᵢ, i, yʹ);
      }

      // Compute new x.
      const UnboundedMatrix<double> new_x =
          integrator_.cx_cα_ * (step / Second) * yʹ + cₓx₀;

      // Check for convergence by computing the ∞-norm.
      double norm = 0.0;
      for (int i = 0; i < x.rows(); i++) {
        for (int j = 0; j < x.columns(); j++) {
          norm = std::max(
              norm, std::abs(new_x(i, j) - x(i, j)) /
                        std::max(std::abs(new_x(i, j)), std::abs(x(i, j))));
        }
      }
      std::cout << "Norm[" << iteration << "]: " << norm << std::endl;
      x = new_x;

      if (std::max(norm, prev_norm) < params.stopping_criterion) {
        converged = true;
        break;
      }

      prev_norm = norm;
    }

    if (converged) {
      // We have successfully converged!
      for (int i = 0; i < x.rows(); i++) {
        append_state(State(t[i], DependentVariablesFromMatrixRow<ODE>(x, i)));
      }

      // Set the current state to the final state we appended.
      current_state =
          State(t[x.rows() - 1],
                DependentVariablesFromMatrixRow<ODE>(x, x.rows() - 1));
    } else {
      // We failed to converge.
      return absl::Status(absl::StatusCode::kFailedPrecondition,
                          "Chebyshev-Picard iteration failed to converge.");
    }
  }

  return absl::OkStatus();
}

template <typename ODE_>
ChebyshevPicardIterator<ODE_> const&
ChebyshevPicardIterator<ODE_>::Instance::integrator() const {
  return integrator_;
}

template <typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ChebyshevPicardIterator<ODE_>::Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template <typename ODE_>
ChebyshevPicardIterator<ODE_>::Instance::Instance(
    InitialValueProblem<ODE> const& problem, AppendState const& append_state,
    Time const& step, ChebyshevPicardIterator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      integrator_(integrator) {}

template <typename ODE_>
ChebyshevPicardIterator<ODE_>::ChebyshevPicardIterator(
    const ChebyshevPicardIterationParams& params)
    : params_(params),
      nodes_(params.M + 1, uninitialized),
      cx_(params.M + 1, params.N + 1, uninitialized),
      cx_cα_(params.M + 1, params.N + 1, uninitialized) {
  // We use the notation from Macomber's thesis, section 1.4.3.
  const int M = params_.M;
  const int N = params_.N;
  CHECK_GE(M, 1);
  CHECK_GE(N, 1);

  // Populate nodes.
  for (int i = 0; i <= M; i++) {
    nodes_[i] = -Cos(π / M * i * Radian);
  }

  // βT is a (M + 1)×(N + 1) matrix of Chebyshev polynomials evaluated at nodes.
  // See Macomber's thesis, equation (1.20).
  UnboundedMatrix<double> βT(M + 1, N + 1, uninitialized);

  for (int i = 0; i <= M; i++) {
    const auto τᵢ = nodes_[i];
    // The 0-degree polynomial is uniformly 1.
    βT(i, 0) = 1;
    // The 0-degree polynomial is the identity.
    βT(i, 1) = τᵢ;

    // We populate the rest of βT using the recurrence relation.
    for (int j = 2; j <= N; j++) {
      βT(i, j) = 2 * τᵢ * βT(i, j - 1) - βT(i, j - 2);

      // Make sure the zeroes are actually zero.
      if (std::abs(βT(i, j)) < 1e-14) βT(i, j) = 0;
    }
  }

  // βW is a diagonal (N + 1)×(N + 1) matrix with diagonal [½, 1, 1, ..., ½].
  // See Macomber's thesis, equation (1.20).
  UnboundedMatrix<double> βW(N + 1, N + 1);
  βW(0, 0) = 0.5;
  βW(N, N) = 0.5;
  for (int i = 1; i < N; i++) {
    βW(i, i) = 1;
  }

  cx_ = βT * βW;

  // r is a diagonal (N + 1)×(N + 1) matrix.
  // See Macomber's thesis, equation (1.25).
  UnboundedMatrix<double> r(N + 1, N + 1);
  r(0, 0) = 1;
  r(N, N) = 1.0 / N;
  for (int i = 1; i < N; i++) {
    r(i, i) = 1.0 / (2 * i);
  }

  // s is an (N + 1)×N matrix.
  // See equation 1.26 in Macomber's thesis.
  UnboundedMatrix<double> s(N + 1, N);
  s(0, 0) = 1;
  s(0, 1) = -0.5;
  for (int k = 2; k < N; k++) {
    s(0, k) = (k % 2 == 1 ? 1 : -1) * (1.0 / (k - 1) - 1.0 / (k + 1));
  }
  for (int i = 0; i < N; i++) {
    s(i + 1, i) = 1;
  }
  for (int i = 1; i + 2 < N; i++) {
    s(i, i + 1) = -1;
  }

  // fT is βTᵀ with the last row removed.
  // See Macomber's thesis, equation (1.22).
  UnboundedMatrix<double> fT(N, M + 1, uninitialized);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j <= M; j++) {
      fT(i, j) = βT(j, i);
    }
  }

  // tV is 1/M * βW (we do not assign it to a variable).
  // cα is r * s * fT * tV (we do not assign it to a variable).

  cx_cα_ = 1.0 / M * cx_ * r * s * fT * βW;
}

template <typename ODE_>
ChebyshevPicardIterationParams const& ChebyshevPicardIterator<ODE_>::params()
    const {
  return params_;
}

template <typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ChebyshevPicardIterator<ODE_>::NewInstance(
    InitialValueProblem<ODE_> const& problem, AppendState const& append_state,
    Time const& step) const {
  // Cannot use `make_not_null_unique` because the constructor of `Instance` is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template <typename ODE_>
void ChebyshevPicardIterator<ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  // TODO(rnlahaye): implement me
}

}  // namespace internal
}  // namespace _chebyshev_picard_iterator
}  // namespace integrators
}  // namespace principia