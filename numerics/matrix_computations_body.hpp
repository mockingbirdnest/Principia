
#pragma once

#include "numerics/matrix_computations.hpp"

#include "base/tags.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_arrays {

using base::uninitialized;
using quantities::Pow;
using quantities::Sqrt;

// [Hig02], Algorithm 10.2.
template<typename UpperTriangularMatrix>
typename CholeskyGenerator<UpperTriangularMatrix>::Result
CholeskyDecomposition(UpperTriangularMatrix const& A) {
  using G = CholeskyGenerator<UpperTriangularMatrix>;
  using Scalar = G::Scalar;
  typename G::Result R(A.columns(), uninitialized);
  for (int j = 0; j < A.columns(); ++j) {
    for (int i = 0; i < j; ++i) {
      Scalar Σrₖᵢrₖⱼ{};
      for (int k = 0; k < i; ++k) {
        Σrₖᵢrₖⱼ += R[k][i] * R[k][j];
      }
      R[i][j] = (A[i][j] - Σrₖᵢrₖⱼ) / R[i][i];
    }
    Scalar Σrₖⱼ²{};
    for (int k = 0; k < j; ++k) {
      Σrₖⱼ² += Pow<2>(R[k][j]);
    }
    // This will produce NaNs if the matrix is not positive definite.
    R[j][j] = Sqrt(A[j][j] - Σrₖⱼ²);
  }
  return R;
}

// [KM13], formulæ (10) and (11).
template<typename Vector, typename UpperTriangularMatrix>
typename ᵗRDRGenerator<Vector, UpperTriangularMatrix>::Result
ᵗRDRDecomposition(UpperTriangularMatrix const& A) {
  using G = ᵗRDRGenerator<Vector, UpperTriangularMatrix>;
  using Scalar = G::Scalar;
  typename G::Result result;
  auto& R = result.R;
  auto& D = result.D;
  for (int i = 0; i < A.columns(); ++i) {
    Scalar Σrₖᵢ²dₖ{};
    for (int k = 0; k < i; ++k) {
      Σrₖᵢ²dₖ += Pow<2>(R[k][i]) * D[k];
    }
    D[i] = A[i][i] - Σrₖᵢ²dₖ;
    for (int j = i + 1; j < A.columns(); ++j) {
      Scalar Σrₖᵢrₖⱼdₖ{};
      for (int k = 0; k < i; ++k) {
        Σrₖᵢrₖⱼdₖ += R[k][i] * R[k][j] * D[k];
      }
      R[i][j] = (A[i][j] - Σrₖᵢrₖⱼdₖ) / D[i];
    }
    R[i][i] = 1;
  }
  return result;
}

// [Hig02], Algorithm 8.1.
template<typename UpperTriangularMatrix, typename Vector>
typename SubstitutionGenerator<UpperTriangularMatrix, Vector>::Result
BackSubstitution(UpperTriangularMatrix const& U,
                 Vector const& b) {
  typename SubstitutionGenerator<UpperTriangularMatrix, Vector>::Result
  x(b.size(), uninitialized);
  int const n = b.size() - 1;
  x[n] = b[n] / U[n][n];
  for (int i = n - 1; i >= 0; --i) {
    auto s = b[i];
    for (int j = i + 1; j <= n; ++j) {
      s -= U[i][j] * x[j];
    }
    x[i] = s / U[i][i];
  }
  return x;
}

// [Hig02] says: "We will not state the analogous algorithm for solving a lower
// triangular system, forward substitution."  So we follow
// https://en.wikipedia.org/wiki/Triangular_matrix#Forward_substitution.
template<typename LowerTriangularMatrix, typename Vector>
typename SubstitutionGenerator<LowerTriangularMatrix, Vector>::Result
ForwardSubstitution(LowerTriangularMatrix const& L,
                    Vector const& b) {
  typename SubstitutionGenerator<LowerTriangularMatrix, Vector>::Result
  x(b.size(), uninitialized);
  x[0] = b[0] / L[0][0];
  for (int i = 1; i < b.size(); ++i) {
    auto s = b[i];
    for (int j = 0; j < i; ++j) {
      s -= L[i][j] * x[j];
    }
    x[i] = s / L[i][i];
  }
  return x;
}

}  // namespace internal_arrays
}  // namespace numerics
}  // namespace principia
