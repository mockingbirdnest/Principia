
#pragma once

#include "numerics/matrix_computations.hpp"

#include <algorithm>

#include "base/tags.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_matrix_computations {

using base::uninitialized;
using quantities::Abs;
using quantities::Pow;
using quantities::Sqrt;

template<typename Scalar_>
struct CholeskyDecompositionGenerator<UnboundedUpperTriangularMatrix<Scalar_>> {
  using Scalar = Scalar_;
  using Result = UnboundedUpperTriangularMatrix<SquareRoot<Scalar>>;
  static Result Uninitialized(UnboundedUpperTriangularMatrix<Scalar> const& u);
};

template<typename Scalar_, int columns>
struct CholeskyDecompositionGenerator<
    FixedUpperTriangularMatrix<Scalar_, columns>> {
  using Scalar = Scalar_;
  using Result = FixedUpperTriangularMatrix<SquareRoot<Scalar>, columns>;
  static Result Uninitialized(
      FixedUpperTriangularMatrix<Scalar, columns> const& u);
};

template<typename Scalar_>
struct ᵗRDRDecompositionGenerator<UnboundedVector<Scalar_>,
                                  UnboundedUpperTriangularMatrix<Scalar_>> {
  using Scalar = Scalar_;
  struct Result {
    UnboundedUpperTriangularMatrix<double> R;
    UnboundedVector<Scalar> D;
  };
  static Result Uninitialized(UnboundedUpperTriangularMatrix<Scalar> const& u);
};

template<typename Scalar_, int columns>
struct ᵗRDRDecompositionGenerator<
    FixedVector<Scalar_, columns>,
    FixedUpperTriangularMatrix<Scalar_, columns>> {
  using Scalar = Scalar_;
  struct Result {
    FixedUpperTriangularMatrix<double, columns> R;
    FixedVector<Scalar, columns> D;
  };
  static Result Uninitialized(
      FixedUpperTriangularMatrix<Scalar, columns> const& u);
};

template<typename LScalar, typename RScalar,
         template<typename S> typename TriangularMatrix>
struct SubstitutionGenerator<TriangularMatrix<LScalar>,
                             UnboundedVector<RScalar>> {
  using Result = UnboundedVector<Quotient<RScalar, LScalar>>;
  static Result Uninitialized(TriangularMatrix<LScalar> const& m);
};

template<typename LScalar, typename RScalar, int dimension,
         template<typename S, int d> typename TriangularMatrix>
struct SubstitutionGenerator<TriangularMatrix<LScalar, dimension>,
                             FixedVector<RScalar, dimension>> {
  using Result = FixedVector<Quotient<RScalar, LScalar>, dimension>;
  static Result Uninitialized(TriangularMatrix<LScalar, dimension> const& m);
};

template<typename MScalar, typename VScalar>
struct RayleighQuotientGenerator<UnboundedMatrix<MScalar>,
                                 UnboundedVector<VScalar>> {
  using Result = MScalar;
};

template<typename MScalar, typename VScalar, int dimension>
struct RayleighQuotientGenerator<FixedMatrix<MScalar, dimension, dimension>,
                                 FixedVector<VScalar, dimension>> {
  using Result = MScalar;
};

template<typename MScalar, typename VScalar>
struct SolveGenerator<UnboundedMatrix<MScalar>, UnboundedVector<VScalar>> {
  using Scalar = MScalar;
  using Result = UnboundedVector<Quotient<VScalar, MScalar>>;
  static UnboundedLowerTriangularMatrix<double> UninitializedL(
      UnboundedMatrix<MScalar> const& m);
  static UnboundedUpperTriangularMatrix<MScalar> UninitializedU(
      UnboundedMatrix<MScalar> const& m);
};

template<typename MScalar, typename VScalar, int rows, int columns>
struct SolveGenerator<FixedMatrix<MScalar, rows, columns>,
                      FixedVector<VScalar, rows>> {
  using Scalar = MScalar;
  using Result = FixedVector<Quotient<VScalar, MScalar>, columns>;
  static FixedLowerTriangularMatrix<double, rows> UninitializedL(
      FixedMatrix<MScalar, rows, columns> const& m);
  static FixedUpperTriangularMatrix<MScalar, columns> UninitializedU(
      FixedMatrix<MScalar, rows, columns> const& m);
};

template<typename Scalar_>
auto CholeskyDecompositionGenerator<UnboundedUpperTriangularMatrix<Scalar_>>::
Uninitialized(UnboundedUpperTriangularMatrix<Scalar> const& u) -> Result {
  return Result(u.columns(), uninitialized);
}

template<typename Scalar_, int columns>
auto CholeskyDecompositionGenerator<
    FixedUpperTriangularMatrix<Scalar_, columns>>::
Uninitialized(FixedUpperTriangularMatrix<Scalar, columns> const& u) -> Result {
  return Result(uninitialized);
}

template<typename Scalar_>
auto ᵗRDRDecompositionGenerator<UnboundedVector<Scalar_>,
                                UnboundedUpperTriangularMatrix<Scalar_>>::
Uninitialized(UnboundedUpperTriangularMatrix<Scalar> const& u) -> Result {
  return {UnboundedUpperTriangularMatrix<double>(u.columns(), uninitialized),
          UnboundedVector<Scalar>(u.columns(), uninitialized)};
}

template<typename Scalar_, int columns>
auto ᵗRDRDecompositionGenerator<FixedVector<Scalar_, columns>,
                                FixedUpperTriangularMatrix<Scalar_, columns>>::
Uninitialized(FixedUpperTriangularMatrix<Scalar, columns> const& u) -> Result {
  return {FixedUpperTriangularMatrix<double, columns>(uninitialized),
          FixedVector<Scalar, columns>(uninitialized)};
}

template<typename LScalar, typename RScalar,
         template<typename S> typename TriangularMatrix>
auto SubstitutionGenerator<TriangularMatrix<LScalar>,
                           UnboundedVector<RScalar>>::
Uninitialized(TriangularMatrix<LScalar> const& m) -> Result {
  return Result(m.columns(), uninitialized);
}

template<typename LScalar, typename RScalar, int dimension,
         template<typename S, int d> typename TriangularMatrix>
auto SubstitutionGenerator<TriangularMatrix<LScalar, dimension>,
                           FixedVector<RScalar, dimension>>::Uninitialized(
    TriangularMatrix<LScalar, dimension> const& m) -> Result {
  return Result(uninitialized);
}

template<typename MScalar, typename VScalar>
UnboundedLowerTriangularMatrix<double>
SolveGenerator<UnboundedMatrix<MScalar>, UnboundedVector<VScalar>>::
UninitializedL(UnboundedMatrix<MScalar> const& m) {
  return UnboundedLowerTriangularMatrix<double>(m.rows(), uninitialized);
}

template<typename MScalar, typename VScalar>
UnboundedUpperTriangularMatrix<MScalar>
SolveGenerator<UnboundedMatrix<MScalar>, UnboundedVector<VScalar>>::
UninitializedU(UnboundedMatrix<MScalar> const& m) {
  return UnboundedUpperTriangularMatrix<MScalar>(m.columns(), uninitialized);
}

template<typename MScalar, typename VScalar, int rows, int columns>
FixedLowerTriangularMatrix<double, rows>
SolveGenerator<FixedMatrix<MScalar, rows, columns>,
               FixedVector<VScalar, rows>>::
UninitializedL(FixedMatrix<MScalar, rows, columns> const& m) {
  return FixedLowerTriangularMatrix<double, rows>(uninitialized);
}

template<typename MScalar, typename VScalar, int rows, int columns>
FixedUpperTriangularMatrix<MScalar, columns>
SolveGenerator<FixedMatrix<MScalar, rows, columns>,
               FixedVector<VScalar, rows>>::
UninitializedU(FixedMatrix<MScalar, rows, columns> const& m) {
  return FixedUpperTriangularMatrix<MScalar, columns>(uninitialized);
}

// [Hig02], Algorithm 10.2.
template<typename UpperTriangularMatrix>
typename CholeskyDecompositionGenerator<UpperTriangularMatrix>::Result
CholeskyDecomposition(UpperTriangularMatrix const& A) {
  using G = CholeskyDecompositionGenerator<UpperTriangularMatrix>;
  using Scalar = typename G::Scalar;
  auto R = G::Uninitialized(A);
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
typename ᵗRDRDecompositionGenerator<Vector, UpperTriangularMatrix>::Result
ᵗRDRDecomposition(UpperTriangularMatrix const& A) {
  using G = ᵗRDRDecompositionGenerator<Vector, UpperTriangularMatrix>;
  using Scalar = typename G::Scalar;
  auto result = G::Uninitialized(A);
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
  using G = SubstitutionGenerator<UpperTriangularMatrix, Vector>;
  auto x = G::Uninitialized(U);
  int const n = x.size() - 1;
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
  using G = SubstitutionGenerator<LowerTriangularMatrix, Vector>;
  auto x = G::Uninitialized(L);
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

template<typename Matrix, typename Vector>
typename RayleighQuotientGenerator<Matrix, Vector>::Result RayleighQuotient(
    Matrix const& A, Vector const& x) {
  // [GV13], section 8.2.3.
  return x.Transpose() * (A * x) / (x.Transpose() * x);
}

template<typename Matrix, typename Vector>
typename RayleighQuotientIterationGenerator<Matrix, Vector>::Result
RayleighQuotientIteration(Matrix const& A, Vector const& x) {
  auto xₖ = x;
  for (;;) {
    auto const μₖ = RayleighQuotient(A, xₖ);
    auto const A_minus_μₖ_I = A;
    for (int i = 0; i < A.size(); ++i) {
      A[i][i] -= μₖ;
    }
    auto const zₖ₊₁ = Solve(A_minus_μₖ_I, xₖ);
    xₖ = zₖ₊₁ / zₖ₊₁.Norm();
  }
  return xₖ;
}


template<typename Matrix, typename Vector>
typename SolveGenerator<Matrix, Vector>::Result
Solve(Matrix A, Vector b) {
  // This implementation follows [Hig02].
  using G = SolveGenerator<Matrix, Vector>;
  using Scalar = typename G::Scalar;

  // The units make it inconvenient to overlay L and U onto A.
  auto L = G::UninitializedL(A);
  auto U = G::UninitializedU(A);

  // Doolittle's method: write P * A = L * U where P is an implicit permutation
  // that is also applied to b.  See [Hig02], Algorithm 9.2 p. 162.
  for (int k = 0; k < A.columns(); ++k) {
    // Partial pivoting.
    int r = -1;
    Scalar max{};
    for (int i = k; i < A.rows(); ++i) {
      if (Abs(A[i][k]) > max) {
        r = i;
        max = Abs(A[i][k]);
      }
    }
    CHECK_LE(0, r) << A << " cannot pivot";
    CHECK_LT(r, A.rows()) << A << " cannot pivot";

    // Swap the rows of A.
    for (int i = 0; i < A.columns(); ++i) {
      std::swap(A[k][i], A[r][i]);
    }

    // Swap the rows of L.
    for (int i = 0; i < k; ++i) {
      std::swap(L[k][i], L[r][i]);
    }

    std::swap(b[k], b[r]);
    CHECK_NE(Scalar{}, A[k][k]) << A << " is singular";

    for (int j = k; j < A.columns(); ++j) {
      auto U_kj = A[k][j];
      for (int i = 0; i < k; ++i) {
        U_kj -= L[k][i] * U[i][j];
      }
      U[k][j] = U_kj;
    }
    for (int i = k + 1; i < A.rows(); ++i) {
      auto L_ik = A[i][k];
      for (int j = 0; j < k; ++j) {
        L_ik -= L[i][j] * U[j][k];
      }
      L[i][k] = L_ik / U[k][k];
    }
    L[k][k] = 1;
  }

  // For the resolution of triangular systems see [Hig02], Algorithm 8.1 p. 140.

  // Find y such that L * y = P * b.
  auto const y = ForwardSubstitution(L, b);
  // Find x such that U * x = y.
  auto const x = BackSubstitution(U, y);

  return x;
}

}  // namespace internal_matrix_computations
}  // namespace numerics
}  // namespace principia
