
#pragma once

#include "numerics/matrix_computations.hpp"

#include "base/tags.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_matrix_computations {

using base::uninitialized;
using quantities::Abs;
using quantities::Pow;
using quantities::Sqrt;

// M may be a vanilla matrix or a lower-triangular matrix.
template<typename M>
void SwapRows(M& m, int const r1, int const r2) {
  for (int i = 0; i < M.columns(); ++i) {
    std::swap(m[r1][i], m[r2][i]);
  }
}

template<typename Scalar_, template<typename S> typename UpperTriangularMatrix>
struct CholeskyDecompositionGenerator<UpperTriangularMatrix<Scalar_>> {
  using Scalar = Scalar_;
  using Result = UpperTriangularMatrix<SquareRoot<Scalar>>;
  static Result Uninitialized(UpperTriangularMatrix<Scalar> const& u);
};

template<typename Scalar_, int columns,
         template<typename S, int c> typename UpperTriangularMatrix>
struct CholeskyDecompositionGenerator<UpperTriangularMatrix<Scalar_, columns>> {
  using Scalar = Scalar_;
  using Result = UpperTriangularMatrix<SquareRoot<Scalar>, columns>;
  static Result Uninitialized(UpperTriangularMatrix<Scalar, columns> const& u);
};

template<typename Scalar_,
         template<typename S> typename Vector,
         template<typename S> typename UpperTriangularMatrix>
struct ᵗRDRDecompositionGenerator<Vector<Scalar_>,
                                  UpperTriangularMatrix<Scalar_>> {
  using Scalar = Scalar_;
  struct Result {
    UpperTriangularMatrix<double> R;
    Vector<Scalar> D;
  };
  static Result Uninitialized(UpperTriangularMatrix<Scalar> const& u);
};

template<typename Scalar_, int columns,
         template<typename S, int c> typename Vector,
         template<typename S, int c> typename UpperTriangularMatrix>
struct ᵗRDRDecompositionGenerator<Vector<Scalar_, columns>,
                                  UpperTriangularMatrix<Scalar_, columns>> {
  using Scalar = Scalar_;
  struct Result {
    UpperTriangularMatrix<double, columns> R;
    Vector<Scalar, columns> D;
  };
  static Result Uninitialized(UpperTriangularMatrix<Scalar, columns> const& u);
};

template<typename LScalar, typename RScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
struct SubstitutionGenerator<Matrix<LScalar>, Vector<RScalar>> {
  using Result = Vector<Quotient<RScalar, LScalar>>;
  static Result Uninitialized(Matrix<LScalar> const& m);
};

template<typename LScalar, typename RScalar, int dimension,
         template<typename S, int d> typename Matrix,
         template<typename S, int d> typename Vector>
struct SubstitutionGenerator<Matrix<LScalar, dimension>,
                             Vector<RScalar, dimension>> {
  using Result = Vector<Quotient<RScalar, LScalar>, dimension>;
  static Result Uninitialized(Matrix<LScalar, dimension> const& m);
};

template<typename MScalar, typename VScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
struct RayleighQuotientGenerator<Matrix<MScalar>, Vector<VScalar>> {
  using Result = MScalar;
};

template<typename MScalar, typename VScalar, int dimension,
         template<typename S, int r, int c> typename Matrix,
         template<typename S, int d> typename Vector>
struct RayleighQuotientGenerator<Matrix<MScalar, dimension, dimension>,
                                 Vector<VScalar, dimension>> {
  using Result = MScalar;
};

template<typename MScalar, typename VScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
struct SolveGenerator<Matrix<MScalar>, Vector<VScalar>> {
  using Result = Vector<Quotient<VScalar, MScalar>>;
  static Result Uninitialized(Matrix<MScalar> const& m);
  static Matrix<double> UninitializedL(Matrix<MScalar> const& m);
  static Matrix<MScalar> UninitializedU(Matrix<MScalar> const& m);
};

template<typename MScalar, typename VScalar, int rows, int columns,
         template<typename S, int r, int c> typename Matrix,
         template<typename S, int d> typename Vector>
struct SolveGenerator<Matrix<MScalar, rows, columns>,
                      Vector<VScalar, rows>> {
  using Result = MScalar;
  static Result Uninitialized(
      Matrix<MScalar, rows, columns> const& m);
  static Matrix<double, rows, columns> UninitializedL(
      Matrix<MScalar, rows, columns> const& m);
  static Matrix<MScalar, columns, columns> UninitializedU(
      Matrix<MScalar, rows, columns> const& m);
};

template<typename Scalar_, template<typename S> typename UpperTriangularMatrix>
auto CholeskyDecompositionGenerator<UpperTriangularMatrix<Scalar_>>::
Uninitialized(UpperTriangularMatrix<Scalar> const& u) -> Result {
  return Result(u.columns(), uninitialized);
}

template<typename Scalar_, int columns,
         template<typename S, int c> typename UpperTriangularMatrix>
auto CholeskyDecompositionGenerator<UpperTriangularMatrix<Scalar_, columns>>::
Uninitialized(UpperTriangularMatrix<Scalar, columns> const& u) -> Result {
  return Result(uninitialized);
}

template<typename Scalar_,
         template<typename S> typename Vector,
         template<typename S> typename UpperTriangularMatrix>
auto ᵗRDRDecompositionGenerator<Vector<Scalar_>,
                                UpperTriangularMatrix<Scalar_>>::
Uninitialized(UpperTriangularMatrix<Scalar> const& u) -> Result {
  return {UpperTriangularMatrix<double>(u.columns(), uninitialized),
          Vector<Scalar>(u.columns(), uninitialized)};
}

template<typename Scalar_, int columns,
         template<typename S, int c> typename Vector,
         template<typename S, int c> typename UpperTriangularMatrix>
auto ᵗRDRDecompositionGenerator<Vector<Scalar_, columns>,
                                UpperTriangularMatrix<Scalar_, columns>>::
Uninitialized(UpperTriangularMatrix<Scalar, columns> const& u) -> Result {
  return {UpperTriangularMatrix<double, columns>(uninitialized),
          Vector<Scalar, columns>(uninitialized)};
}

template<typename LScalar, typename RScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
auto SubstitutionGenerator<Matrix<LScalar>, Vector<RScalar>>::Uninitialized(
    Matrix<LScalar> const& m) -> Result {
  return Result(m.columns(), uninitialized);
}

template<typename LScalar, typename RScalar, int dimension,
         template<typename S, int d> typename Matrix,
         template<typename S, int d> typename Vector>
auto SubstitutionGenerator<Matrix<LScalar, dimension>,
                           Vector<RScalar, dimension>>::Uninitialized(
    Matrix<LScalar, dimension> const& m) -> Result {
  return Result(uninitialized);
}

template<typename MScalar, typename VScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
auto SolveGenerator<Matrix<MScalar>, Vector<VScalar>>::Uninitialized(
    Matrix<MScalar> const& m) -> Result {
  return Result(m.columns(), uninitialized);
}

template<typename MScalar, typename VScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
auto SolveGenerator<Matrix<MScalar>, Vector<VScalar>>::UninitializedL(
    Matrix<MScalar> const& m) -> Matrix<double> {
  return Matrix<double>(m.rows(), m.columns(), uninitialized);
}

template<typename MScalar, typename VScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
auto SolveGenerator<Matrix<MScalar>, Vector<VScalar>>::UninitializedU(
    Matrix<MScalar> const& m) -> Matrix<MScalar> {
  return Matrix<MScalar>(m.columns(), m.columns(), uninitialized);
}

template<typename MScalar, typename VScalar, int rows, int columns,
         template<typename S, int r, int c> typename Matrix,
         template<typename S, int d> typename Vector>
auto SolveGenerator<Matrix<MScalar, rows, columns>,
                    Vector<VScalar, rows>>::Uninitialized(
    Matrix<MScalar, rows, columns> const& m) -> Result {
  return Result(uninitialized);
}

template<typename MScalar, typename VScalar, int rows, int columns,
         template<typename S, int r, int c> typename Matrix,
         template<typename S, int d> typename Vector>
auto SolveGenerator<Matrix<MScalar, rows, columns>,
                    Vector<VScalar, rows>>::UninitializedL(
    Matrix<MScalar, rows, columns> const& m) -> Matrix<double, rows, columns> {
  return Matrix<double, rows, columns>(uninitialized);
}

template<typename MScalar, typename VScalar, int rows, int columns,
         template<typename S, int r, int c> typename Matrix,
         template<typename S, int d> typename Vector>
auto SolveGenerator<Matrix<MScalar, rows, columns>,
                    Vector<VScalar, rows>>::UninitializedU(
   Matrix<MScalar, rows, columns> const& m)
        -> Matrix<MScalar, columns, columns> {
  return Matrix<MScalar, columns, columns>(uninitialized);
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
typename SolveGenerator<Matrix, Vector>::Result
Solve(Matrix const& A, Vector const& b) {
  // This implementation follows [Hig02].
  using G = LUDecompositionGenerator<Matrix,
                                     LowerTriangularMatrix,
                                     UpperTriangularMatrix>;
  using Scalar = typename G::Scalar;
  auto x = G::Uninitialized(A);

  // The units make it inconvenient to overlay L and U onto A.
  Matrix L = G::UninitializedL(A);
  Matrix U = G::UninitializedU(A);

  // Doolittle's method: write P * A = L * U where P is an implicit permutation
  // that is also applied to b.  See [Hig02], Algorithm 9.2 p. 162.
  for (int k = 0; k < A.columns(); ++k) {
    // Partial pivoting.
    int r = -1;
    Scalar max{};
    for (int i = k; i < A.rows(); ++i) {
      if (Abs(A[i][k]) > max) {
        r = i;
        max = Abs(A[i][k]));
      }
    }
    CHECK_LE(0, r) << A << " cannot pivot";
    CHECK_GT(A.size(), r) << A << " cannot pivot";
    SwapRows(A, k, r);
    SwapRows(L, k, r);
    std::swap(b[k], b[r]);
    CHECK_NE(Scalar{}, A[k][k])) << *this << " is singular";

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
}

}  // namespace internal_matrix_computations
}  // namespace numerics
}  // namespace principia
