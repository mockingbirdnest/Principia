
#pragma once

#include "numerics/matrix_computations.hpp"

#include "base/tags.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_matrix_computations {

using base::uninitialized;
using quantities::Pow;
using quantities::Sqrt;

template<typename Scalar_, template<typename S> typename UpperTriangularMatrix>
struct CholeskyGenerator<UpperTriangularMatrix<Scalar_>> {
  using Scalar = Scalar_;
  using Result = UpperTriangularMatrix<SquareRoot<Scalar>>;
  static Result Uninitialized(UpperTriangularMatrix<Scalar> const& u);
};

template<typename Scalar_, int columns,
         template<typename S, int c> typename UpperTriangularMatrix>
struct CholeskyGenerator<UpperTriangularMatrix<Scalar_, columns>> {
  using Scalar = Scalar_;
  using Result = UpperTriangularMatrix<SquareRoot<Scalar>, columns>;
  static Result Uninitialized(UpperTriangularMatrix<Scalar, columns> const& u);
};

template<typename Scalar_,
         template<typename S> typename Vector,
         template<typename S> typename UpperTriangularMatrix>
struct ᵗRDRGenerator<Vector<Scalar_>, UpperTriangularMatrix<Scalar_>> {
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
struct ᵗRDRGenerator<Vector<Scalar_, columns>,
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

};


template<typename Scalar_, template<typename S> typename UpperTriangularMatrix>
auto CholeskyGenerator<UpperTriangularMatrix<Scalar_>>::Uninitialized(
    UpperTriangularMatrix<Scalar> const& u) -> Result {
  return Result(u.columns(), uninitialized);
}

template<typename Scalar_, int columns,
         template<typename S, int c> typename UpperTriangularMatrix>
auto CholeskyGenerator<UpperTriangularMatrix<Scalar_, columns>>::Uninitialized(
    UpperTriangularMatrix<Scalar, columns> const& u) -> Result {
  return Result(uninitialized);
}

template<typename Scalar_,
         template<typename S> typename Vector,
         template<typename S> typename UpperTriangularMatrix>
auto ᵗRDRGenerator<Vector<Scalar_>, UpperTriangularMatrix<Scalar_>>::
Uninitialized(UpperTriangularMatrix<Scalar> const& u) -> Result {
  return {UpperTriangularMatrix<double>(u.columns(), uninitialized),
          Vector<Scalar>(u.columns(), uninitialized)};
}

template<typename Scalar_, int columns,
         template<typename S, int c> typename Vector,
         template<typename S, int c> typename UpperTriangularMatrix>
auto ᵗRDRGenerator<Vector<Scalar_, columns>,
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

// [Hig02], Algorithm 10.2.
template<typename UpperTriangularMatrix>
typename CholeskyGenerator<UpperTriangularMatrix>::Result
CholeskyDecomposition(UpperTriangularMatrix const& A) {
  using G = CholeskyGenerator<UpperTriangularMatrix>;
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
typename ᵗRDRGenerator<Vector, UpperTriangularMatrix>::Result
ᵗRDRDecomposition(UpperTriangularMatrix const& A) {
  using G = ᵗRDRGenerator<Vector, UpperTriangularMatrix>;
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
typename RayleighQuotientGenerator<Matrix, Vector> RayleighQuotient(
    Matrix const& A, Vector const& x) {
  // [GV13], section 8.2.3.
  return x * (A * x) / (x * x);
}

}  // namespace internal_matrix_computations
}  // namespace numerics
}  // namespace principia
