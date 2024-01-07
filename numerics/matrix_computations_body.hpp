#pragma once

#include "numerics/matrix_computations.hpp"

#include <algorithm>
#include <limits>
#include <utility>

#include "base/tags.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/transposed_view.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _matrix_computations {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

// TODO(phl): The view stuff should be (1) made completed, i.e., have all the
// operations that exist for fixed/unbounded vectors/matrices; (2) moved to a
// common place (probably together with TransposedView); (3) unified with
// fixed/unbounded arrays so that we don't have to write each algorithm N times.

// A view of a column of a matrix.
// TODO(phl): Could we extract Scalar from Matrix?
template<typename Scalar, typename Matrix>
struct ColumnView {
  Matrix& matrix;
  int first_row;
  int last_row;
  int column;

  Scalar Norm() const;
  Square<Scalar> Norm²() const;
  constexpr int size() const;

  // Constructs an unbounded vector by copying data from the view.  Note that
  // the result is unbounded even if the matrix being viewed is a FixedMatrix.
  explicit operator UnboundedVector<Scalar>() const;

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  ColumnView& operator/=(double right);
};

template<typename Scalar, typename Matrix>
Square<Scalar> ColumnView<Scalar, Matrix>::Norm²() const {
  Square<Scalar> result{};
  for (int i = first_row; i <= last_row; ++i) {
    result += Pow<2>(matrix(i, column));
  }
  return result;
}

template<typename Scalar, typename Matrix>
Scalar ColumnView<Scalar, Matrix>::Norm() const {
  return Sqrt(Norm²());
}

template<typename Scalar, typename Matrix>
constexpr int ColumnView<Scalar, Matrix>::size() const {
  return last_row - first_row + 1;
}

template<typename Scalar, typename Matrix>
ColumnView<Scalar, Matrix>::operator UnboundedVector<Scalar>() const {
  UnboundedVector<Scalar> result(size(), uninitialized);
  for (int i = first_row; i <= last_row; ++i) {
    result[i - first_row] = matrix(i, column);
  }
  return result;
}

template<typename Scalar, typename Matrix>
constexpr Scalar& ColumnView<Scalar, Matrix>::operator[](int const index) {
  CONSTEXPR_DCHECK(index <= last_row - first_row);
  return matrix(first_row + index, column);
}

template<typename Scalar, typename Matrix>
constexpr Scalar const& ColumnView<Scalar, Matrix>::operator[](
    int const index) const {
  CONSTEXPR_DCHECK(index <= last_row - first_row);
  return matrix(first_row + index, column);
}

template<typename Scalar, typename Matrix>
ColumnView<Scalar, Matrix>& ColumnView<Scalar, Matrix>::operator/=(
    double const right) {
  for (int i = first_row; i < last_row; ++i) {
    matrix(i, column) /= right;
  }
}

template<typename Scalar, typename Matrix>
UnboundedVector<double> Normalize(ColumnView<Scalar, Matrix> const& view) {
  return UnboundedVector<Scalar>(view) / view.Norm();
}

template<typename Scalar, typename Matrix>
std::ostream& operator<<(std::ostream& out,
                         ColumnView<Scalar, Matrix> const& view) {
  std::stringstream s;
  for (int i = 0; i < view.size(); ++i) {
    s << (i == 0 ? "{" : "") << view[i]
      << (i == view.size() - 1 ? "}" : ", ");
  }
  out << s.str();
  return out;
}

template<typename Scalar, typename Matrix>
struct BlockView {
  Matrix& matrix;
  int first_row;
  int last_row;
  int first_column;
  int last_column;

  constexpr int rows() const;
  constexpr int columns() const;

  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  BlockView& operator-=(UnboundedMatrix<Scalar> const& right);
};

template<typename Scalar, typename Matrix>
constexpr int BlockView<Scalar, Matrix>::rows() const {
  return last_row - first_row + 1;
}

template<typename Scalar, typename Matrix>
constexpr int BlockView<Scalar, Matrix>::columns() const {
  return last_column - first_column + 1;
}

template<typename Scalar, typename Matrix>
constexpr Scalar& BlockView<Scalar, Matrix>::operator()(int const row,
                                                        int const column) {
  CONSTEXPR_DCHECK(row <= last_row - first_row);
  CONSTEXPR_DCHECK(column <= last_column - first_column);
  return matrix(first_row + row, first_column + column);
}

template<typename Scalar, typename Matrix>
constexpr Scalar const& BlockView<Scalar, Matrix>::operator()(
    int const row,
    int const column) const {
  CONSTEXPR_DCHECK(row <= last_row - first_row);
  CONSTEXPR_DCHECK(column <= last_column - first_column);
  return matrix(first_row + row, first_column + column);
}

template<typename Scalar, typename Matrix>
BlockView<Scalar, Matrix>& BlockView<Scalar, Matrix>::operator-=(
    UnboundedMatrix<Scalar> const& right) {
  CHECK_EQ(rows(), right.rows());
  CHECK_EQ(columns(), right.columns());
  for (int i = 0; i < right.rows(); ++i) {
    for (int j = 0; j < right.columns(); ++j) {
      matrix(first_row + i, first_column + j) -= right(i, j);
    }
  }
  return *this;
}


template<typename LScalar, typename RScalar, typename Matrix>
UnboundedVector<Product<LScalar, RScalar>> operator*(
    BlockView<LScalar, Matrix> const& left,
    UnboundedVector<RScalar> const& right) {
  CHECK_EQ(left.columns(), right.size());
  UnboundedVector<Product<LScalar, RScalar>> result(left.rows());
  for (int i = 0; i < left.rows(); ++i) {
    auto& result_i = result[i];
    for (int j = 0; j < left.columns(); ++j) {
      result_i += left(i, j) * right[j];
    }
  }
  return result;
}

template<typename LScalar, typename RScalar, typename Matrix>
UnboundedVector<Product<LScalar, RScalar>> operator*(
    TransposedView<BlockView<LScalar, Matrix>> const& left,
    UnboundedVector<RScalar> const& right) {
  CHECK_EQ(left.transpose.rows(), right.size());
  UnboundedVector<Product<LScalar, RScalar>> result(left.transpose.columns());
  for (int j = 0; j < left.transpose.columns(); ++j) {
    auto& result_j = result[j];
    for (int i = 0; i < left.transpose.rows(); ++i) {
      result_j += left.transpose(i, j) * right[i];
    }
  }
  return result;
}

// As mentioned in [GV13] section 5.1.4, "It is critical to exploit structure
// when applying [the Householder reflection] to a matrix".
struct HouseholderReflection {
  UnboundedVector<double> v;
  double β;
};

template<typename Vector>
HouseholderReflection ComputeHouseholderReflection(Vector const& x) {
  // In order to avoid issues with quantities, we start by normalizing x.  This
  // implies that μ is 1.
  auto const normalized_x = Normalize(x);
  HouseholderReflection result{.v = normalized_x, .β = 0};
  double const& x₁ = normalized_x[0];
  double& v₁ = result.v[0];
  auto x₂ₘ = normalized_x;
  x₂ₘ[0] = 0;
  auto const σ = x₂ₘ.Norm²();
  v₁ = 1;
  if (σ == 0) {
    if (x₁ < 0) {
      result.β = -2;
    }
  } else {
    static constexpr double μ = 1;
    if (x₁ <= 0) {
      v₁ = x₁ - μ;
    } else {
      v₁ = -σ / (x₁ + μ);
    }
    double const v₁² = Pow<2>(v₁);
    result.β = 2 * v₁² / (σ + v₁²);
    result.v /= v₁;
  }
  return result;
}

// A becomes P A.
template<typename Matrix>
void Premultiply(HouseholderReflection const& P, Matrix& A) {
  // We don't have a multiplication TransposedView<Vector> * Matrix because the
  // ownership of the result is problematic.  Instead, we transpose twice.  That
  // costs essentially nothing.
  auto const ᵗAv = TransposedView{A} * P.v;  // NOLINT
  auto const ᵗvA = TransposedView{ᵗAv};      // NOLINT
  auto const βv = P.β * P.v;
  A -= βv * ᵗvA;
}

// A becomes A P.
template<typename Matrix>
void PostMultiply(Matrix& A, HouseholderReflection const& P) {
  auto const βv = P.β * P.v;
  auto const βᵗv = TransposedView{βv};  // NOLINT
  auto const Av = A * P.v;
  A -= Av * βᵗv;
}

// This is J(p, q, θ) in [GV13] section 8.5.1.  This matrix is also called a
// Givens rotation.  As mentioned in [GV13] section 5.1.9, "It is critical that
// the special structure of a Givens rotation matrix be exploited".
struct JacobiRotation {
  double cos;
  double sin;
  int p;
  int q;
};

// See [GV13] section 8.5.2, algorithm 8.5.1.
template<typename Scalar, typename Matrix>
JacobiRotation SymmetricSchurDecomposition2By2(Matrix const& A,
                                               int const p,
                                               int const q) {
  DCHECK_LE(0, p);
  DCHECK_LT(p, q);
  DCHECK_LT(q, A.rows());
  constexpr Scalar const zero{};
  JacobiRotation J = {.p = p, .q = q};
  auto& c = J.cos;
  auto& s = J.sin;
  if (A(p, q) != zero) {
    double const τ = (A(q, q) - A(p, p)) / (2 * A(p, q));
    double t;
    if (τ >= 0) {
      t = 1 / (τ + Sqrt(1 + τ * τ));
    } else {
      t = 1 / (τ - Sqrt(1 + τ * τ));
    }
    c = 1 / Sqrt(1 + t * t);
    s = t * c;
  } else {
    c = 1;
    s = 0;
  }
  return J;
};

// For these two functions, see [GV13] section 5.1.9.

// A becomes ᵗJ A.
template<typename Matrix>
void PremultiplyByTranspose(JacobiRotation const& J, Matrix& A) {
  auto const& [c, s, p, q] = J;
  for (int j = 0; j < A.columns(); ++j) {
    auto const τ₁ = A(p, j);
    auto const τ₂ = A(q, j);
    A(p, j) = c * τ₁ - s * τ₂;
    A(q, j) = s * τ₁ + c * τ₂;
  }
}

// A becomes A J.
template<typename Matrix>
void PostMultiply(Matrix& A, JacobiRotation const& J) {
  auto const& [c, s, p, q] = J;
  for (int j = 0; j < A.rows(); ++j) {
    auto const τ₁ = A(j, p);
    auto const τ₂ = A(j, q);
    A(j, p) = c * τ₁ - s * τ₂;
    A(j, q) = s * τ₁ + c * τ₂;
  }
}

// [GV13] algorithm 7.5.1.
template<typename Scalar, typename Matrix>
void FrancisQRStep(Matrix& H) {
  int const n = H.rows();
  int const m = n - 1;
  auto const s = H(m - 1, m - 1) + H(n - 1, n - 1);
  auto const t = H(m - 1, m - 1) * H(n - 1, n - 1) -
                 H(m - 1, n - 1) * H(n - 1, m - 1);
  FixedVector<Scalar, 3> xyz(uninitialized);
  auto& x = xyz[0];
  auto& y = xyz[1];
  auto& z = xyz[2];
  x = Pow<2>(H(0, 0)) + H(0, 1) * H(1, 0) - s * H(0, 0) + t;
  y = H(1, 0) * (H(0, 0) + H(1, 1) - s);
  z = H(1, 0) * H(2, 1);
  for (int k = 0; k < n - 2; ++k) {
    auto const P = ComputeHouseholderReflection(xyz);
    int const q = std::max(1, k);
    {
      auto block = BlockView<Scalar, Matrix>{.matrix = H,
                                             .first_row = k,
                                             .last_row = k + 2,
                                             .first_column = q - 1,
                                             .last_column = n - 1};
      Premultiply(P, block);
    }
    int const r = std::min(k + 4, n);
    {
      auto block = BlockView<Scalar, Matrix>{.matrix = H,
                                             .first_row = 0,
                                             .last_row = r - 1,
                                             .first_column = k,
                                             .last_column = k + 2};
      PostMultiply(block, P);
    }
    x = H(k + 1, k);
    y = H(k + 2, k);
    if (k < n - 3) {
      z = H(k + 3, k);
    }
  }
  FixedVector<Scalar, 2> xy({x, y});
  auto const P = ComputeHouseholderReflection(xy);
  {
    auto block = BlockView<Scalar, Matrix>{.matrix = H,
                                           .first_row = n - 2,
                                           .last_row = n - 1,
                                           .first_column = n - 3,
                                           .last_column = n - 1};
    Premultiply(P, block);
  }
  {
    auto block = BlockView<Scalar, Matrix>{.matrix = H,
                                           .first_row = 0,
                                           .last_row = n - 1,
                                           .first_column = n - 2,
                                           .last_column = n - 1};
    PostMultiply(block, P);
  }
  // TODO(phl): Accumulate Z.
}


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

template<typename Scalar_>
struct HessenbergDecompositionGenerator<UnboundedMatrix<Scalar_>> {
  using Scalar = Scalar_;
  struct Result {
    UnboundedMatrix<Scalar> H;
  };
};

template<typename Scalar_, int dimension>
struct HessenbergDecompositionGenerator<
    FixedMatrix<Scalar_, dimension, dimension>> {
  using Scalar = Scalar_;
  struct Result {
    FixedMatrix<Scalar, dimension, dimension> H;
  };
};

template<typename Scalar_>
struct QRDecompositionGenerator<UnboundedMatrix<Scalar_>> {
  using Scalar = Scalar_;
  struct Result {
    UnboundedMatrix<Scalar> R;
  };
};

template<typename Scalar_, int dimension>
struct QRDecompositionGenerator<
    FixedMatrix<Scalar_, dimension, dimension>> {
  using Scalar = Scalar_;
  struct Result {
    FixedMatrix<Scalar, dimension, dimension> R;
  };
};

template<typename Scalar_>
struct ClassicalJacobiGenerator<UnboundedMatrix<Scalar_>> {
  using Scalar = Scalar_;
  using Rotation = UnboundedMatrix<double>;
  struct Result {
    UnboundedMatrix<double> rotation;
    UnboundedVector<Scalar> eigenvalues;
  };
  static Rotation Identity(UnboundedMatrix<Scalar> const& m);
  static Result Uninitialized(UnboundedMatrix<Scalar> const& m);
};

template<typename Scalar_, int dimension>
struct ClassicalJacobiGenerator<FixedMatrix<Scalar_, dimension, dimension>> {
  using Scalar = Scalar_;
  using Rotation = FixedMatrix<double, dimension, dimension>;
  struct Result {
    FixedMatrix<double, dimension, dimension> rotation;
    FixedVector<Scalar, dimension> eigenvalues;
  };
  static Rotation Identity(FixedMatrix<Scalar, dimension, dimension> const& m);
  static Result Uninitialized(
      FixedMatrix<Scalar, dimension, dimension> const& m);
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
struct RayleighQuotientIterationGenerator<
    UnboundedMatrix<MScalar>,
    UnboundedVector<VScalar>> {
  struct Result {
    UnboundedVector<VScalar> eigenvector;
    MScalar eigenvalue;
  };
  static Result Uninitialized(UnboundedVector<VScalar> const& v);
};

template<typename MScalar, typename VScalar, int dimension>
struct RayleighQuotientIterationGenerator<
    FixedMatrix<MScalar, dimension, dimension>,
    FixedVector<VScalar, dimension>> {
  struct Result {
    FixedVector<VScalar, dimension> eigenvector;
    MScalar eigenvalue;
  };
  static Result Uninitialized(FixedVector<VScalar, dimension> const& v);
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

template<typename Scalar_>
auto ClassicalJacobiGenerator<UnboundedMatrix<Scalar_>>::Identity(
    UnboundedMatrix<Scalar_> const& m) -> Rotation {
  return UnboundedMatrix<Scalar>::Identity(m.rows(), m.columns());
}

template<typename Scalar_>
auto ClassicalJacobiGenerator<UnboundedMatrix<Scalar_>>::
Uninitialized(UnboundedMatrix<Scalar> const& m) -> Result {
  return {.rotation = UnboundedMatrix<double>(m.rows(), m.columns()),
          .eigenvalues = UnboundedVector<Scalar>(m.columns())};
}

template<typename Scalar_, int dimension>
auto ClassicalJacobiGenerator<FixedMatrix<Scalar_, dimension, dimension>>::
Identity(FixedMatrix<Scalar_, dimension, dimension> const& m) -> Rotation {
  return Rotation::Identity();
}

template<typename Scalar_, int dimension>
auto ClassicalJacobiGenerator<FixedMatrix<Scalar_, dimension, dimension>>::
Uninitialized(FixedMatrix<Scalar, dimension, dimension> const& m) -> Result {
  return {.rotation = FixedMatrix<double, dimension, dimension>(),
          .eigenvalues = FixedVector<Scalar, dimension>()};
}

template<typename MScalar, typename VScalar>
auto RayleighQuotientIterationGenerator<
    UnboundedMatrix<MScalar>,
    UnboundedVector<VScalar>>::
Uninitialized(UnboundedVector<VScalar> const& v) -> Result {
  return {UnboundedVector<VScalar>(v.size(), uninitialized), MScalar()};
}

template<typename MScalar, typename VScalar, int dimension>
auto RayleighQuotientIterationGenerator<
    FixedMatrix<MScalar, dimension, dimension>,
    FixedVector<VScalar, dimension>>::
Uninitialized(FixedVector<VScalar, dimension> const& v) -> Result {
  return {FixedVector<VScalar, dimension>(uninitialized), MScalar()};
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
        Σrₖᵢrₖⱼ += R(k, i) * R(k, j);
      }
      R(i, j) = (A(i, j) - Σrₖᵢrₖⱼ) / R(i, i);
    }
    Scalar Σrₖⱼ²{};
    for (int k = 0; k < j; ++k) {
      Σrₖⱼ² += Pow<2>(R(k, j));
    }
    // This will produce NaNs if the matrix is not positive definite.
    R(j, j) = Sqrt(A(j, j) - Σrₖⱼ²);
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
      Σrₖᵢ²dₖ += Pow<2>(R(k, i)) * D[k];
    }
    D[i] = A(i, i) - Σrₖᵢ²dₖ;
    for (int j = i + 1; j < A.columns(); ++j) {
      Scalar Σrₖᵢrₖⱼdₖ{};
      for (int k = 0; k < i; ++k) {
        Σrₖᵢrₖⱼdₖ += R(k, i) * R(k, j) * D[k];
      }
      R(i, j) = (A(i, j) - Σrₖᵢrₖⱼdₖ) / D[i];
    }
    R(i, i) = 1;
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
  x[n] = b[n] / U(n, n);
  for (int i = n - 1; i >= 0; --i) {
    auto s = b[i];
    for (int j = i + 1; j <= n; ++j) {
      s -= U(i, j) * x[j];
    }
    x[i] = s / U(i, i);
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
  x[0] = b[0] / L(0, 0);
  for (int i = 1; i < b.size(); ++i) {
    auto s = b[i];
    for (int j = 0; j < i; ++j) {
      s -= L(i, j) * x[j];
    }
    x[i] = s / L(i, i);
  }
  return x;
}

template<typename Matrix>
typename HessenbergDecompositionGenerator<Matrix>::Result
HessenbergDecomposition(Matrix const& A) {
  using G = HessenbergDecompositionGenerator<Matrix>;
  using Scalar = typename G::Scalar;
  typename HessenbergDecompositionGenerator<Matrix>::Result result{.H = A};
  auto& H = result.H;
  int const n = A.rows();

  // [GV13], Algorithm 7.4.2.
  for (int k = 0; k < n - 2; ++k) {
    auto const P = ComputeHouseholderReflection(
        ColumnView<Scalar, Matrix>{.matrix = H,
                                   .first_row = k + 1,
                                   .last_row = n - 1,
                                   .column = k});
    {
      auto block = BlockView<Scalar, Matrix>{.matrix = H,
                                             .first_row = k + 1,
                                             .last_row = n - 1,
                                             .first_column = k,
                                             .last_column = n - 1};
      Premultiply(P, block);
    }
    {
      auto block = BlockView<Scalar, Matrix>{.matrix = H,
                                             .first_row = 0,
                                             .last_row = n - 1,
                                             .first_column = k + 1,
                                             .last_column = n - 1};
      PostMultiply(block, P);
    }
  }
  return result;
}

template<typename Matrix>
typename QRDecompositionGenerator<Matrix>::Result
QRDecomposition(Matrix const& A, double const ε) {
  using G = QRDecompositionGenerator<Matrix>;
  using Scalar = typename G::Scalar;
  static const auto zero = Scalar{};
  auto hessenberg = HessenbergDecomposition(A);
  auto& H = hessenberg.H;
  int const n = H.rows();
  for (;;){
    for (int i = 1; i < n - 1; ++i) {
      if (H(i, i - 1) <= ε * (H(i, i) + H(i - 1, i - 1))) {
        H(i, i - 1) = zero;
      }
    }

    // Upper quasi-triangular means that we don't have consecutive nonzero
    // subdiagonal elements.
    bool has_subdiagonal_element = false;
    int q = 1;
    for (; q < n; ++q) {
      if (H(n - q, n - q - 1) == zero) {
        has_subdiagonal_element = false;
      } else {
        if (has_subdiagonal_element) {
          break;
        } else {
          has_subdiagonal_element = true;
        }
      }
    }

    if (q == n) {
      break;
    }

    int p = n - q;
    for (; p > 0; --p) {
      if (H(p, p - 1) == zero) {
        break;
      }
    }

    auto H₂₂ = BlockView<Scalar, Matrix>{.matrix = H,
                                         .first_row = p,
                                         .last_row = n - q - 1,
                                         .first_column = p,
                                         .last_column = n - q - 1};
    FrancisQRStep<Scalar>(H₂₂);
  }
  //TODO(phl)Diagonal blocks?
  typename QRDecompositionGenerator<Matrix>::Result result{.R = H};
  return result;
}

template<typename Matrix>
typename ClassicalJacobiGenerator<Matrix>::Result
ClassicalJacobi(Matrix const& A,  int max_iterations, double const ε) {
  using G = ClassicalJacobiGenerator<Matrix>;
  using Scalar = typename G::Scalar;
  auto result = G::Uninitialized(A);
  auto& V = result.rotation;

  // [GV13], Algorithm 8.5.2.
  auto const identity = G::Identity(A);
  auto const A_frobenius_norm = A.FrobeniusNorm();
  V = identity;
  auto diagonalized_A = A;
  for (int k = 0; k < max_iterations; ++k) {
    Scalar max_Apq{};
    int max_p = -1;
    int max_q = -1;

    // Find the largest off-diagonal element and exit if it's small.
    for (int p = 0; p < diagonalized_A.rows(); ++p) {
      for (int q = p + 1; q < diagonalized_A.columns(); ++q) {
        Scalar const abs_Apq = Abs(diagonalized_A(p, q));
        if (abs_Apq >= max_Apq) {
          max_Apq = abs_Apq;
          max_p = p;
          max_q = q;
        }
      }
    }
    if (max_Apq <= ε * A_frobenius_norm) {
      break;
    }

    auto const J =
        SymmetricSchurDecomposition2By2<Scalar>(diagonalized_A, max_p, max_q);

    // A = ᵗJ A J
    PostMultiply(diagonalized_A, J);
    PremultiplyByTranspose(J, diagonalized_A);

    // V = V J
    PostMultiply(V, J);
    if (k == max_iterations - 1) {
      LOG(ERROR) << "Difficult diagonalization: " << A
                 << ", stopping with: " << diagonalized_A;
    }
  }

  for (int i = 0; i < A.rows(); ++i) {
    result.eigenvalues[i] = diagonalized_A(i, i);
  }
  return result;
}

template<typename Matrix, typename Vector>
typename RayleighQuotientGenerator<Matrix, Vector>::Result
RayleighQuotient(Matrix const& A, Vector const& x) {
  // [GV13], section 8.2.3.
  return x.Transpose() * (A * x) / (x.Transpose() * x);
}

template<typename Matrix, typename Vector>
typename RayleighQuotientIterationGenerator<Matrix, Vector>::Result
RayleighQuotientIteration(Matrix const& A, Vector const& x) {
  using G = RayleighQuotientIterationGenerator<Matrix, Vector>;
  auto result = G::Uninitialized(x);
  auto& xₖ = result.eigenvector;
  auto& μₖ = result.eigenvalue;

  // [GV13], section 8.2.3.
  xₖ = x / x.Norm();
  for (int iteration = 0; iteration < 10; ++iteration) {
    μₖ = RayleighQuotient(A, xₖ);
    auto A_minus_μₖ_I = A;
    for (int i = 0; i < A.rows(); ++i) {
      A_minus_μₖ_I(i, i) -= μₖ;
    }
    auto const residual = (A_minus_μₖ_I * xₖ).Norm();
    // TODO(phl): This test is cheezy.  It should be based on the norm of the
    // matrix.
    if (residual < 2 * std::numeric_limits<double>::epsilon() *
                       si::Unit<decltype(residual)>) {
      return result;
    }
    auto const zₖ₊₁ = Solve(A_minus_μₖ_I, xₖ);
    xₖ = zₖ₊₁ / zₖ₊₁.Norm();
  }
  LOG(WARNING) << "Unconverged Rayleigh quotient iteration: " << A;
  return result;
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
      if (Abs(A(i, k)) >= max) {
        r = i;
        max = Abs(A(i, k));
      }
    }
    CHECK_LE(0, r) << A << " cannot pivot";
    CHECK_LT(r, A.rows()) << A << " cannot pivot";

    // Swap the rows of A.
    for (int i = 0; i < A.columns(); ++i) {
      std::swap(A(k, i), A(r, i));
    }
    // Swap the rows of L.
    for (int i = 0; i < k; ++i) {
      std::swap(L(k, i), L(r, i));
    }
    // Swap the rows of b.
    std::swap(b[k], b[r]);

    LOG_IF(WARNING, A(k, k) == Scalar{})
        << A << " does not have a unique LU decomposition";

    for (int j = k; j < A.columns(); ++j) {
      auto U_kj = A(k, j);
      for (int i = 0; i < k; ++i) {
        U_kj -= L(k, i) * U(i, j);
      }
      U(k, j) = U_kj;
    }
    for (int i = k + 1; i < A.rows(); ++i) {
      auto L_ik = A(i, k);
      for (int j = 0; j < k; ++j) {
        L_ik -= L(i, j) * U(j, k);
      }
      L(i, k) = L_ik / U(k, k);
    }
    L(k, k) = 1;
  }

  // For the resolution of triangular systems see [Hig02], Algorithm 8.1 p. 140.

  // Find y such that L * y = P * b.
  auto const y = ForwardSubstitution(L, b);
  // Find x such that U * x = y.
  auto const x = BackSubstitution(U, y);

  return x;
}

}  // namespace internal
}  // namespace _matrix_computations
}  // namespace numerics
}  // namespace principia
