
#pragma once

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_arrays {

using quantities::Quotient;
using quantities::SquareRoot;

template<typename T>
struct CholeskyGenerator;

template<typename Scalar, template<typename S> typename UpperTriangularMatrix>
struct CholeskyGenerator<UpperTriangularMatrix<Scalar>> {
  using type = UpperTriangularMatrix<SquareRoot<Scalar>>;
};

template<typename Scalar, int columns,
         template<typename S, int c> typename UpperTriangularMatrix>
struct CholeskyGenerator<UpperTriangularMatrix<Scalar, columns>> {
  using type = UpperTriangularMatrix<SquareRoot<Scalar>, columns>;
};

template<typename V, typename U>
struct ᵗRDRGenerator;

template<typename Scalar,
         template<typename S> typename Vector,
         template<typename S> typename UpperTriangularMatrix>
struct ᵗRDRGenerator<Vector<Scalar>, UpperTriangularMatrix<Scalar>> {
  struct type {
    UpperTriangularMatrix<double> R;
    Vector<Scalar> D;
  };
};

template<typename Scalar, int columns,
         template<typename S, int c> typename Vector,
         template<typename S, int c> typename UpperTriangularMatrix>
struct ᵗRDRGenerator<Vector<Scalar, columns>,
                     UpperTriangularMatrix<Scalar, columns>> {
  struct type {
    UpperTriangularMatrix<double, columns> R;
    Vector<Scalar, columns>& D;
  };
};

template<typename M, typename V>
struct SubstitutionGenerator;

template<typename LScalar, typename RScalar,
         template<typename S> typename Matrix,
         template<typename S> typename Vector>
struct SubstitutionGenerator<Matrix<LScalar>, Vector<RScalar>> {
  using type = Vector<Quotient<RScalar, LScalar>>;
};

template<typename LScalar, typename RScalar, int dimension,
         template<typename S, int d> typename Matrix,
         template<typename S, int d> typename Vector>
struct SubstitutionGenerator<Matrix<LScalar, dimension>,
                             Vector<RScalar, dimension>> {
  using type = Vector<Quotient<RScalar, LScalar>, dimension>;
};


// If A is the upper half of a symmetric positive definite matrix, returns R so
// that A = ᵗR R.
template<typename UpperTriangularMatrix>
typename CholeskyGenerator<UpperTriangularMatrix>::type
CholeskyDecomposition(UpperTriangularMatrix const& A);

// If A is the upper half of a symmetric matrix, returns R and D so that
// A = ᵗR D R.  The diagonal matrix is represented as a vector.
template<typename Vector, typename UpperTriangularMatrix>
typename ᵗRDRGenerator<Vector, UpperTriangularMatrix>::type
ᵗRDRDecomposition(UpperTriangularMatrix const& A);

// Returns x such that U x = b.
template<typename UpperTriangularMatrix, typename Vector>
typename SubstitutionGenerator<UpperTriangularMatrix, Vector>::type
BackSubstitution(UpperTriangularMatrix const& U,
                 Vector const& b);

// Return x such that L x = b.
template<typename LowerTriangularMatrix, typename Vector>
typename SubstitutionGenerator<LowerTriangularMatrix, Vector>::type
ForwardSubstitution(LowerTriangularMatrix const& L,
                    Vector const& b);

}  // namespace internal_arrays

using internal_arrays::BackSubstitution;
using internal_arrays::CholeskyDecomposition;
using internal_arrays::ForwardSubstitution;
using internal_arrays::ᵗRDRDecomposition;

}  // namespace numerics
}  // namespace principia

#include "numerics/matrix_computations_body.hpp"
