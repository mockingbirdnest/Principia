
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


// If A is the upper half of a symmetric positive definite matrix, returns R so
// that A = ᵗR R.
template<typename UpperTriangularMatrix>
typename CholeskyGenerator<UpperTriangularMatrix>::type
CholeskyDecomposition(UpperTriangularMatrix const& A);

// If A is the upper half of a symmetric matrix, returns R and D so that
// A = ᵗR D R.  The diagonal matrix is represented as a vector.
template<typename Scalar,
         template<typename S> typename UpperTriangularMatrix,
         template<typename S> typename Vector>
void ᵗRDRDecomposition(UpperTriangularMatrix<Scalar> const& A,
                       UpperTriangularMatrix<double>& R,
                       Vector<Scalar>& D);

// Returns x such that U x = b.
template<typename LScalar, typename RScalar,
         template<typename S> typename UpperTriangularMatrix,
         template<typename S> typename Vector>
Vector<Quotient<RScalar, LScalar>> BackSubstitution(
    UpperTriangularMatrix<LScalar> const& U,
    Vector<RScalar> const& b);

// Return x such that L x = b.
template<typename LScalar, typename RScalar,
         template<typename S> typename LowerTriangularMatrix,
         template<typename S> typename Vector>
Vector<Quotient<RScalar, LScalar>> ForwardSubstitution(
    LowerTriangularMatrix<LScalar> const& L,
    Vector<RScalar> const& b);

}  // namespace internal_arrays

using internal_arrays::BackSubstitution;
using internal_arrays::CholeskyDecomposition;
using internal_arrays::ForwardSubstitution;
using internal_arrays::ᵗRDRDecomposition;

}  // namespace numerics
}  // namespace principia

#include "numerics/matrix_computations_body.hpp"
