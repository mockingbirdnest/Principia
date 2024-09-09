#pragma once

#include <limits>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _matrix_computations {
namespace internal {

using namespace principia::quantities::_named_quantities;

// Declares:
//   using Result = ⟨upper triangular matrix⟩;
template<typename U>
struct CholeskyDecompositionGenerator;

// Declares:
//   struct Result {
//     ⟨upper triangular matrix⟩ R;
//     ⟨vector⟩ D;
//   };
template<typename V, typename U>
struct ᵗRDRDecompositionGenerator;

// Declares:
//   using Result = ⟨vector⟩;
template<typename M, typename V>
struct SubstitutionGenerator;

// Declares:
//   struct Result {
//     ⟨upper triangular matrix⟩ R;
//     ⟨matrix⟩ Q;
//   };
template<typename M>
struct GramSchmidtGenerator;

// Declares:
//   struct Result {
//     ⟨upper triangular matrix⟩ R;
//     ⟨matrix⟩ Q;
//   };
template<typename M>
struct UnitriangularGramSchmidtGenerator;

// Declares:
//   struct Result {
//     (matrix) H;
//     (matrix) U;
//   }
// TODO(phl): Add support for U.
template<typename M>
struct HessenbergDecompositionGenerator;

// Declares:
//   struct Result {
//     (matrix) T;
//     (matrix) Q;
//   };
// TODO(phl): Add support for Q.
template<typename M>
struct RealSchurDecompositionGenerator;

// Declares:
//   struct Result {
//     ⟨matrix⟩ rotation;
//     ⟨vector⟩ eigenvalues;
//   };
// Note that in `rotation` the eigenvectors appear in column.  They are
// normalized.
template<typename M>
struct ClassicalJacobiGenerator;

// Declares:
//   using Result = ⟨scalar⟩;
template<typename M, typename V>
struct RayleighQuotientGenerator;

// Declares:
//   struct Result {
//     ⟨vector⟩ eigenvector;
//     ⟨scalar⟩ eigenvalue;
//   };
template<typename M, typename V>
struct RayleighQuotientIterationGenerator;

// Declares:
//   using Result = ⟨vector⟩;
template<typename M, typename V>
struct SolveGenerator;

// If A is the upper half of a symmetric positive definite matrix, returns R so
// that A = ᵗR R.
template<typename UpperTriangularMatrix>
typename CholeskyDecompositionGenerator<UpperTriangularMatrix>::Result
CholeskyDecomposition(UpperTriangularMatrix const& A);

// If A is the upper half of a symmetric matrix, returns R and D so that
// A = ᵗR D R.  The diagonal matrix is represented as a vector.
template<typename Vector, typename UpperTriangularMatrix>
typename ᵗRDRDecompositionGenerator<Vector, UpperTriangularMatrix>::Result
ᵗRDRDecomposition(UpperTriangularMatrix const& A);

// Returns x such that U x = b.
template<typename UpperTriangularMatrix, typename Vector>
typename SubstitutionGenerator<UpperTriangularMatrix, Vector>::Result
BackSubstitution(UpperTriangularMatrix const& U,
                 Vector const& b);

// Returns x such that L x = b.
template<typename LowerTriangularMatrix, typename Vector>
typename SubstitutionGenerator<LowerTriangularMatrix, Vector>::Result
ForwardSubstitution(LowerTriangularMatrix const& L,
                    Vector const& b);

// Returns Q and R such that A = Q R where R is upper triangular and Q
// orthogonal.
template<typename Matrix>
typename GramSchmidtGenerator<Matrix>::Result
ClassicalGramSchmidt(Matrix const& A);

// Returns Q and R such that A = Q R where R is upper unitriangular and Q
// has orthogonal columns.
template<typename Matrix>
typename UnitriangularGramSchmidtGenerator<Matrix>::Result
UnitriangularGramSchmidt(Matrix const& A);

// If A is a square matrix, returns U and H so that A = ᵗU H U, where H is an
// upper Hessenberg matrix.
// TODO(phl): Add support for returning U.
template<typename Matrix>
typename HessenbergDecompositionGenerator<Matrix>::Result
HessenbergDecomposition(Matrix const& A);

// If A is a square matrix, returns Q and T so that A = Q T ᵗQ, where T is upper
// quasi-triangular.
template<typename Matrix>
typename RealSchurDecompositionGenerator<Matrix>::Result
RealSchurDecomposition(Matrix const& A, double ε);

// Returns the eigensystem of A, which must be symmetric.
// As a safety measure we limit the number of iterations.  We prefer to exit
// when the matrix is nearly diagonal, though.
template<typename Matrix>
typename ClassicalJacobiGenerator<Matrix>::Result ClassicalJacobi(
    Matrix const& A,
    std::int64_t max_iterations = 16,
    double ε = std::numeric_limits<double>::epsilon() / 128);

// Returns the Rayleigh quotient r(x) = ᵗx A x / ᵗx x.
template<typename Matrix, typename Vector>
typename RayleighQuotientGenerator<Matrix, Vector>::Result
RayleighQuotient(Matrix const& A, Vector const& x);

// Returns the eigenvector closest to x and its eigenvalue.  A must be
// symmetric.
template<typename Matrix, typename Vector>
typename RayleighQuotientIterationGenerator<Matrix, Vector>::Result
RayleighQuotientIteration(Matrix const& A, Vector const& x);

// Returns x such that A x = b.
template<typename Matrix, typename Vector>
typename SolveGenerator<Matrix, Vector>::Result
Solve(Matrix A, Vector b);

}  // namespace internal

using internal::BackSubstitution;
using internal::CholeskyDecomposition;
using internal::ClassicalGramSchmidt;
using internal::ClassicalJacobi;
using internal::ForwardSubstitution;
using internal::HessenbergDecomposition;
using internal::RayleighQuotient;
using internal::RayleighQuotientIteration;
using internal::RealSchurDecomposition;
using internal::Solve;
using internal::UnitriangularGramSchmidt;
using internal::ᵗRDRDecomposition;

}  // namespace _matrix_computations
}  // namespace numerics
}  // namespace principia

#include "numerics/matrix_computations_body.hpp"
