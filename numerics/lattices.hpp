#pragma once

#include "numerics/concepts.hpp"

namespace principia {
namespace numerics {
namespace _lattices {
namespace internal {

using namespace principia::numerics::_concepts;

// In this function, the base of the lattice is in *columns*.  This differs from
// usual conventions, and is due to the fact that our QR factorization returns
// orthogonal *columns*.
template<typename Matrix>
  requires two_dimensional<Matrix>
Matrix LenstraLenstraLovász(Matrix const& L);

// Same convention as above, but using a more efficient algorithm.
template<typename Matrix>
  requires two_dimensional<Matrix>
Matrix NguyễnStehlé(Matrix const& L);

}  // namespace internal

using internal::LenstraLenstraLovász;
using internal::NguyễnStehlé;

}  // namespace _lattices
}  // namespace numerics
}  // namespace principia

#include "numerics/lattices_body.hpp"
