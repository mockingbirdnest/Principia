#pragma once

#include "numerics/lattices.hpp"

#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _lattices {
namespace internal {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;

template<typename Matrix>
struct LenstraLenstraLovászGenerator;

template<typename Scalar>
struct LenstraLenstraLovászGenerator<UnboundedMatrix<Scalar>> {};

template<typename Scalar, int dimension>
struct LenstraLenstraLovászGenerator<
    FixedMatrix<Scalar, dimension, dimension>> {};

template<typename Matrix>
  requires two_dimensional<Matrix> &&
           std::is_integral_v<typename Matrix::Scalar>
Matrix LenstraLenstraLovász(Matrix const& L) {
  //TODO(phl):Rows/columns confusion.
  using G = LenstraLenstraLovászGenerator<Matrix>;
  auto const n = L.columns();
  auto const m = L.rows();
  auto v = L;
  for (int k = 1; k < n;) {
    auto const qr = UnitriangularGramSchmidt(L);
    auto vₖ = ColumnView{.matrix = v,
                        .first_row = 0,
                        .last_row = m,
                        .column = k};
    for (int j = k - 1; j >= 0; --j) {
      auto const μₖⱼ = qr.R(j, k);
      auto vⱼ = ColumnView{.matrix = v,
                          .first_row = 0,
                          .last_row = m,
                          .column = j};
      vₖ -= std::round(μₖⱼ) * typename G::Lvector(vⱼ);
    }
    auto const μₖₖ₋₁ = qr.R(k - 1, k);
    auto v𐌟ₖ = ColumnView{matrix = qr.Q,
                         .first_row = 0,
                         .last_row = m,
                         .column = k};
    auto v𐌟ₖ₋₁ = ColumnView{matrix = qr.Q,
                         .first_row = 0,
                         .last_row = m,
                         .column = k - 1};
    if (v𐌟ₖ.Norm²() >= (0.75 - Pow<2>(μₖₖ₋₁)) * v𐌟ₖ₋₁.Norm²()) {
      ++k;
    } else {
      auto vₖ₋₁ = ColumnView{.matrix = v,
                            .first_row = 0,
                            .last_row = m,
                            .column = k - 1};
      std::swap(vₖ₋₁, vₖ);
      k = std::max(k - 1, 1);
    }
  }
  return v;
}

}  // namespace internal
}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
