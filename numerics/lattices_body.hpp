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
struct LenstraLenstraLovászGenerator<UnboundedMatrix<Scalar>> {
  using Vector = UnboundedVector<Scalar>;
};

template<typename Scalar, int rows, int columns>
struct LenstraLenstraLovászGenerator<
    FixedMatrix<Scalar, rows, columns>> {
  using Vector = FixedVector<Scalar, rows>;
};

// This implements [HPS], theorem 7.71, figure 7.8.  Note that figures 7.9 and
// 7.10 are supposedly more efficient, but they are significantly more
// complicated.  If performance is an issue, we should look into recent
// improvements of LLL.
template<typename Matrix>
  requires two_dimensional<Matrix>
Matrix LenstraLenstraLovász(Matrix const& L) {
  using G = LenstraLenstraLovászGenerator<Matrix>;
  auto const n = L.columns();
  auto const m = L.rows();
  auto v = L;
  for (int k = 1; k < n;) {
    auto qr = UnitriangularGramSchmidt(v);
    auto vₖ = ColumnView{.matrix = v,
                        .first_row = 0,
                        .last_row = m - 1,
                        .column = k};
    for (int j = k - 1; j >= 0; --j) {
      auto const μₖⱼ = qr.R(j, k);
      auto vⱼ = ColumnView{.matrix = v,
                           .first_row = 0,
                           .last_row = m - 1,
                           .column = j};
      auto const round_μₖⱼ = std::round(μₖⱼ);
      if (round_μₖⱼ != 0) {
        vₖ -= round_μₖⱼ * typename G::Vector(vⱼ);
        qr = UnitriangularGramSchmidt(v);
      }
    }
    auto const μₖₖ₋₁ = qr.R(k - 1, k);
    auto v𐌟ₖ = ColumnView{.matrix = qr.Q,
                         .first_row = 0,
                         .last_row = m - 1,
                         .column = k};
    auto v𐌟ₖ₋₁ = ColumnView{.matrix = qr.Q,
                           .first_row = 0,
                           .last_row = m - 1,
                           .column = k - 1};
    if (v𐌟ₖ.Norm²() >= (0.75 - Pow<2>(μₖₖ₋₁)) * v𐌟ₖ₋₁.Norm²()) {
      ++k;
    } else {
      auto vₖ₋₁ = ColumnView{.matrix = v,
                            .first_row = 0,
                            .last_row = m - 1,
                            .column = k - 1};
      SwapColumns(vₖ₋₁, vₖ);
      k = std::max(k - 1, 1);
    }
  }
  return v;
}

}  // namespace internal
}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
