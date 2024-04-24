#pragma once

#include "numerics/lattices.hpp"

#include "base/tags.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/transposed_view.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _lattices {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

template<typename Matrix>
struct LenstraLenstraLovászGenerator;

template<typename Scalar>
struct LenstraLenstraLovászGenerator<UnboundedMatrix<Scalar>> {
  using Vector = UnboundedVector<Scalar>;
  using DoubleMatrix = UnboundedMatrix<double>;
  using Norm²Vector = UnboundedVector<Square<Scalar>>;
};

template<typename Scalar, int rows, int columns>
struct LenstraLenstraLovászGenerator<
    FixedMatrix<Scalar, rows, columns>> {
  using Vector = FixedVector<Scalar, rows>;
  using DoubleMatrix = FixedMatrix<double, rows, columns>;
  using Norm²Vector = FixedVector<Square<Scalar>, rows>;
};

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
      vₖ -= std::round(μₖⱼ) * typename G::Vector(vⱼ);
      qr = UnitriangularGramSchmidt(v);
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

template<typename Matrix>
  requires two_dimensional<Matrix>
Matrix LenstraLenstraLovász2(Matrix const& L) {
  using G = LenstraLenstraLovászGenerator<Matrix>;
  auto const n = L.columns();
  auto const rows = L.rows();
  Matrix v = L;
  Matrix v𐌟;
  typename G::DoubleMatrix μ(uninitialized);
  typename G::Norm²Vector B(uninitialized);
  int k_max = 0;

  auto Red = [&](int const k, int const l) {
    if (Abs(μ(k, l)) <= 0.5) {
      return;
    }
    auto const m = std::round(μ(k, l));
    auto vₖ = ColumnView{.matrix = v,
                        .first_row = 0,
                        .last_row = rows - 1,
                        .column = k};
    auto vₗ = ColumnView{.matrix = v,
                        .first_row = 0,
                        .last_row = rows - 1,
                        .column = l};
    vₖ -= m * typename G::Vector(vₗ);
    μ(k, l) -= m;
    for (int i = 0; i < l - 1; ++i) {
      μ(k, i) -= m * μ(l, i);
    }
  };

  auto Swap = [&](int const k) {
    auto vₖ = ColumnView{.matrix = v,
                        .first_row = 0,
                        .last_row = rows - 1,
                        .column = k};
    auto vₖ₋₁ = ColumnView{.matrix = v,
                          .first_row = 0,
                          .last_row = rows - 1,
                          .column = k - 1};
    SwapColumns(vₖ, vₖ₋₁);
    for (int j = 0; j < k - 2; ++j) {
      std::swap(μ(k - 1, j), μ(k, j));
    }
    auto const μ_value = μ(k, k - 1);
    auto const B_value = B[k] + Pow<2>(μ_value) * B[k - 1];
    μ(k, k - 1) = μ_value * B[k - 1] / B_value;
    B[k] = B[k - 1] * B[k] / B_value;
    B[k - 1] = B_value;
    for (int i = k + 1; i < k_max; ++i) {
      auto const m = μ(i, k);
      μ(i, k) = μ(i, k - 1) - μ_value * m;
      μ(i, k - 1) = m + μ(k, k - 1) * μ(i, k);
    }
  };

  ColumnView<Matrix> v𐌟₀ = ColumnView{.matrix = v𐌟,
                        .first_row = 0,
                        .last_row = rows - 1,
                        .column = 0};
  ColumnView<Matrix> v₀ = ColumnView{.matrix = v,
                       .first_row = 0,
                       .last_row = rows - 1,
                       .column = 0};
  v𐌟₀ = v₀;
  B[0] = v₀.Norm²();
  for (int k = 1; k < n;) {
    if (k > k_max) {
      k_max = k;
      auto v𐌟ₖ = ColumnView{.matrix = v𐌟,
                           .first_row = 0,
                           .last_row = rows - 1,
                           .column = k};
      auto vₖ = ColumnView{.matrix = v,
                           .first_row = 0,
                           .last_row = rows - 1,
                           .column = k};
      v𐌟ₖ = vₖ;
      for (int j = 0; j < k - 1; ++j) {
        auto v𐌟ⱼ = ColumnView{.matrix = v𐌟,
                              .first_row = 0,
                              .last_row = rows - 1,
                              .column = j};
        μ(k, j) = TransposedView{vₖ} * v𐌟ⱼ / B[j];
        v𐌟ₖ -= μ(k, j) * typename G::Vector(v𐌟ⱼ);
      }
      B[k] = v𐌟ₖ.Norm²();
    }
    for (;;) {
      Red(k, k - 1);
      if (B[k] < (0.75 - Pow<2>(μ(k, k - 1))) * B[k - 1]) {
        Swap(k);
        k = std::max(k - 1, 1);
      } else {
        for (int l = k - 3; l >= 0; --l) {
          Red(k, l);
        }
        ++k;
        break;
      }
    }
  }
  return v;
}

}  // namespace internal
}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
