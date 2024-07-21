#pragma once

#include "numerics/lattices.hpp"

#include <algorithm>

#include "boost/multiprecision/cpp_int.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _lattices {
namespace internal {

using namespace boost::multiprecision;
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

template<typename Matrix>
struct GramGenerator;

template<typename Scalar>
struct GramGenerator<UnboundedMatrix<Scalar>> {
  using Result = UnboundedMatrix<Scalar>;
  static Result Uninitialized(UnboundedMatrix<Scalar> const& m);
};

template<typename Scalar, int rows, int columns>
struct GramGenerator<FixedMatrix<Scalar, rows, columns>> {
  using Result = Fixed<Scalar, columns, columns>;
  static Result Uninitialized(FixedMatrix<Scalar, rows, columns> const& m);
};

template<typename Matrix>
struct NguyễnStehléGenerator;

template<typename Scalar>
struct NguyễnStehléGenerator<UnboundedMatrix<Scalar>> {
  using R = UnboundedMatrix<Scalar>;
};

template<>
struct NguyễnStehléGenerator<UnboundedMatrix<cpp_int>> {
  using R = UnboundedMatrix<double>;
};

template<typename Scalar, int rows, int columns>
struct NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>> {
  using R = FixedMatrix<Scalar, rows, columns>;
};

template<int rows, int columns>
struct NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>> {
  using R = FixedMatrix<cpp_int, rows, columns>;
};


template<typename Scalar>
auto GramGenerator<UnboundedMatrix<Scalar>>::Uninitialized(
UnboundedMatrix<Scalar> const& m) -> Result {
  return Result(m.rows(), m.rows(), uninitialized);
}

template<typename Scalar, int rows, int columns>
auto GramGenerator<FixedMatrix<Scalar, rows, columns>>::Uninitialized(
FixedMatrix<Scalar, rows, columns> const& m) -> Result {
  return Result(uninitialized);
}


template<typename Matrix>
typename GramGenerator<Matrix>::Result Gram(Matrix const& L) {
  using G = GramGenerator<Matrix>;
  auto result = G::Uninitialized(L);
  for (std::int64_t i = 0; i < columns) {
    auto const bᵢ = ColumnView{.matrix = L,
                               .first_row = 0,
                               .last_row = rows - 1,
                               .column = i};
    for (std::int64_t j = 0; j <= i) {
      auto const bⱼ = ColumnView{.matrix = L,
                                 .first_row = 0,
                                 .last_row = rows - 1,
                                 .column = j};
      auto const bᵢbⱼ = bᵢ * bⱼ;
      result(i, j) = bᵢbⱼ;
      result(j, i) = bᵢbⱼ;
    }
  }
  return result;
}

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
      auto const round_μₖⱼ = Round(μₖⱼ);
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

template<typename Matrix>
  requires two_dimensional<Matrix>
Matrix NguyễnStehlé(Matrix const& L) {
  using G = NguyễnStehléGenerator<Matrix>;
  //[NS09] figure 7.
  double const ẟ = 0.75;
  double const η = 0.55;
  //[NS09] Section 6.
  auto b = L;
  // Step 1.
  auto const G = Gram(b);
  // Step 2.
  double const ẟ = (ẟ + 1) / 2;
  auto const b₀ = ColumnView{.matrix = b,
                             .first_row = 0,
                             .last_row = b.rows(),
                             .column = 0};
  typename G::R r = G::Uninitialized();
  r(0, 0) = static_cast<typename G::R::ElementType>(b₀.Norm²());
  std::int64_t κ = 1;
  std::int64_t ζ = -1;
  while (κ < b.columns()) {
    // Step 3.
    SizeReduce(b, ζ + 1, κ);
    // Step 4.
    std::int64_t κʹ = κ;
    while (κ >= ζ + 2 &&  ẟ * r(κ - 1, κ - 1) >= s(κ - 1)) {
      --κ;
    }
    // Step 5.
    for (std::int64_t i = ζ + 1; i < κ; ++i) {
      μ(κ, i) = μ(κʹ, i);
      r(κ, i) = r(κʹ, i);
    }
    r(κ, κ) = s(κ);
    // Step 6.
    Insert(b, κʹ, κ);
    UpdateGram(G, b);
    // Step 7.
    auto const bκ = ColumnView{.matrix = b,
                               .first_row = 0,
                               .last_row = b.rows(),
                               .column = κ};
    if (bκ == zero) {
      ++ζ;
    }
    // Step 8.
    κ = std::max(ζ + 2, κ + 1);
  }
}

}  // namespace internal
}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
