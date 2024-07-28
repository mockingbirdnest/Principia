#pragma once

#include "numerics/lattices.hpp"

#include <algorithm>

#include "base/tags.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _lattices {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::base::_tags;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

// In the terminology of [NS09], our vectors are in columns, so |d| is |columns|
// and |n| is |rows|.
template<typename Matrix>
struct GramGenerator;

template<typename Scalar>
struct GramGenerator<UnboundedMatrix<Scalar>> {
  using Result = UnboundedMatrix<Square<Scalar>>;
  static Result Uninitialized(UnboundedMatrix<Scalar> const& m);
};

template<typename Scalar, int rows, int columns>
struct GramGenerator<FixedMatrix<Scalar, rows, columns>> {
  using Result = FixedMatrix<Square<Scalar>, columns, columns>;
  static Result Uninitialized(FixedMatrix<Scalar, rows, columns> const& m);
};

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
struct NguyễnStehléGenerator;

template<typename Scalar>
struct NguyễnStehléGenerator<UnboundedMatrix<Scalar>> {
  using R = UnboundedMatrix<Scalar>;
  static R Uninitialized(UnboundedMatrix<Scalar> const& m);
};

template<>
struct NguyễnStehléGenerator<UnboundedMatrix<cpp_int>> {
  using R = UnboundedMatrix<double>;
  static R Uninitialized(UnboundedMatrix<cpp_int> const& m);
};

template<typename Scalar, int rows, int columns>
struct NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>> {
  using R = FixedMatrix<Scalar, rows, columns>;
  static R Uninitialized(FixedMatrix<Scalar, rows, columns> const& m);
};

template<int rows, int columns>
struct NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>> {
  using R = FixedMatrix<cpp_int, rows, columns>;
  static R Uninitialized(FixedMatrix<cpp_int, rows, columns> const& m);
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

template<typename Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::Uninitialized(
    UnboundedMatrix<Scalar> const& m) -> R {
  return R(m.rows(), m.columns(), uninitialized);
}

auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::Uninitialized(
    UnboundedMatrix<cpp_int> const& m) -> R {
  return R(m.rows(), m.columns(), uninitialized);
}

template<typename Scalar, int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::Uninitialized(
    FixedMatrix<Scalar, rows, columns> const& m) -> R {
  return R(m.rows(), m.columns(), uninitialized);
}

template<int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>>::Uninitialized(
    FixedMatrix<cpp_int, rows, columns> const& m) -> R {
  return R(m.rows(), m.columns(), uninitialized);
}


template<typename Matrix>
typename GramGenerator<Matrix>::Result Gram(Matrix const& L) {
  using G = GramGenerator<Matrix>;
  std::int64_t const rows = L.rows();
  std::int64_t const columns = L.columns();
  auto result = G::Uninitialized(L);
  for (std::int64_t i = 0; i < columns; ++i) {
    auto const bᵢ = ColumnView{.matrix = L,
                               .first_row = 0,
                               .last_row = rows - 1,
                               .column = i};
    for (std::int64_t j = 0; j <= i; ++j) {
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
  using Gen = NguyễnStehléGenerator<Matrix>;
  auto b = L;
  std::int64_t const d = b.columns();

  //[NS09] figure 7.
  double const ẟ = 0.75;
  double const η = 0.55;
  //[NS09] figure 9.
  // Step 1.
  auto const G = Gram(b);
  // Step 2.
  double const ẟ = (ẟ + 1) / 2;
  auto const b₀ = ColumnView{.matrix = b,
                             .first_row = 0,
                             .last_row = b.rows(),
                             .column = 0};
  typename Gen::R r = Gen::Uninitialized(b);
  r(0, 0) = static_cast<typename Gen::R::ElementType>(b₀.Norm²());
  std::int64_t κ = 1;
  std::int64_t ζ = -1;
  while (κ < d) {
    // Step 3.
    SizeReduce(b, ζ + 1, κ);
    // Step 4.
    //TODO(phl)high index probably useless
    std::int64_t κʹ = κ;
    while (κ >= ζ + 2 &&  ẟ * r(κ - 1, κ - 1) >= s(κ - 1)) {
      --κ;
    }
    // Step 5.
    for (std::int64_t i = ζ + 1; i < κ - 1; ++i) {
      μ(κ, i) = μ(κʹ, i);
      r(κ, i) = r(κʹ, i);
    }
    r(κ, κ) = s(κ);
    // Step 6.
    Insert(b, κʹ, κ, G);
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

  return b;
}

}  // namespace internal
}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
