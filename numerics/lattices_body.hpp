#pragma once

#include "numerics/lattices.hpp"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "base/tags.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/transposed_view.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/concepts.hpp"
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
using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

// The largest `int64_t` that can be converted to and from `double` without loss
// of accuracy.
constexpr std::int64_t largest_convertible_integer =
    2LL << std::numeric_limits<double>::digits;

// In the terminology of [NS09], our vectors are in columns, so `d` is `columns`
// and `n` is `rows`.
template<typename Matrix>
struct GramGenerator;

template<typename Scalar>
struct GramGenerator<UnboundedMatrix<Scalar>> {
  using G = UnboundedMatrix<Square<Scalar>>;
  static G Uninitialized(UnboundedMatrix<Scalar> const& m);
};

template<typename Scalar, int rows, int columns>
struct GramGenerator<FixedMatrix<Scalar, rows, columns>> {
  using G = FixedMatrix<Square<Scalar>, columns, columns>;
  static G Uninitialized(FixedMatrix<Scalar, rows, columns> const& m);
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
  requires quantity<Scalar>
struct NguyễnStehléGenerator<UnboundedMatrix<Scalar>> {
  using R = UnboundedUpperTriangularMatrix<Square<Scalar>>;
  using Μ = UnboundedStrictlyUpperTriangularMatrix<double>;
  using S = UnboundedVector<Square<Scalar>>;
  using Vector = UnboundedVector<Scalar>;
  static R UninitializedR(UnboundedMatrix<Scalar> const& m);
  static Μ UninitializedΜ(UnboundedMatrix<Scalar> const& m);
  static S UninitializedS(UnboundedMatrix<Scalar> const& m);
  static Vector Zero(UnboundedMatrix<Scalar> const& m);
};

template<>
struct NguyễnStehléGenerator<UnboundedMatrix<cpp_int>> {
  using R = UnboundedUpperTriangularMatrix<double>;
  using Μ = UnboundedStrictlyUpperTriangularMatrix<double>;
  using S = UnboundedVector<double>;
  using Vector = UnboundedVector<cpp_int>;
  static R UninitializedR(UnboundedMatrix<cpp_int> const& m);
  static Μ UninitializedΜ(UnboundedMatrix<cpp_int> const& m);
  static S UninitializedS(UnboundedMatrix<cpp_int> const& m);
  static Vector Zero(UnboundedMatrix<cpp_int> const& m);
};

template<typename Scalar, int rows, int columns>
  requires quantity<Scalar>
struct NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>> {
  using R = FixedUpperTriangularMatrix<Square<Scalar>, columns>;
  using Μ = FixedStrictlyUpperTriangularMatrix<double, columns>;
  using S = FixedVector<Square<Scalar>, columns>;
  using Vector = FixedVector<Scalar, rows>;
  static R UninitializedR(FixedMatrix<Scalar, rows, columns> const& m);
  static Μ UninitializedΜ(FixedMatrix<Scalar, rows, columns> const& m);
  static S UninitializedS(FixedMatrix<Scalar, rows, columns> const& m);
  static Vector Zero(FixedMatrix<Scalar, rows, columns> const& m);
};

template<int rows, int columns>
struct NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>> {
  using R = FixedUpperTriangularMatrix<double, columns>;
  using Μ = FixedStrictlyUpperTriangularMatrix<double, columns>;
  using S = FixedVector<double, columns>;
  using Vector = FixedVector<cpp_int, rows>;
  static R UninitializedR(FixedMatrix<cpp_int, rows, columns> const& m);
  static Μ UninitializedΜ(FixedMatrix<cpp_int, rows, columns> const& m);
  static S UninitializedS(FixedMatrix<cpp_int, rows, columns> const& m);
  static Vector Zero(FixedMatrix<cpp_int, rows, columns> const& m);
};


template<typename Scalar>
auto GramGenerator<UnboundedMatrix<Scalar>>::Uninitialized(
UnboundedMatrix<Scalar> const& m) -> G {
  return G(m.rows(), m.rows(), uninitialized);
}

template<typename Scalar, int rows, int columns>
auto GramGenerator<FixedMatrix<Scalar, rows, columns>>::Uninitialized(
FixedMatrix<Scalar, rows, columns> const& m) -> G {
  return G(uninitialized);
}

template<typename Scalar>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::UninitializedR(
    UnboundedMatrix<Scalar> const& m) -> R {
  return R(m.rows(), m.columns(), uninitialized);
}

template<typename Scalar>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::UninitializedΜ(
    UnboundedMatrix<Scalar> const& m) -> Μ {
  return Μ(m.rows(), m.columns(), uninitialized);
}

template<typename Scalar>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::UninitializedS(
    UnboundedMatrix<Scalar> const& m) -> S {
  return S(m.rows(), uninitialized);
}

template<typename Scalar>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::Zero(
    UnboundedMatrix<Scalar> const& m) -> Vector {
  return Vector(m.rows());
}

inline auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::UninitializedR(
    UnboundedMatrix<cpp_int> const& m) -> R {
  return R(m.columns(), uninitialized);
}

inline auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::UninitializedΜ(
    UnboundedMatrix<cpp_int> const& m) -> Μ {
  return Μ(m.columns(), uninitialized);
}

inline auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::UninitializedS(
    UnboundedMatrix<cpp_int> const& m) -> S {
  return S(m.columns(), uninitialized);
}

inline auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::Zero(
    UnboundedMatrix<cpp_int> const& m) -> Vector {
  return Vector(m.rows());
}

template<typename Scalar, int rows, int columns>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::UninitializedR(
    FixedMatrix<Scalar, rows, columns> const& m) -> R {
  return R(uninitialized);
}

template<typename Scalar, int rows, int columns>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::UninitializedΜ(
    FixedMatrix<Scalar, rows, columns> const& m) -> Μ {
  return Μ(uninitialized);
}

template<typename Scalar, int rows, int columns>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::UninitializedS(
    FixedMatrix<Scalar, rows, columns> const& m) -> S {
  return S(uninitialized);
}

template<typename Scalar, int rows, int columns>
  requires quantity<Scalar>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::Zero(
    FixedMatrix<Scalar, rows, columns> const& m) -> Vector {
  return Vector();
}

template<int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>>::UninitializedR(
    FixedMatrix<cpp_int, rows, columns> const& m) -> R {
  return R(uninitialized);
}

template<int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>>::UninitializedΜ(
    FixedMatrix<cpp_int, rows, columns> const& m) -> Μ {
  return Μ(uninitialized);
}

template<int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>>::UninitializedS(
    FixedMatrix<cpp_int, rows, columns> const& m) -> S {
  return S(uninitialized);
}

template<int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<cpp_int, rows, columns>>::Zero(
    FixedMatrix<cpp_int, rows, columns> const& m) -> Vector {
  return Vector();
}


template<typename Matrix>
typename GramGenerator<Matrix>::G Gram(Matrix const& L) {
  using G = GramGenerator<Matrix>;
  std::int64_t const rows = L.rows();
  std::int64_t const columns = L.columns();
  auto g = G::Uninitialized(L);
  for (std::int64_t i = 0; i < columns; ++i) {
    auto const bᵢ = TransposedView{ColumnView{.matrix = L,
                                              .first_row = 0,
                                              .last_row = rows - 1,
                                              .column = i}};
    for (std::int64_t j = 0; j <= i; ++j) {
      auto const bⱼ = ColumnView{.matrix = L,
                                 .first_row = 0,
                                 .last_row = rows - 1,
                                 .column = j};
      auto const bᵢbⱼ = bᵢ * bⱼ;
      g(i, j) = bᵢbⱼ;
      g(j, i) = bᵢbⱼ;
    }
  }
  return g;
}

// In `NguyễnStehlé` and related functions, the indices for the matrices are
// swapped with respect to [NS09] because their vectors `b` are in rows but ours
// are in columns.  This doesn't matter for the Gram matrix, which is symmetric.

// Moves the column of `b` with index `from_column` to index `to_column`,
// nudging everything after `to_column` by one column (so this is really a
// "rotate").  Updates `G` accordingly.
template<typename Matrix,
         typename GG = GramGenerator<Matrix>>
void Insert(std::int64_t const from_column,
            std::int64_t const to_column,
            Matrix& b,
            typename GG::G& G) {
  CHECK_LT(to_column, from_column);
  std::int64_t const d = b.columns();
  std::int64_t const n = b.rows();

  for (std::int64_t i = 0; i < n; ++i) {
    auto from = std::move(b(i, from_column));
    for (std::int64_t j = from_column; j > to_column; --j) {
      b(i, j) = std::move(b(i, j - 1));
    }
    b(i, to_column) = std::move(from);
  }

  for (std::int64_t i = 0; i < d; ++i) {
    auto from = std::move(G(i, from_column));
    for (std::int64_t j = from_column; j > to_column; --j) {
      G(i, j) = std::move(G(i, j - 1));
    }
    G(i, to_column) = std::move(from);
  }
  for (std::int64_t j = 0; j < d; ++j) {
    auto from = std::move(G(from_column, j));
    for (std::int64_t i = from_column; i > to_column; --i) {
      G(i, j) = std::move(G(i - 1, j));
    }
    G(to_column, j) = std::move(from);
  }

#if _DEBUG
  for (std::int64_t i = 0; i < d ; ++i) {
    auto const bᵢ = ColumnView{.matrix = b,
                               .first_row = 0,
                               .last_row = n - 1,
                               .column = i};
    for (std::int64_t j = 0; j < d; ++j) {
      CHECK_EQ(G(i, j),
               (TransposedView{bᵢ} *                                   // NOLINT
                ColumnView{.matrix = b,
                           .first_row = 0,
                           .last_row = n - 1,
                           .column = j})) << i << " " << j;
    }
  }
#endif
}

// This is [NS09] figure 4, steps 2 to 7.
template<typename Matrix,
         typename GG = GramGenerator<Matrix>,
         typename NSG = NguyễnStehléGenerator<Matrix>>
void CholeskyFactorization(std::int64_t const κ,
                           typename GG::G& G,
                           typename NSG::R& r,
                           typename NSG::Μ& μ,
                           typename NSG::S& s) {
  std::int64_t const i = κ;
  // Step 2.
  for (std::int64_t j = 0; j <= i - 1; ++j) {
    // Step 3.
    r(j, i) = static_cast<double>(G(i, j));
    // Step 4.
    for (std::int64_t k = 0; k <= j - 1; ++k) {
      r(j, i) -= μ(k, j) * r(k, i);
    }
    // Step 5.
    μ(j, i) = r(j, i) / r(j, j);
    // Step 6.
    s[0] = static_cast<double>(G(i, i));
    // Calling this index `j` is a tad confusing, but that's what the article
    // does.
    for (std::int64_t j = 1; j <= i; ++j) {
      s[j] = s[j - 1] - μ(j - 1, i) * r(j - 1, i);
    }
    // Step 7.
    r(i, i) = s[i];
  }
}

// Unless otherwise indicated, this is [NS09] figure 5.
template<typename Matrix,
         typename GG = GramGenerator<Matrix>,
         typename NSG = NguyễnStehléGenerator<Matrix>>
void SizeReduce(std::int64_t const κ,
                Matrix& b,
                typename GG::G& G,
                typename NSG::R& r,
                typename NSG::Μ& μ,
                typename NSG::S& s) {
  // [NS09] figure 7.
  double const η = 0.55;

  std::int64_t const d = b.columns();
  std::int64_t const n = b.rows();
  // Step 1.
  double const η̄ = (η + 0.5) / 2;
  for (;;) {
    // Step 2.
    CholeskyFactorization<Matrix>(κ, G, r, μ, s);
    // Step 3.
    bool terminate = true;
    for (std::int64_t j = 0; j < κ; ++j) {
      if (Abs(μ(j, κ)) > η̄) {
        terminate = false;
        break;
      }
    }
    if (terminate) {
      return;
    }

    std::vector<std::int64_t> X(κ);
    bool making_progress = false;
    for (std::int64_t i = κ - 1; i >= 0; --i) {
      // Step 4.  A large `double` value may overflow `int64_t`, and cause
      // `llround` to return junk.  To avoid this, we make the conversion
      // saturating and trust that the reduction will iterate.
      double const μ_iκ = μ(i, κ);
      if (μ_iκ > largest_convertible_integer) {
        X[i] = largest_convertible_integer;
      } else if (μ_iκ < -largest_convertible_integer) {
        X[i] = -largest_convertible_integer;
      } else {
        X[i] = std::llround(μ_iκ);
      }
      making_progress |= X[i] != 0;
      // Step 5.
      for (std::int64_t j = 0; j <= i - 1; ++j) {
        μ(j, κ) -= X[i] * μ(j, i);
      }
    }
    CHECK(making_progress) << b;
    // Step 6.
    auto b_κ = ColumnView{.matrix = b,
                          .first_row = 0,
                          .last_row = n - 1,
                          .column = κ};
    for (std::int64_t i = 0; i < κ; ++i) {
      auto const bᵢ = typename NSG::Vector(ColumnView{.matrix = b,
                                                      .first_row = 0,
                                                      .last_row = n - 1,
                                                      .column = i});
      // The sum is associated differently from [NS09], which is legitimate
      // because the elements of `b` are of type `cpp_int`.
      b_κ -= X[i] * bᵢ;
    }

    // [NS09], below figure 6.  G is symmetric, so we can write the indices just
    // like in the article.  The article is at best confusing, at worst
    // incorrect.  In particular the summations must stop at `κ - 1` (and not
    // `≠ κ`) and the last sum must have `i < j`.
    // Note that we cannot compute terms like `X[i] * X[j]` as they could
    // overflow.
    typename GG::G::Scalar ΣⱼXⱼ²bⱼ²{};
    typename GG::G::Scalar ΣⱼXⱼbⱼb_κ{};
    typename GG::G::Scalar ΣᵢΣⱼXᵢXⱼbᵢbⱼ{};
    for (std::int64_t j = 0; j < κ; ++j) {
      ΣⱼXⱼ²bⱼ² += X[j] * (X[j] * G(j, j));
      ΣⱼXⱼbⱼb_κ += X[j] * G(j, κ);
      for (std::int64_t i = 0; i < j; ++i) {
        ΣᵢΣⱼXᵢXⱼbᵢbⱼ += X[i] * (X[j] * G(i, j));
      }
    }
    G(κ, κ) += ΣⱼXⱼ²bⱼ² - 2 * ΣⱼXⱼbⱼb_κ + 2 * ΣᵢΣⱼXᵢXⱼbᵢbⱼ;
    DCHECK_EQ(G(κ, κ), TransposedView{b_κ} * b_κ);                     // NOLINT

    for (std::int64_t i = 0; i < d; ++i) {
      if (i != κ) {
        typename GG::G::Scalar ΣⱼXⱼbᵢbⱼ{};
        for (std::int64_t j = 0; j < κ; ++j) {
          ΣⱼXⱼbᵢbⱼ += X[j] * G(i, j);
        }
        G(i, κ) -= ΣⱼXⱼbᵢbⱼ;
        G(κ, i) = G(i, κ);
        DCHECK_EQ(G(i, κ),
                  (TransposedView{ColumnView{.matrix = b,              // NOLINT
                                             .first_row = 0,
                                             .last_row = n - 1,
                                             .column = i}} *
                   b_κ)) << i;
      }
    }
  }
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

// Unless otherwise indicated, this is [NS09] figure 9.
template<typename Matrix>
  requires two_dimensional<Matrix>
Matrix NguyễnStehlé(Matrix const& L) {
  // [NS09] figure 7.
  double const ẟ = 0.75;

  using NSG = NguyễnStehléGenerator<Matrix>;
  auto b = L;
  std::int64_t const d = b.columns();
  std::int64_t const n = b.rows();
  auto const zero = NSG::Zero(b);

  // Step 1.
  auto G = Gram(b);
  // Step 2.
  // Note that the combining macron doesn't work well here so we use a modifier
  // macron.
  double const δ̄ = (ẟ + 1) / 2;
  typename NSG::R r = NSG::UninitializedR(b);
  typename NSG::Μ μ = NSG::UninitializedΜ(b);
  typename NSG::S s = NSG::UninitializedS(b);
  r(0, 0) = static_cast<typename NSG::R::Scalar>(G(0, 0));
  std::int64_t κ = 1;
  std::int64_t ζ = -1;
  while (κ < d) {
    // Step 3.
    SizeReduce(κ, b, G, r, μ, s);
    // Step 4.
    std::int64_t κʹ = κ;
    while (κ >= ζ + 2 && δ̄ * r(κ - 1, κ - 1) >= s[κ - 1]) {
      --κ;
    }
    // Step 5.
    for (std::int64_t i = ζ + 1; i <= κ - 1; ++i) {
      μ(i, κ) = μ(i, κʹ);
      r(i, κ) = r(i, κʹ);
    }
    r(κ, κ) = s[κ];
    // Step 6.
    if (κʹ != κ) {
      Insert(/*from_column=*/κʹ, /*to_column=*/κ, b, G);
    }
    // Step 7.
    auto const bκ = ColumnView{.matrix = b,
                               .first_row = 0,
                               .last_row = n - 1,
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
