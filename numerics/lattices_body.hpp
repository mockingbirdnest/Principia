#pragma once

#include "numerics/lattices.hpp"

#include <algorithm>

#include "base/tags.hpp"
#include "boost/multiprecision/cpp_int.hpp"
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

using namespace boost::multiprecision;
using namespace principia::base::_tags;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

// In the terminology of [NS09], our vectors are in columns, so |d| is |columns|
// and |n| is |rows|.
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
struct NguyễnStehléGenerator<UnboundedMatrix<Scalar>> {
  using R = UnboundedUpperTriangularMatrix<Square<Scalar>>;
  using Μ = UnboundedUpperTriangularMatrix<double>;  //TODO(phl)Strictly
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
  using Μ = UnboundedUpperTriangularMatrix<double>;  //TODO(phl)Strictly
  using S = UnboundedVector<double>;
  using Vector = UnboundedVector<cpp_int>;
  static R UninitializedR(UnboundedMatrix<cpp_int> const& m);
  static Μ UninitializedΜ(UnboundedMatrix<cpp_int> const& m);
  static S UninitializedS(UnboundedMatrix<cpp_int> const& m);
  static Vector Zero(UnboundedMatrix<cpp_int> const& m);
};

template<typename Scalar, int rows, int columns>
struct NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>> {
  using R = FixedUpperTriangularMatrix<Square<Scalar>, columns>;
  using Μ = FixedUpperTriangularMatrix<double, columns>;  //TODO(phl)Strictly
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
  using Μ = FixedUpperTriangularMatrix<double, columns>;  //TODO(phl)Strictly
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
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::UninitializedR(
    UnboundedMatrix<Scalar> const& m) -> R {
  return R(m.rows(), m.columns(), uninitialized);
}

template<typename Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::UninitializedΜ(
    UnboundedMatrix<Scalar> const& m) -> Μ {
  return Μ(m.rows(), m.columns(), uninitialized);
}

template<typename Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::UninitializedS(
    UnboundedMatrix<Scalar> const& m) -> S {
  return S(m.rows(), uninitialized);
}

template<typename Scalar>
auto NguyễnStehléGenerator<UnboundedMatrix<Scalar>>::Zero(
    UnboundedMatrix<Scalar> const& m) -> Vector {
  return Vector(m.rows());
}

auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::UninitializedR(
    UnboundedMatrix<cpp_int> const& m) -> R {
  return R(m.columns(), uninitialized);
}

auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::UninitializedΜ(
    UnboundedMatrix<cpp_int> const& m) -> Μ {
  return Μ(m.columns(), uninitialized);
}

auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::UninitializedS(
    UnboundedMatrix<cpp_int> const& m) -> S {
  return S(m.columns(), uninitialized);
}

auto NguyễnStehléGenerator<UnboundedMatrix<cpp_int>>::Zero(
    UnboundedMatrix<cpp_int> const& m) -> Vector {
  return Vector(m.rows());
}

template<typename Scalar, int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::UninitializedR(
    FixedMatrix<Scalar, rows, columns> const& m) -> R {
  return R(uninitialized);
}

template<typename Scalar, int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::UninitializedΜ(
    FixedMatrix<Scalar, rows, columns> const& m) -> Μ {
  return Μ(uninitialized);
}

template<typename Scalar, int rows, int columns>
auto NguyễnStehléGenerator<FixedMatrix<Scalar, rows, columns>>::UninitializedS(
    FixedMatrix<Scalar, rows, columns> const& m) -> S {
  return S(uninitialized);
}

template<typename Scalar, int rows, int columns>
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


//TODO(phl)comment
template<typename Matrix,
         typename GG = GramGenerator<Matrix>>
void Insert(std::int64_t const from_column,
            std::int64_t const to_column,
            Matrix& b,
            typename GG::G& G) {
  CHECK_LT(to_column, from_column);

  //LOG(ERROR)<<"from "<<from_column<<" to "<<to_column;
  //LOG(ERROR)<<"b before "<<b;

  for (std::int64_t i = 0; i < b.rows(); ++i) {
    auto from = std::move(b(i, from_column));
    for (std::int64_t j = from_column; j > to_column; --j) {
      b(i, j) = std::move(b(i, j - 1));
    }
    b(i, to_column) = std::move(from);
  }
  //LOG(ERROR)<<"b after "<<b;

  std::int64_t const to_row = to_column;
  std::int64_t const from_row = from_column;

  //LOG(ERROR)<<"G before "<<G;

  // Squirrel away the row `from_row` (the last one) since it will be
  // overwritten below.
  std::vector<typename GG::G::Scalar> t;
  t.reserve(from_column - to_column + 1);
  for (std::int64_t j = to_column; j <= from_column; ++j) {
    t.push_back(std::move(G(from_row, j)));
  }

  // Move the bulk of the elements.
  for (std::int64_t i = from_row; i > to_row; --i) {
    auto from = std::move(G(i - 1, from_column));
    for (std::int64_t j = from_column; j > to_column; --j) {
      G(i, j) = std::move(G(i - 1, j - 1));
    }
    G(i, to_column) = std::move(from);
  }

  // Restore the row that we saved at the beginning.  It now goes to `to_row`
  // (the first one) with a suitable rotation.
  for (std::int64_t j = to_column + 1; j <= from_column; ++j) {
    G(to_row, j) = std::move(t[j - to_column - 1]);
  }
  G(to_row, to_column) = std::move(t[t.size() - 1]);

  //LOG(ERROR)<<"G after "<<G;
}

// This is [NS09] figure 4.
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
  for (std::int64_t j = 0; j < i; ++j) {
    // Step 3.
    r(j, i) = static_cast<double>(G(i, j));
    // Step 4.
    for (std::int64_t k = 0; k < j; ++k) {
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

  std::int64_t const rows = b.rows();
  // Step 1.
  double const ηˉ = (η + 0.5) / 2;
  for (;;) {
    // Step 2.
    CholeskyFactorization<Matrix>(κ, G, r, μ, s);
LOG(ERROR)<<"b:\n"<<b;
LOG(ERROR)<<"G:\n"<<G;
LOG(ERROR)<<"r:\n"<<r;
LOG(ERROR)<<"μ:\n"<<μ;
LOG(ERROR)<<"s:\n"<<s;
    // Step 3.
    bool terminate = true;
    for (std::int64_t j = 0; j < κ; ++j) {
      if (Abs(μ(j, κ)) > ηˉ) {
        terminate = false;
        break;
      }
    }
    if (terminate) {
      return;
    }

    std::vector<std::int64_t> X(κ);
    for (std::int64_t i = κ - 1; i >= 0; --i) {
      // Step 4.
      X[i] = std::llround(μ(i, κ));
      // Step 5.
      for (std::int64_t j = 0; j < i; ++j) {
        μ(j, κ) -= X[i] * μ(j, i);
      }
    }
    // Step 6.
    auto b_κ = ColumnView{.matrix = b,
                          .first_row = 0,
                          .last_row = rows - 1,
                          .column = κ};
LOG(ERROR)<<"kappa: "<<κ;
LOG(ERROR)<<"before: "<<b_κ;
    for (std::int64_t i = 0; i < κ; ++i) {
      auto const bᵢ = typename NSG::Vector(ColumnView{.matrix = b,
                                                      .first_row = 0,
                                                      .last_row = rows - 1,
                                                      .column = i});
      b_κ -= X[i] * bᵢ;
LOG(ERROR)<<"i: "<<i<<" Xi: "<<X[i]<<" b: "<<b_κ;
    }
LOG(ERROR)<<"after: "<<b_κ;

    // [NS09], below figure 6.  Note that G is symmetric, so we can write the
    // indices just like in the paper.
    // Note that we cannot compute terms like X[i] * X[j] as they could
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
LOG(ERROR)<<G(κ, κ);
LOG(ERROR)<<ΣⱼXⱼ²bⱼ²;
LOG(ERROR)<<ΣⱼXⱼbⱼb_κ;
LOG(ERROR)<<ΣᵢΣⱼXᵢXⱼbᵢbⱼ;
    G(κ, κ) += ΣⱼXⱼ²bⱼ² - 2 * ΣⱼXⱼbⱼb_κ + 2 * ΣᵢΣⱼXᵢXⱼbᵢbⱼ;
    DCHECK_EQ(G(κ, κ), TransposedView{b_κ} * b_κ);

    for (std::int64_t i = 0; i < κ; ++i) {
      typename GG::G::Scalar ΣⱼXⱼbᵢbⱼ{};
      for (std::int64_t j = 0; j < κ; ++j) {
        ΣⱼXⱼbᵢbⱼ += X[j] * G(i, j);
      }
      G(i, κ) -= ΣⱼXⱼbᵢbⱼ;
      G(κ, i) = G(i, κ);
      DCHECK_EQ(G(i, κ),
                (TransposedView{ColumnView{.matrix = b,
                                           .first_row = 0,
                                           .last_row = rows - 1,
                                           .column = i}} *
                 b_κ));
    }
LOG(ERROR)<<"after: "<<G;
  }
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
  double const δˉ = (ẟ + 1) / 2;
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
    //TODO(phl)high index probably useless
    std::int64_t κʹ = κ;
LOG(ERROR)<<"Lovacs: "<<r(κ - 1, κ - 1) << " vs " << s[κ - 1];
    while (κ >= ζ + 2 && δˉ * r(κ - 1, κ - 1) >= s[κ - 1]) {
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
