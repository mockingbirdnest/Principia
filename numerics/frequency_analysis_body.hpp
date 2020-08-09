
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <functional>

#include "numerics/fixed_arrays.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using quantities::Inverse;
using quantities::Square;
using quantities::SquareRoot;
namespace si = quantities::si;

// A helper struct for generating the Poisson series tⁿ sin ω t and tⁿ cos ω t.
template<typename Series, int n>
struct SeriesGenerator {
  // The series tⁿ sin ω t.
  static Series Sin(AngularFrequency const& ω, Instant const& origin);
  // The series tⁿ cos ω t.
  static Series Cos(AngularFrequency const& ω, Instant const& origin);

 private:
  // The polynomial tⁿ.
  static typename Series::Polynomial Unit(Instant const& origin);
};

// A helper struct for generating the Кудрявцев basis, i.e., functions of the
// form tⁿ sin ω t and tⁿ cos ω t properly ordered.
template<typename Series,
         typename = std::make_index_sequence<Series::degree>>
struct BasisGenerator;

template<typename Series, std::size_t... indices>
struct BasisGenerator<Series, std::index_sequence<indices...>> {
  //TODO(phl): omega=0
  static std::array<Series, 2 * Series::degree> Basis(AngularFrequency const& ω,
                                                      Instant const& origin);
};


template<typename Series, int n>
Series SeriesGenerator<Series, n>::Sin(AngularFrequency const& ω,
                                       Instant const& origin) {
  typename Series::Polynomial::Coefficients const zeros;
  typename Series::Polynomial const zero{zeros, origin};
  return Series(zero,
                {{ω,
                  {/*sin=*/Unit(origin),
                   /*cos=*/zero}}});
}

template<typename Series, int n>
Series SeriesGenerator<Series, n>::Cos(AngularFrequency const& ω,
                                       Instant const& origin) {
  typename Series::Polynomial::Coefficients const zeros;
  typename Series::Polynomial const zero{zeros, origin};
  return Series(zero,
                {{ω,
                  {/*sin=*/zero,
                   /*cos=*/Unit(origin)}}});
}

template<typename Series, int n>
typename Series::Polynomial SeriesGenerator<Series, n>::Unit(
    Instant const& origin) {
  typename Series::Polynomial::Coefficients coefficients;
  std::get<n>(coefficients) = si::Unit<
      std::tuple_element_t<n, typename Series::Polynomial::Coefficients>>;
  return Series::Polynomial(coefficients, origin);
}

template<typename Series, std::size_t... indices>
std::array<Series, 2 * Series::degree>
BasisGenerator<Series, std::index_sequence<indices...>>::Basis(
    AngularFrequency const& ω,
    Instant const& origin) {
  // This has the elements {Sin(ωt), t Sin(ωt), t² Sin(ωt), ..., Cos(ωt), ...}
  // which is not the order we want (we want lower-degree polynomials first).
  std::array<Series, 2 * Series::degree> all_series = {
      SeriesGenerator<Series, indices>::Sin(ω, origin)...,
      SeriesGenerator<Series, indices>::Cos(ω, origin)...};

  // Order all_series by repeatedly swapping its elements.
  if (Series::degree >= 2) {
    // The index of this array is the current index of a series in all_series.
    // The value is the index of the final resting place of that series in
    // all_series.  The elements at indices 0 and 2 * Series::degree - 1 are
    // unused.
    std::array<int, 2 * Series::degree> permutation;
    for (int i = 1; i < 2 * Series::degree - 1; ++i) {
      permutation[i] =
          i < Series::degree ? 2 * i : 2 * (i - Series::degree) + 1;
    }
    for (int i = 1; i < 2 * Series::degree - 1;) {
      // Swap the series currently at index i to its final resting place.
      // Iterate until the series at index i is at its final resting place
      // (i.e., after we have executed an entire cycle of the permutation).
      // Then move to the next series.
      if (i == permutation[i]) {
        ++i;
      } else {
        int const j = permutation[i];
        std::swap(all_series[i], all_series[j]);
        std::swap(permutation[i], permutation[j]);
      }
    }
  }
  return all_series;
}


template<typename Function,
         typename RValue, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    DotProduct<Function, RValue, rdegree_, wdegree_, Evaluator> const& dot) {
  using DotResult = Product<std::invoke_result_t<Function, Instant>, RValue>;
  using Degree0 = PoissonSeries<double, 0, Evaluator>;

  auto amplitude = [&dot, &function, &weight](AngularFrequency const& ω) {
    Instant const& t0 = weight.origin();
    Degree0 const sin(typename Degree0::Polynomial({0}, t0),
                      {{ω,
                        {/*sin=*/typename Degree0::Polynomial({1}, t0),
                         /*cos=*/typename Degree0::Polynomial({0}, t0)}}});
    Degree0 const cos(typename Degree0::Polynomial({0}, t0),
                      {{ω,
                        {/*sin=*/typename Degree0::Polynomial({0}, t0),
                         /*cos=*/typename Degree0::Polynomial({1}, t0)}}});
    auto const sin_amplitude = dot(function, sin, weight);
    auto const cos_amplitude = dot(function, cos, weight);
    return sin_amplitude * sin_amplitude + cos_amplitude * cos_amplitude;
  };

  return GoldenSectionSearch(amplitude,
                             fft_mode.min,
                             fft_mode.max,
                             std::greater<Square<DotResult>>());
}

template<typename Function,
         typename RValue, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, rdegree_, Evaluator>
Projection(
    AngularFrequency const& ω,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    DotProduct<Function, RValue, rdegree_, wdegree_, Evaluator> const& dot) {
  using Value = std::invoke_result_t<Function, Instant>;
  using DotProductResult = Product<Value, RValue>;
  using Series = PoissonSeries<Value, degree_, Evaluator>;

  Instant const& t0 = weight.origin();
  auto const basis = BasisGenerator<Series>::Basis(ω, t0);
  constexpr int basis_size = std::tuple_size_v<decltype(basis);

  // Our indices start at 0, unlike those of Кудрявцев which start at 1.
  //TODO(phl):Check the indices!
  std::array<DotProductResult, basis_size> F;
  FixedLowerTriangularMatrix<DotProductResult, basis_size> Q;
  FixedLowerTriangularMatrix<Inverse<SquareRoot<DotProductResult>>,
                             basis_size> α;
  FixedLowerTriangularMatrix<double, basis_size> A;
  FixedStrictlyLowerTriangularMatrix<SquareRoot<DotProductResult>,
                                     basis_size> B;
  std::array<Function, basis_size> f;

  F[0] = dot(function, basis[0], weight);
  Q[0][0] = dot(basis[0], basis[0], weight);
  α[0][0] = 1 / Sqrt(Q[0][0]);
  A[0][0] = α[0][0] * α[0][0] * F[0];
  f[0] = function - A[0][0] * basis[0];
  for (int m = 1; m < basis.size(); ++m) {
    F[m] = dot(f[m - 1], basis[m], weight);
    for (int j = 0; j <= m; ++j) {
      Q[m][j] = dot(basis[m], basis[j], weight);
    }

    for (int j = 0; j < m; ++j) {
      SquareRoot<DotProductResult> accumulator;
      for (int s = 0; s <= j; ++s) {
        accumulator -= α[j][s] * Q[m][s];
      }
      B[j][m] = accumulator;
    }

    {
      DotProductResult accumulator = Q[m][m];
      for (int s = 0; s < m; ++s) {
        accumulator -= B[s][m] * B[s][m];
      }
      α[m][m] = 1 / Sqrt(accumulator);
    }

    for (int j = 0; j < m; ++j) {
      double accumulator = 0;
      for (int s = j; s < m; ++s) {
        accumulator += B[s][m] * α[s][j]
      }
      α[m][j] = α[m][m] * accumulator;
    }

    A[m][m] = α[m][m] * α[m][m] * F[m];

    for (int j = 0; j < m; ++j) {
      A[j][m] = A[j][m - 1] + α[m][m] * α[m][j] * F[m];
    }

    {
      PoissonSeries<double, degree_, Evaluator> accumulator =
          α[m][0] * basis[0];
      for (int i = 1; i <= m; ++i) {
        accumulator += α[m][i] * basis[i];
      }
      f[m] = f[m - 1] - α[m][m] * F[m] * accumulator;
    }
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
