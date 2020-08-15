
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "numerics/fixed_arrays.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using quantities::Inverse;
using quantities::Sqrt;
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
         typename = std::make_index_sequence<Series::degree + 1>>
struct BasisGenerator;

template<typename Series, std::size_t... indices>
struct BasisGenerator<Series, std::index_sequence<indices...>> {
  // TODO(phl): Will need to properly handle ω = 0.
  static std::array<Series, 2 * Series::degree + 2> Basis(
      AngularFrequency const& ω,
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
  return typename Series::Polynomial(coefficients, origin);
}

template<typename Series, std::size_t... indices>
std::array<Series, 2 * Series::degree + 2>
BasisGenerator<Series, std::index_sequence<indices...>>::Basis(
    AngularFrequency const& ω,
    Instant const& origin) {
  // This has the elements {Sin(ωt), t Sin(ωt), t² Sin(ωt), ..., Cos(ωt), ...}
  // which is not the order we want (we want lower-degree polynomials first).
  std::array<Series, 2 * Series::degree + 2> all_series = {
      SeriesGenerator<Series, indices>::Sin(ω, origin)...,
      SeriesGenerator<Series, indices>::Cos(ω, origin)...};

  // Order all_series by repeatedly swapping its elements.
  if (Series::degree >= 2) {
    // The index of this array is the current index of a series in all_series.
    // The value is the index of the final resting place of that series in
    // all_series.  The elements at indices 0 and 2 * Series::degree + 1 are
    // unused.
    std::array<int, 2 * Series::degree + 2> permutation;
    for (int i = 1; i < 2 * Series::degree + 1; ++i) {
      permutation[i] =
          i <= Series::degree ? 2 * i : 2 * (i - Series::degree) - 1;
    }
    for (int i = 1; i < 2 * Series::degree + 1;) {
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
         int wdegree_,
         typename Dot,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Dot const& dot) {
  using Value = std::invoke_result_t<Function, Instant>;
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
                             std::greater<Square<Value>>());
}

template<int degree_,
         typename Function,
         int wdegree_, typename Dot,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
Projection(AngularFrequency const& ω,
           Function const& function,
           PoissonSeries<double, wdegree_, Evaluator> const& weight,
           Dot const& dot) {
  std::optional optional_ω = ω;

  // A calculator that returns optional_ω once and then stops.
  auto angular_frequency_calculator = [&optional_ω](auto const& residual) {
    auto const result = optional_ω;
    optional_ω = std::nullopt;
    return result;
  };

  return IncrementalProjection<degree_>(function,
                                        angular_frequency_calculator,
                                        weight,
                                        dot);
}

template<int degree_,
         typename Function,
         typename AngularFrequencyCalculator, int wdegree_, typename Dot,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double, wdegree_, Evaluator> const& weight,
                      Dot const& dot) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Series = PoissonSeries<Value, degree_, Evaluator>;

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<Series> basis;

  auto const ω_basis = BasisGenerator<Series>::Basis(ω.value(), t0);
  int basis_size = std::tuple_size_v<decltype(ω_basis)>;
  std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));

  FixedLowerTriangularMatrix<Inverse<Value>, basis_size> α;

  // Only indices 0 to m - 1 are used in this array.  At the beginning of
  // iteration m it contains Aⱼ⁽ᵐ⁻¹⁾.
  std::vector<double> A;
  A.resize(basis_size, 0);

  auto const F₀ = dot(function, basis[0], weight);
  auto const Q₀₀ = dot(basis[0], basis[0], weight);
  α[0][0] = 1 / Sqrt(Q₀₀);
  A[0] = F₀ / Q₀₀;

  // At the beginning of iteration m this contains fₘ₋₁.
  auto f = function - A[0] * basis[0];

  int m_begin = 1;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      // Contains Fₘ.
      auto const F = dot(f, basis[m], weight);

      // Only indices 0 to m are used in this array.  It contains Qₘⱼ.
      std::array<Square<Value>, basis_size> Q;
      for (int j = 0; j <= m; ++j) {
        Q[j] = dot(basis[m], basis[j], weight);
      }

      // Only indices 0 to m - 1 are used in this array.  It contains Bⱼ⁽ᵐ⁾.
      std::array<Value, basis_size> B;
      for (int j = 0; j < m; ++j) {
        Value Σ_αⱼₛ_Qₘₛ{};
        for (int s = 0; s <= j; ++s) {
          Σ_αⱼₛ_Qₘₛ += α[j][s] * Q[s];
        }
        B[j] = -Σ_αⱼₛ_Qₘₛ;
      }

      {
        Square<Value> Σ_Bₛ⁽ᵐ⁾²{};
        for (int s = 0; s < m; ++s) {
          Σ_Bₛ⁽ᵐ⁾² += B[s] * B[s];
        }
        DCHECK_LE(Σ_Bₛ⁽ᵐ⁾², Q[m]);
        α[m][m] = 1 / Sqrt(Q[m] - Σ_Bₛ⁽ᵐ⁾²);
      }

      for (int j = 0; j < m; ++j) {
        double Σ_Bₛ⁽ᵐ⁾_αₛⱼ = 0;
        for (int s = j; s < m; ++s) {
          Σ_Bₛ⁽ᵐ⁾_αₛⱼ += B[s] * α[s][j];
        }
        α[m][j] = α[m][m] * Σ_Bₛ⁽ᵐ⁾_αₛⱼ;
      }

      A[m] = α[m][m] * α[m][m] * F;

      for (int j = 0; j < m; ++j) {
        A[j] += α[m][m] * α[m][j] * F;
      }

      {
        PoissonSeries<double, degree_, Evaluator> Σ_αₘᵢ_eᵢ = α[m][0] * basis[0];
        for (int i = 1; i <= m; ++i) {
          Σ_αₘᵢ_eᵢ += α[m][i] * basis[i];
        }
        f -= α[m][m] * F * Σ_αₘᵢ_eᵢ;
      }
    }

    PoissonSeries<Value, degree_, Evaluator> result = A[0] * basis[0];
    for (int i = 1; i < basis_size; ++i) {
      result += A[i] * basis[i];
    }

    ω = calculator(f);
    if (!ω.has_value()) {
      return result;
    }

    m_begin = basis_size;
    auto const ω_basis = BasisGenerator<Series>::Basis(ω.value(), t0);
    basis_size += std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    α.Extend(basis_size);
    A.resize(basis_size, 0);
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
