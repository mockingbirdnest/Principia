
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/poisson_series_basis.hpp"
#include "numerics/root_finders.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using base::uninitialized;
using geometry::Hilbert;
using geometry::Vector;
using quantities::Inverse;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;

template<typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight) {
  auto const weighted_function = weight * function;
  auto const weighted_function_spectrum = weighted_function.FourierTransform();

  auto power =
      [&weighted_function_spectrum](AngularFrequency const& ω) {
        return weighted_function_spectrum(ω).Norm²();
      };

  return Brent(power,
               fft_mode.min,
               fft_mode.max,
               std::greater<>());
}

template<int degree_,
         typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
Projection(Function const& function,
           AngularFrequency const& ω,
           PoissonSeries<double, wdegree_, Evaluator> const& weight,
           Instant const& t_min,
           Instant const& t_max) {
  std::optional<AngularFrequency> optional_ω = ω;

  // A calculator that returns optional_ω once and then stops.
  auto angular_frequency_calculator = [&optional_ω](auto const& residual) {
    auto const result = optional_ω;
    optional_ω = std::nullopt;
    return result;
  };

  return IncrementalProjection<degree_>(function,
                                        angular_frequency_calculator,
                                        weight,
                                        t_min, t_max);
}

template<int degree_,
         typename Function,
         typename AngularFrequencyCalculator, int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double, wdegree_, Evaluator> const& weight,
                      Instant const& t_min,
                      Instant const& t_max) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Norm = typename Hilbert<Value>::NormType;
  using Norm² = typename Hilbert<Value>::InnerProductType;
  using Normalized = typename Hilbert<Value>::NormalizedType;
  using Series = PoissonSeries<Value, degree_, Evaluator>;

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<Series> basis;

  int basis_size;
  if (ω.value() == AngularFrequency{}) {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
            t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
  } else {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
            ω.value(), t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
  }

  UnboundedUpperTriangularMatrix<Norm> r(basis_size, uninitialized);

  std::vector<PoissonSeries<Normalized, degree_, Evaluator>> q;

  auto const a₀ = basis[0];
  r[0][0] = a₀.Norm(weight, t_min, t_max);
  q.push_back(a₀ / r[0][0]);

  auto const A₀ = InnerProduct(function, q[0], weight, t_min, t_max);

  PoissonSeries<Value, degree_, Evaluator> F = A₀ * q[0];
  auto f = function - F;

  int m_begin = 1;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      //a^(k)_n+1
      auto aₘ⁽ᵏ⁾ = basis[m];
      for (int k = 0; k < m; ++k) {
        r[k][m] = InnerProduct(q[k], aₘ⁽ᵏ⁾, weight, t_min, t_max);
        aₘ⁽ᵏ⁾ -= r[k][m] * q[k];
      }

      r[m][m] = aₘ⁽ᵏ⁾.Norm(weight, t_min, t_max);
      q.push_back(aₘ⁽ᵏ⁾ / r[m][m]);
      DCHECK_EQ(m + 1, q.size());

      auto A = InnerProduct(f, q[m], weight, t_min, t_max);

      f -= A * q[m];
      F += A * q[m];
    }

    ω = calculator(f);
    if (!ω.has_value()) {
      return F;
    }

    int ω_basis_size;
    if (ω.value() == AngularFrequency{}) {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
              t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    } else {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
              ω.value(), t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    }
    r.Extend(ω_basis_size, uninitialized);
    m_begin = basis_size;
    basis_size += ω_basis_size;
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
