#pragma once

#include "numerics/poisson_series_basis.hpp"

#include <algorithm>

#include "quantities/si.hpp"
#include "quantities/traits.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series_basis {

using quantities::is_quantity_v;
namespace si = quantities::si;

// A helper struct for generating the Poisson series tⁿ sin ω t and tⁿ cos ω t.
// d is either 0 (for a 1-dimensional value type) or 0, 1, 2 (for a
// 3-dimensional value type).
template<typename Series, int n, int d>
struct SeriesGenerator {
  // The series tⁿ.
  static Series Aperiodic(Instant const& origin);
  // The series tⁿ sin ω t.
  static Series Sin(AngularFrequency const& ω, Instant const& origin);
  // The series tⁿ cos ω t.
  static Series Cos(AngularFrequency const& ω, Instant const& origin);

 private:
  // The polynomial tⁿ.
  static typename Series::Polynomial Unit(Instant const& origin);
};

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Aperiodic(Instant const& origin) {
  return Series(Unit(origin), {});
}

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Sin(AngularFrequency const& ω,
                                          Instant const& origin) {
  typename Series::Polynomial::Coefficients const zeros;
  typename Series::Polynomial const zero{zeros, origin};
  return Series(zero,
                {{ω,
                  {/*sin=*/Unit(origin),
                   /*cos=*/zero}}});
}

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Cos(AngularFrequency const& ω,
                                          Instant const& origin) {
  typename Series::Polynomial::Coefficients const zeros;
  typename Series::Polynomial const zero{zeros, origin};
  return Series(zero,
                {{ω,
                  {/*sin=*/zero,
                   /*cos=*/Unit(origin)}}});
}

template<typename Series, int n, int d>
typename Series::Polynomial SeriesGenerator<Series, n, d>::Unit(
    Instant const& origin) {
  typename Series::Polynomial::Coefficients coefficients;
  using Coefficient =
      std::tuple_element_t<n, typename Series::Polynomial::Coefficients>;
  Coefficient& coefficient = std::get<n>(coefficients);
  if constexpr (is_quantity_v<Coefficient>) {
    coefficient = si::Unit<Coefficient>;
  } else {
    using Scalar = typename Hilbert<Coefficient>::NormType;
    if constexpr (d == 0) {
      coefficient = Coefficient({si::Unit<Scalar>, Scalar{}, Scalar{}});
    } else if constexpr (d == 1) {
      coefficient = Coefficient({Scalar{}, si::Unit<Scalar>, Scalar{}});
    } else if constexpr (d == 2) {
      coefficient = Coefficient({Scalar{}, Scalar{}, si::Unit<Scalar>});
    }
  }
  return typename Series::Polynomial(coefficients, origin);
}

template<typename Series, int dimension, std::size_t... indices>
std::array<Series, dimension*(Series::degree + 1)> PoissonSeriesBasisGenerator<
    Series, dimension,
    std::index_sequence<indices...>>::Basis(Instant const& origin) {
  return {SeriesGenerator<Series, indices / dimension, indices % dimension>::
              Aperiodic(origin)...};
}

template<typename Series, int dimension, std::size_t... indices>
std::array<Series, 2 * dimension*(Series::degree + 1)>
PoissonSeriesBasisGenerator<
    Series, dimension,
    std::index_sequence<indices...>>::Basis(AngularFrequency const& ω,
                                            Instant const& origin) {
  // This has the elements {Sin(ωt), t Sin(ωt), t² Sin(ωt), ..., Cos(ωt), ...}
  // in the scalar case and {x Sin(ωt), y Sin(ωt), z Sin(ωt), x t Sin(ωt), ...}
  // in the vector case.  This is not the order we want (we want lower-degree
  // polynomials first) so we'll need to reorder the terms.
  std::array<Series, 2 * dimension * (Series::degree + 1)> all_series = {
      SeriesGenerator<Series, indices / dimension, indices % dimension>::Sin(
          ω, origin)...,
      SeriesGenerator<Series, indices / dimension, indices % dimension>::Cos(
          ω, origin)...};

  // Order all_series by repeatedly swapping its elements.
  if constexpr (all_series.size() > 2) {
    // The index of this array is the current index of a series in all_series.
    // The value is the index of the final resting place of that series in
    // all_series.  The elements at indices 0 and
    // 2 * dimension * (Series::degree + 1) are unused.
    std::array<int, all_series.size()> permutation;
    for (int i = 1; i < all_series.size() - 1; ++i) {
      permutation[i] = i < dimension * (Series::degree + 1)
                           ? 2 * i
                           : 2 * (i - dimension * (Series::degree + 1)) + 1;
    }
    for (int i = 1; i < all_series.size() - 1;) {
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

}  // namespace internal_poisson_series_basis
}  // namespace numerics
}  // namespace principia
