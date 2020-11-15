#pragma once

#include "numerics/poisson_series_basis.hpp"

#include <algorithm>

#include "geometry/hilbert.hpp"
#include "quantities/si.hpp"
#include "quantities/traits.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series_basis {

using geometry::Hilbert;
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
  template<typename Polynomial>
  static Polynomial Unit(Instant const& origin);
};

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Aperiodic(Instant const& origin) {
  return Series(Unit<typename Series::AperiodicPolynomial>(origin), {});
}

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Sin(AngularFrequency const& ω,
                                          Instant const& origin) {
  typename Series::AperiodicPolynomial const aperiodic_zero{{}, origin};
  typename Series::PeriodicPolynomial const periodic_zero{{}, origin};
  return Series(aperiodic_zero,
                {{ω,
                  {/*sin=*/Unit<typename Series::PeriodicPolynomial>(origin),
                   /*cos=*/periodic_zero}}});
}

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Cos(AngularFrequency const& ω,
                                          Instant const& origin) {
  typename Series::AperiodicPolynomial const aperiodic_zero{{}, origin};
  typename Series::PeriodicPolynomial const periodic_zero{{}, origin};
  return Series(
      aperiodic_zero,
      {{ω,
        {/*sin=*/periodic_zero,
         /*cos=*/Unit<typename Series::PeriodicPolynomial>(origin)}}});
}

template<typename Series, int n, int d>
template<typename Polynomial>
Polynomial SeriesGenerator<Series, n, d>::Unit(Instant const& origin) {
  typename Polynomial::Coefficients coefficients;
  using Coefficient =
      std::tuple_element_t<n, typename Polynomial::Coefficients>;
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
  return Polynomial(coefficients, origin);
}

inline bool PoissonSeriesSubspace::orthogonal(PoissonSeriesSubspace const v,
                                              PoissonSeriesSubspace const w) {
  return v.coordinate_ != w.coordinate_;
}

inline PoissonSeriesSubspace::PoissonSeriesSubspace(Coordinate coordinate,
                                                    Parity parity)
    : coordinate_(coordinate), parity_(parity) {}

template<typename Series,
         int dimension,
         int degree,
         std::size_t... indices,
         std::size_t... periodic_indices>
std::array<Series, dimension*(degree + 1)> PoissonSeriesBasisGenerator<
    Series,
    dimension,
    degree,
    std::index_sequence<indices...>,
    std::index_sequence<periodic_indices...>>::Basis(Instant const& origin) {
  return {SeriesGenerator<Series, indices / dimension, indices % dimension>::
              Aperiodic(origin)...};
}

template<typename Series,
         int dimension,
         int degree,
         std::size_t... indices,
         std::size_t... periodic_indices>
inline std::array<PoissonSeriesSubspace, dimension*(degree + 1)>
PoissonSeriesBasisGenerator<Series,
                            dimension,
                            degree,
                            std::index_sequence<indices...>,
                            std::index_sequence<periodic_indices...>>::
    Subspaces(Instant const& origin) {
  return {PoissonSeriesSubspace{
      static_cast<PoissonSeriesSubspace::Coordinate>(indices % dimension),
      static_cast<PoissonSeriesSubspace::Parity>(indices / dimension)}...};
}

template<typename Series,
         int dimension,
         int degree,
         std::size_t... indices,
         std::size_t... periodic_indices>
std::array<Series, 2 * dimension * (degree + 1)> PoissonSeriesBasisGenerator<
    Series,
    dimension,
    degree,
    std::index_sequence<indices...>,
    std::index_sequence<periodic_indices...>>::Basis(AngularFrequency const& ω,
                                                     Instant const& origin) {
  // This has the elements {Sin(ωt), t Sin(ωt), t² Sin(ωt), ..., Cos(ωt), ...}
  // in the scalar case and {x Sin(ωt), y Sin(ωt), z Sin(ωt), x t Sin(ωt), ...}
  // in the vector case.  This is not the order we want (we want lower-degree
  // polynomials first) so we'll need to reorder the terms.
  std::array<Series, 2 * dimension * (degree + 1)> all_series = {
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
      permutation[i] = i < dimension * (degree + 1)
                           ? 2 * i
                           : 2 * (i - dimension * (degree + 1)) + 1;
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

template<typename Series,
         int dimension,
         int degree,
         std::size_t... indices,
         std::size_t... periodic_indices>
inline std::array<PoissonSeriesSubspace, 2 * dimension*(degree + 1)>
PoissonSeriesBasisGenerator<Series,
                            dimension,
                            degree,
                            std::index_sequence<indices...>,
                            std::index_sequence<periodic_indices...>>::
    Subspaces(AngularFrequency const& ω, Instant const& origin) {
  return {PoissonSeriesSubspace{
      static_cast<PoissonSeriesSubspace::Coordinate>((periodic_indices / 2) %
                                                     dimension),
      static_cast<PoissonSeriesSubspace::Parity>(
          (1 + periodic_indices + (periodic_indices / 2) / dimension) % 2)}...};
}

}  // namespace internal_poisson_series_basis
}  // namespace numerics
}  // namespace principia
