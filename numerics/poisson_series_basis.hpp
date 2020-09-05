#pragma once

#include <array>
#include <utility>

#include "geometry/named_quantities.hpp"
#include "numerics/poisson_series.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series_basis {

using geometry::Instant;
using quantities::AngularFrequency;

// A helper struct for generating the Кудрявцев basis, i.e., functions of the
// form tⁿ sin ω t and tⁿ cos ω t properly ordered.
template<typename Series,
         int dimension,
         typename = std::make_index_sequence<dimension * (Series::degree + 1)>>
struct PoissonSeriesBasisGenerator;

// In this template, the indices encode the degree and the dimension of the
// basis term so that, in the terminology of SeriesGenerator, n (the degree) is
// indices / dimension and d (the dimension index) is indices % dimension.
template<typename Series, int dimension, std::size_t... indices>
struct PoissonSeriesBasisGenerator<Series,
                                   dimension,
                                   std::index_sequence<indices...>> {
  // Basis of aperiodic terms.
  static std::array<Series, dimension * (Series::degree + 1)> Basis(
      Instant const& origin);

  // Basis of periodic terms.
  static std::array<Series, 2 * dimension * (Series::degree + 1)> Basis(
      AngularFrequency const& ω,
      Instant const& origin);
};

}  // namespace internal_poisson_series_basis

using internal_poisson_series_basis::PoissonSeriesBasisGenerator;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_basis_body.hpp"
