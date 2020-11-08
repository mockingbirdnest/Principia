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
// form tⁿ sin ω t and tⁿ cos ω t properly ordered.  |dimension| is the number
// of multivector dimensions to produce.  |degree| is the maximum degree of tⁿ.
template<typename Series,
         int dimension,
         int degree,
         typename = std::make_index_sequence<dimension * (degree + 1)>>
struct PoissonSeriesBasisGenerator;

// A |PoissonSeriesSubspace| represents a linear subspace of the space of
// Poisson series.
// The type |PoissonSeriesSubspace| defines an orthogonal decomposition of the
// space of Poisson series, i.e., The space of Poisson series is the orthogonal
// sum of all values of |PoissonSeriesSubspace|.
class PoissonSeriesSubspace {
 public:
  PoissonSeriesSubspace(PoissonSeriesSubspace const&) = default;
  PoissonSeriesSubspace(PoissonSeriesSubspace&&) = default;

  PoissonSeriesSubspace& operator=(PoissonSeriesSubspace const&) = default;
  PoissonSeriesSubspace& operator=(PoissonSeriesSubspace&&) = default;

  // Whether the subspaces |v| and |w| are orthogonal.
  // TODO(egg): when we take parity into account, orthogonality will be defined
  // with respect to an inner product over an interval centred on the origin of
  // the Poisson series, with an even apodization.
  static bool orthogonal(PoissonSeriesSubspace v, PoissonSeriesSubspace w);

 private:
  enum class Coordinate { X = 0, Y = 1, Z = 2 };
  enum class Parity { Even = 0, Odd = 1 };

  PoissonSeriesSubspace(Coordinate coordinate, Parity parity);

  Coordinate coordinate_;
  // The parity of the Poisson series is defined with respect to its origin.
  Parity parity_;

  template<typename Series, int dimension, int degree, typename>
  friend struct PoissonSeriesBasisGenerator;
};

// In this template, the indices encode the degree and the dimension of the
// basis term so that, in the terminology of SeriesGenerator, n (the degree) is
// indices / dimension and d (the dimension index) is indices % dimension.
template<typename Series, int dimension, int degree, std::size_t... indices>
struct PoissonSeriesBasisGenerator<Series,
                                   dimension,
                                   degree,
                                   std::index_sequence<indices...>> {
  // Basis of aperiodic terms.
  static std::array<Series, dimension * (degree + 1)> Basis(
      Instant const& origin);
  // The subspaces to which the above terms belong.
  static std::array<PoissonSeriesSubspace, dimension * (degree + 1)> Subspaces(
      Instant const& origin);

  // Basis of periodic terms.
  static std::array<Series, 2 * dimension * (degree + 1)> Basis(
      AngularFrequency const& ω,
      Instant const& origin);
  // The subspaces to which the above terms belong.
  static std::array<PoissonSeriesSubspace, 2 * dimension * (degree + 1)>
  Subspaces(AngularFrequency const& ω, Instant const& origin);
};

}  // namespace internal_poisson_series_basis

using internal_poisson_series_basis::PoissonSeriesBasisGenerator;
using internal_poisson_series_basis::PoissonSeriesSubspace;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_basis_body.hpp"
