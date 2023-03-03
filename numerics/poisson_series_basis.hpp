#pragma once

#include <array>
#include <type_traits>
#include <utility>

#include "geometry/hilbert.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/poisson_series.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series_basis {

using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_named_quantities;
using namespace principia::quantities::_named_quantities;

// A |PoissonSeriesSubspace| represents a linear subspace of the space of
// Poisson series.
// The type |PoissonSeriesSubspace| defines an orthogonal decomposition of the
// space of Poisson series, i.e., The space of Poisson series is the orthogonal
// sum of all values of |PoissonSeriesSubspace|.
class PoissonSeriesSubspace {
 public:
  // Whether the subspaces |v| and |w| are orthogonal.
  // TODO(egg): When we take parity into account, orthogonality will be defined
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

  template<typename Series, int degree, int dimension, typename>
  friend struct AperiodicSeriesGenerator;
  template<typename Series, int degree, int dimension, typename>
  friend struct PeriodicSeriesGenerator;
  friend std::ostream& operator<<(std::ostream& out,
                                  PoissonSeriesSubspace const& subspace);
};

// A generator for the Кудрявцев basis, i.e., functions of the
// form tⁿ sin ω t and tⁿ cos ω t properly ordered.  The basis elements have
// type |Series|, which must be free from Quantity.  |degree| is the maximum
// degree of tⁿ.  The basis elements are valid over the interval [t_min, t_max]
// and symmetrical (odd, even) around the midpoint of that interval.
template<typename Series, int degree>
class PoissonSeriesBasisGenerator {
  using Value = std::invoke_result_t<Series, Instant>;
  static constexpr int dimension = Hilbert<Value>::dimension;
  static_assert(std::is_same_v<Value, typename Hilbert<Value>::NormalizedType>,
                "Value type must be free from Quantity");

 public:
  // Basis of aperiodic terms.
  static std::array<Series, dimension * (degree + 1)> Basis(
      Instant const& t_min,
      Instant const& t_max);
  // The subspaces to which the above terms belong.
  static std::array<PoissonSeriesSubspace, dimension * (degree + 1)>
  Subspaces();

  // Basis of periodic terms.
  static std::array<Series, 2 * dimension * (degree + 1)> Basis(
      AngularFrequency const& ω,
      Instant const& t_min,
      Instant const& t_max);
  // The subspaces to which the above terms belong.
  static std::array<PoissonSeriesSubspace, 2 * dimension * (degree + 1)>
  Subspaces(AngularFrequency const& ω);
};

std::ostream& operator<<(std::ostream& out,
                         PoissonSeriesSubspace const& subspace);

}  // namespace internal_poisson_series_basis

using internal_poisson_series_basis::PoissonSeriesBasisGenerator;
using internal_poisson_series_basis::PoissonSeriesSubspace;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_basis_body.hpp"
