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

using geometry::Hilbert;
using geometry::Instant;
using quantities::AngularFrequency;

//TODO(phl):Fix the comments
// A helper struct for generating the Кудрявцев basis, i.e., functions of the
// form tⁿ sin ω t and tⁿ cos ω t properly ordered.  |dimension| is the number
// of multivector dimensions to produce.  |degree| is the maximum degree of tⁿ.

// In this template, the indices encode the degree and the dimension of the
// basis term so that, in the terminology of SeriesGenerator, n (the degree) is
// indices / dimension and d (the dimension index) is indices % dimension.
template<typename Series, int degree>
class PoissonSeriesBasisGenerator {
  using Value = std::invoke_result_t<Series, Instant>;
  static constexpr int dimension = Hilbert<Value>::dimension;

 public:
  // Basis of aperiodic terms.
  static std::array<Series, dimension * (degree + 1)> Basis(
      Instant const& origin);

  // Basis of periodic terms.
  static std::array<Series, 2 * dimension * (degree + 1)> Basis(
      AngularFrequency const& ω,
      Instant const& origin);
};

}  // namespace internal_poisson_series_basis

using internal_poisson_series_basis::PoissonSeriesBasisGenerator;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_basis_body.hpp"
