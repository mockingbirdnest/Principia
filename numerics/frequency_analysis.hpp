
#pragma once

#include <type_traits>

#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/poisson_series.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using geometry::Instant;
using geometry::Interval;
using quantities::AngularFrequency;
using quantities::Primitive;
using quantities::Product;
using quantities::Time;

// In this file Dot is a templated functor that implements the dot product
// between two functors and a weight.  Its declaration must look like:
//
// class Dot {
//  public:
//   ...
//   template<typename LFunction, typename RFunction, typename Weight>
//   Product<std::invoke_result_t<LFunction, Instant>,
//           std::invoke_result_t<RFunction, Instant>>
//   operator()(LFunction const& left,
//              RFunction const& right,
//              Weight const& weight) const;
//   ...
//};
//
// Where the implementation of Dot may assume that Weight is a Poisson series
// returning a double.

// Computes the precise mode of a quasi-periodic |function|, assuming that the
// mode is over the interval |fft_mode| (so named because it has presumably been
// obtained using FFT).  See [Cha95].
template<typename Function,
         int wdegree_,
         typename Dot,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Dot const& dot);

// Computes the Кудрявцев projection of |function| on a basis with angular
// frequency ω and maximum degree |degree_|.  See [Kud07].
// TODO(phl): We really need multiple angular frequencies.
template<int degree_,
         typename Function,
         int wdegree_,
         typename Dot,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
Projection(AngularFrequency const& ω,
           Function const& function,
           PoissonSeries<double, wdegree_, Evaluator> const& weight,
           Dot const& dot);

}  // namespace internal_frequency_analysis

using internal_frequency_analysis::PreciseMode;
using internal_frequency_analysis::Projection;

}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia

#include "numerics/frequency_analysis_body.hpp"
