
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

//TODO(phl):comments
template<typename Function,
         typename RValue, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    std::function<Primitive<Product<decltype(std::declval<Function>().Evaluate(
                                        std::declval<Instant>())),
                                    RValue>,
                            Time>(
        Function const& left,
        PoissonSeries<RValue, rdegree_, Evaluator> const& right,
        PoissonSeries<double, wdegree_, Evaluator> const& weight)> const& dot);

}  // namespace internal_frequency_analysis

using internal_frequency_analysis::PreciseMode;

}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia

#include "numerics/frequency_analysis_body.hpp"