
#pragma once

#include "numerics/poisson_series.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>::PoissonSeries(
    Polynomial const& aperiodic,
    std::map<AngularFrequency, AngularFrequencyPolynomials> const& periodic)
    : aperiodic_(aperiodic), periodic_(periodic) {}

}  // namespace internal_poisson_series
}  // namespace numerics
}  // namespace principia
