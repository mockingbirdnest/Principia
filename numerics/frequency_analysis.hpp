
#pragma once

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

template<typename Function, typename RValue,
         int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
struct Dotty {
  Primitive<Product<decltype(std::declval<Function>()()), RValue>, Time>
  Dot(Function const& left,
      PoissonSeries<RValue, rdegree_, Evaluator> const& right,
      PoissonSeries<double, wdegree_, Evaluator> const& weight,
      Instant const& t_min,
      Instant const& t_max);
};

template<typename Function,
         template<typename,
                  typename,
                  int,
                  int,
                  template<typename, typename, int>
                  class>
         class Dot>
AngularFrequency PreciseMode(Interval<AngularFrequency> const& fft_mode,
                             Function const& function);

}  // namespace internal_frequency_analysis


}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia

#include "numerics/frequency_analysis_body.hpp"