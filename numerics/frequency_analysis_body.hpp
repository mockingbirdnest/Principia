
#pragma once

#include "numerics/frequency_analysis.hpp"

#include "numerics/root_finders.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using quantities::Square;

template<typename Function,
         typename RValue, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    std::function<Primitive<Product<std::invoke_result_t<Function, Instant>,
                                    RValue>,
                            Time>(
        Function const& left,
        PoissonSeries<RValue, rdegree_, Evaluator> const& right,
        PoissonSeries<double, wdegree_, Evaluator> const& weight)> const& dot) {
  using DotResult =
      Primitive<Product<std::invoke_result_t<Function, Instant>, RValue>, Time>;
  using Degree0 = PoissonSeries<double, 0, Evaluator>;

  auto amplitude = [&dot, &function, &weight](AngularFrequency const& ω) {
    Instant const& t0 = weight.origin();
    Degree0 const sin(typename Degree0::Polynomial({0}, t0),
                      {{ω,
                        {/*sin=*/typename Degree0::Polynomial({1}, t0),
                         /*cos=*/typename Degree0::Polynomial({0}, t0)}}});
    Degree0 const cos(typename Degree0::Polynomial({0}, t0),
                      {{ω,
                        {/*sin=*/typename Degree0::Polynomial({0}, t0),
                         /*cos=*/typename Degree0::Polynomial({1}, t0)}}});
    auto const sin_amplitude = dot(function, sin, weight);
    auto const cos_amplitude = dot(function, cos, weight);
    return sin_amplitude * sin_amplitude + cos_amplitude * cos_amplitude;
  };

  return GoldenSectionSearch(amplitude,
                             fft_mode.min,
                             fft_mode.max,
                             std::greater<Square<DotResult>>());
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia