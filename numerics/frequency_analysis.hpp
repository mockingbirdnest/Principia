
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

// A function that implements the dot product between a |Function| and a Poisson
// series with an apodization function |weight| function that is itself a
// Poisson series.  |Function| must be a functor taking an Instant.
//TODO(phl):comment
//template<typename Function,
//         typename RValue, int rdegree_, int wdegree_,
//         template<typename, typename, int> class Evaluator>
//using DotProduct =
//    std::function<Product<std::invoke_result_t<Function, Instant>, RValue>(
//        Function const& left,
//        PoissonSeries<RValue, rdegree_, Evaluator> const& right,
//        PoissonSeries<double, wdegree_, Evaluator> const& weight)>;

// Computes the precise mode of a quasi-periodic |function|, assuming that the
// mode is over the interval |fft_mode| (so named because it has presumably been
// obtained using FFT).  See [Cha95].
template<template<typename, typename, typename> class Dot,
         typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight);

// Computes the Кудрявцев projection of |function| on a basis with angular
// frequency ω and maximum degree |degree_|.  See [Kud07].
// TODO(phl): We really need multiple angular frequencies.
template<template<typename, typename, typename> class Dot,
         int degree_,
         typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
Projection(
    AngularFrequency const& ω,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight);

}  // namespace internal_frequency_analysis

using internal_frequency_analysis::PreciseMode;
using internal_frequency_analysis::Projection;

}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia

#include "numerics/frequency_analysis_body.hpp"
