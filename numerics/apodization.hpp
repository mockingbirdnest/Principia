#pragma once

#include "geometry/instant.hpp"
#include "numerics/poisson_series.hpp"

// The order and terminology in this file follows
// https://en.wikipedia.org/wiki/Window_function.

namespace principia {
namespace numerics {
namespace _apodization {
namespace internal {

using namespace principia::geometry::_instant;
using namespace principia::numerics::_poisson_series;

// ISO 18431-2:2004, section 5.4.
template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> Dirichlet(Instant const& t_min,
                                                 Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> Sine(Instant const& t_min,
                                            Instant const& t_max);

// ISO 18431-2:2004, section 5.2.
template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> Hann(Instant const& t_min,
                                            Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> Hamming(Instant const& t_min,
                                               Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> Blackman(Instant const& t_min,
                                                Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> ExactBlackman(Instant const& t_min,
                                                     Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> Nuttall(Instant const& t_min,
                                               Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> BlackmanNuttall(Instant const& t_min,
                                                       Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> BlackmanHarris(Instant const& t_min,
                                                      Instant const& t_max);

// The flat-top window in Wikipedia is not normalized and comes from Matlab (?).
// We use the normalized function from ISO 18431-2:2004/Cor.1:2008, section 5.3
// instead.
template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, 0, Evaluator> ISO18431_2(Instant const& t_min,
                                                  Instant const& t_max);

}  // namespace internal

using internal::Blackman;
using internal::BlackmanHarris;
using internal::BlackmanNuttall;
using internal::Dirichlet;
using internal::ExactBlackman;
using internal::Hann;
using internal::Hamming;
using internal::ISO18431_2;
using internal::Nuttall;
using internal::Sine;

}  // namespace _apodization
}  // namespace numerics
}  // namespace principia

#include "numerics/apodization_body.hpp"
