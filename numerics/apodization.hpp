
#pragma once

#include "geometry/named_quantities.hpp"
#include "numerics/poisson_series.hpp"

// The order and terminology in this file follows
// https://en.wikipedia.org/wiki/Window_function.

namespace principia {
namespace numerics {
namespace apodization {
namespace internal_apodization {

using geometry::Instant;

// ISO 18431-2:2004, section 5.4.
template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Dirichlet(Instant const& t_min,
                                              Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Sine(Instant const& t_min,
                                         Instant const& t_max);

// ISO 18431-2:2004, section 5.2.
template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Hann(Instant const& t_min,
                                         Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Hamming(Instant const& t_min,
                                            Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Blackman(Instant const& t_min,
                                             Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> ExactBlackman(Instant const& t_min,
                                                  Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Nuttall(Instant const& t_min,
                                            Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> BlackmanNuttall(Instant const& t_min,
                                                    Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> BlackmanHarris(Instant const& t_min,
                                                   Instant const& t_max);

// The flat-top window in Wikipedia is not normalized and comes from Matlab (?).
// We use the normalized function from ISO 18431-2:2004/Cor.1:2008, section 5.3
// instead.
template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> ISO18431_2(Instant const& t_min,
                                               Instant const& t_max);

}  // namespace internal_apodization

using internal_apodization::Blackman;
using internal_apodization::BlackmanHarris;
using internal_apodization::BlackmanNuttall;
using internal_apodization::Dirichlet;
using internal_apodization::ExactBlackman;
using internal_apodization::Hann;
using internal_apodization::Hamming;
using internal_apodization::ISO18431_2;
using internal_apodization::Nuttall;
using internal_apodization::Sine;

}  // namespace apodization
}  // namespace numerics
}  // namespace principia

#include "numerics/apodization_body.hpp"
