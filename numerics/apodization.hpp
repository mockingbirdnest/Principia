
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

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Dirichlet(Instant const& t_min,
                                              Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Sine(Instant const& t_min,
                                         Instant const& t_max);

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Hahn(Instant const& t_min,
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

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> FlatTop(Instant const& t_min,
                                            Instant const& t_max);

}  // namespace internal_apodization

using internal_apodization::Blackman;
using internal_apodization::BlackmanHarris;
using internal_apodization::BlackmanNuttall;
using internal_apodization::Dirichlet;
using internal_apodization::ExactBlackman;
using internal_apodization::FlatTop;
using internal_apodization::Hahn;
using internal_apodization::Hamming;
using internal_apodization::Nuttall;
using internal_apodization::Sine;

}  // namespace apodization
}  // namespace numerics
}  // namespace principia

#include "numerics/apodization_body.hpp"
