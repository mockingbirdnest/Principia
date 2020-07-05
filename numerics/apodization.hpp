
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

}  // namespace internal_apodization

using internal_apodization::Dirichlet;
using internal_apodization::Hahn;
using internal_apodization::Sine;

}  // namespace apodization
}  // namespace numerics
}  // namespace principia

#include "numerics/apodization_body.hpp"
