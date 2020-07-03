
#pragma once

#include "geometry/named_quantities.hpp"
#include "numerics/poisson_series.hpp"

namespace principia {
namespace numerics {
namespace apodization {
namespace internal_apodization {

using geometry::Instant;

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Dirichlet(Instant const& t_min,
                                              Instant const& t_max);

}  // namespace internal_apodization

using internal_apodization::Dirichlet;

}  // namespace apodization
}  // namespace numerics
}  // namespace principia

#include "numerics/apodization_body.hpp"
