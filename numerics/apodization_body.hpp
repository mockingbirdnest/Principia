
#pragma once

#include "numerics/apodization.hpp"

#include "geometry/barycentre_calculator.hpp"

namespace principia {
namespace numerics {
namespace apodization {
namespace internal_apodization {

using geometry::Barycentre;

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Dirichlet(Instant const& t_min,
                                              Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  auto const midpoint = Barycentre<Instant, double>({t_min, t_max}, {0.5, 0.5});
  return Result(typename Result::Polynomial({1}, midpoint), {});
}

}  // namespace internal_apodization
}  // namespace apodization
}  // namespace numerics
}  // namespace principia
