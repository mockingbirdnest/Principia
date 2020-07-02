
#pragma once

#include "numerics/apodization.hpp"

namespace principia {
namespace numerics {
namespace apodization {
namespace internal_apodization {

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Dirichlet(Instant const& t_min,
                                              Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  return Result(typename Result::Polynomial({1}), {});
}

}  // namespace internal_apodization
}  // namespace apodization
}  // namespace numerics
}  // namespace principia
