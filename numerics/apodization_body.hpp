
#pragma once

#include "numerics/apodization.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace apodization {
namespace internal_apodization {

using geometry::Barycentre;
using quantities::AngularFrequency;
using quantities::si::Radian;

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Dirichlet(Instant const& t_min,
                                              Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  return Result(typename Result::Polynomial({1}, t_min), {});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Sine(Instant const& t_min,
                                         Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({0}, t_min),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({1}, t_min),
                   /*cos=*/typename Result::Polynomial({0}, t_min)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Hahn(Instant const& t_min,
                                         Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({0.5}, t_min),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.5}, t_min)}}});
}

}  // namespace internal_apodization
}  // namespace apodization
}  // namespace numerics
}  // namespace principia
