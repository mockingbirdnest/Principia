
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
PoissonSeries<double, 0, Evaluator> Hann(Instant const& t_min,
                                         Instant const& t_max) {
  //TODO(phl): Centre them all.
  using Result = PoissonSeries<double, 0, Evaluator>;
  auto const t_mid = Barycentre<Instant, double>({t_min, t_max}, {0.5, 0.5});
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({0.5}, t_mid),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_mid),
                   /*cos=*/typename Result::Polynomial({0.5}, t_mid)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Hamming(Instant const& t_min,
                                            Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(
      typename Result::Polynomial({25.0 / 46.0}, t_min),
      {{ω,
        {/*sin=*/typename Result::Polynomial({0}, t_min),
         /*cos=*/typename Result::Polynomial({-21.0 / 46.0}, t_min)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Blackman(Instant const& t_min,
                                             Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({0.42}, t_min),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.5}, t_min)}},
                 {2 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({0.08}, t_min)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> ExactBlackman(Instant const& t_min,
                                                  Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(
      typename Result::Polynomial({3969.0 / 9304.0}, t_min),
      {{ω,
        {/*sin=*/typename Result::Polynomial({0}, t_min),
         /*cos=*/typename Result::Polynomial({-1155.0 / 2326.0}, t_min)}},
       {2 * ω,
        {/*sin=*/typename Result::Polynomial({0}, t_min),
         /*cos=*/typename Result::Polynomial({715.0 / 9304.0}, t_min)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> Nuttall(Instant const& t_min,
                                            Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({0.355768}, t_min),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.487396}, t_min)}},
                 {2 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({0.144232}, t_min)}},
                 {3 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.012604}, t_min)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> BlackmanNuttall(Instant const& t_min,
                                                    Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({0.3635819}, t_min),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.4891775}, t_min)}},
                 {2 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({0.1365995}, t_min)}},
                 {3 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.0106411}, t_min)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> BlackmanHarris(Instant const& t_min,
                                                   Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({0.35875}, t_min),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.48829}, t_min)}},
                 {2 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({0.14128}, t_min)}},
                 {3 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial({-0.01168}, t_min)}}});
}

template<template<typename, typename, int> class Evaluator>
PoissonSeries<double, 0, Evaluator> ISO18431_2(Instant const& t_min,
                                               Instant const& t_max) {
  using Result = PoissonSeries<double, 0, Evaluator>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  return Result(typename Result::Polynomial({1.0 / 4.63867187}, t_min),
                {{ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial(
                       {-1.93261719 / 4.63867187}, t_min)}},
                 {2 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial(
                       {1.28613281 / 4.63867187}, t_min)}},
                 {3 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial(
                       {-0.38769531 / 4.63867187}, t_min)}},
                 {4 * ω,
                  {/*sin=*/typename Result::Polynomial({0}, t_min),
                   /*cos=*/typename Result::Polynomial(
                       {0.03222656 / 4.63867187}, t_min)}}});
}

}  // namespace internal_apodization
}  // namespace apodization
}  // namespace numerics
}  // namespace principia
