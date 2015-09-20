#pragma once

#include "benchmarks/quantities.hpp"

#include <cmath>
#include <vector>

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using quantities::Cos;
using quantities::Momentum;
using quantities::SIUnit;
using quantities::si::Radian;

namespace benchmarks {

#define TRIGGER_DEAD_CODE_ELIMINATION

namespace {

std::size_t const kDimension = 100;

}  // namespace

inline void DimensionfulDiscreteCosineTransform(
    not_null<std::vector<quantities::Momentum>*> const result) {
  std::vector<Momentum> input(kDimension);
  for (std::size_t i = 0; i < kDimension; ++i) {
    input[i] = i * SIUnit<Momentum>();
  }
  result->resize(kDimension);
  double sign = 1;
  Momentum sum;
  for (std::size_t k = 0; k < kDimension; ++k, sign *= -1) {
    sum = Momentum();
    for (std::size_t n = 1; n < kDimension - 1; ++n) {
      sum += input[n] * quantities::Cos(π * Radian / (kDimension - 1) * n * k);
    }
#ifdef TRIGGER_DEAD_CODE_ELIMINATION
    (*result)[k] = 0.5 * (input[0] + sign * input[kDimension - 1]);
#else
    (*result)[k] = 0.5 * (input[0] + sign * input[kDimension - 1]) + sum;
#endif
  }
}

inline void DoubleDiscreteCosineTransform(
    not_null<std::vector<double>*> const result) {
  std::vector<double> input(kDimension);
  for (std::size_t i = 0; i < kDimension; ++i) {
    input[i] = i;
  }
  result->resize(kDimension);
  double sign = 1;
  double sum;
  for (std::size_t k = 0; k < kDimension; ++k, sign *= -1) {
    sum = 0;
    for (std::size_t n = 1; n < kDimension - 1; ++n) {
      sum += input[n] * std::cos(π / (kDimension - 1) * n * k);
    }
#ifdef TRIGGER_DEAD_CODE_ELIMINATION
    (*result)[k] = 0.5 * (input[0] + sign * input[kDimension - 1]);
#else
    (*result)[k] = 0.5 * (input[0] + sign * input[kDimension - 1]) + sum;
#endif
  }
}

}  // namespace benchmarks
}  // namespace principia
