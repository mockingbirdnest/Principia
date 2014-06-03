#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include "quantities/dimensionless.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace benchmarks {

#define TRIGGER_DEAD_CODE_ELIMINATION

namespace {

std::size_t const kDimension = 100;

}  // namespace

inline void DimensionfulDiscreteCosineTransform(
    std::vector<quantities::Momentum>* result) {
  using quantities::Cos;
  using quantities::Dimensionless;
  using quantities::Momentum;
  using si::Radian;
  std::vector<Momentum> input(kDimension);
  for (std::size_t i = 0; i < kDimension; ++i) {
    input[i] = i * Momentum::SIUnit();
  }
  result->resize(kDimension);
  Dimensionless sign = 1;
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
    std::vector<double>* result) {
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
      sum += input[n] * std::cos(M_PI / (kDimension - 1) * n * k);
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
