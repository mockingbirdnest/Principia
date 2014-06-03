#pragma once

#include <vector>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace benchmarks {

inline void DimensionfulDiscreteCosineTransform(
  std::vector<quantities::Momentum>* result);

inline void DoubleDiscreteCosineTransform(
  std::vector<double>* result);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/quantities_body.hpp"
