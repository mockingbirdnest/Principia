#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace quantities {

inline void DimensionfulDiscreteCosineTransform(std::vector<Momentum>& result);

inline void DoubleDiscreteCosineTransform(std::vector<double>& result);

}  // namespace quantities
}  // namespace principia

#include "benchmarks/quantities_body.hpp"
