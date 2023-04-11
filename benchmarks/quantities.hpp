#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace benchmarks {
namespace _quantities {
namespace internal {

using namespace principia::quantities::_named_quantities;

inline void DimensionfulDiscreteCosineTransform(std::vector<Momentum>& result);

inline void DoubleDiscreteCosineTransform(std::vector<double>& result);

}  // namespace internal

using internal::DimensionfulDiscreteCosineTransform;
using internal::DoubleDiscreteCosineTransform;

}  // namespace _quantities
}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/quantities_body.hpp"
