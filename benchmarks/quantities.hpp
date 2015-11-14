#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;

namespace quantities {

inline void DimensionfulDiscreteCosineTransform(
    not_null<std::vector<Momentum>*> const result);

inline void DoubleDiscreteCosineTransform(
    not_null<std::vector<double>*> const result);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/quantities_body.hpp"
