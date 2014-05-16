#pragma once

#include <cstdint>

#include "quantities/dimensionless.hpp"

namespace principia {
namespace testing_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar);

template<typename T, typename Norm>
Dimensionless RelativeError(T const& expected, T const& actual,
                            Norm const norm);

std::int64_t ULPDistance(double const x, double const y);

}  // testing_utilities
}  // namespace principia

#include "testing_utilities/numerics_body.hpp"
