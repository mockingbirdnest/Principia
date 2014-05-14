#pragma once

#include <stdint.h>

namespace principia {
namespace testing_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar);

double RelativeError(double const expected, double const actual);

int64_t ULPDistance(double const x, double const y);

}  // testing_utilities
}  // namespace principia

#include "testing_utilities/numerics_body.hpp"
