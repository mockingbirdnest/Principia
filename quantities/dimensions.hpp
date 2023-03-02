#pragma once

#include <cstdint>

namespace principia {
namespace quantities {
namespace _dimensions {
namespace internal {

// Dimensionality of physical quantities.  Note that we strongly type angles.
template<std::int64_t LengthExponent,
         std::int64_t MassExponent,
         std::int64_t TimeExponent,
         std::int64_t CurrentExponent,
         std::int64_t TemperatureExponent,
         std::int64_t AmountExponent,
         std::int64_t LuminousIntensityExponent,
         std::int64_t AngleExponent>
struct Dimensions;

// A double by any other name...
using NoDimensions = Dimensions<0, 0, 0, 0, 0, 0, 0, 0>;

// Instantiating this struct asserts at compile time that the template parameter
// can be serialized.
template<typename Dimensions>
struct DimensionsAreSerializable;

// These structs have a |Type| member that is a |Dimensions| suitable for
// the result of the operation applied to argument(s) having the |Dimensions|
// given as template parameter(s).

template<typename Dimensions, int n>
struct DimensionsExponentiationGenerator;

// Only legal if |n| divides the dimensions.
template<typename Dimensions, int n>
struct DimensionsNthRootGenerator;

template<typename LDimensions, typename RDimensions>
struct DimensionsProductGenerator;

template<typename LDimensions, typename RDimensions>
struct DimensionsQuotientGenerator;

}  // namespace internal
}  // namespace _dimensions
}  // namespace quantities
}  // namespace principia

namespace principia::quantities {
using namespace principia::quantities::_dimensions;
}  // namespace principia::quantities

#include "quantities/dimensions_body.hpp"
