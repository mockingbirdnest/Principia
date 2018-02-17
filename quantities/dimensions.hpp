
#pragma once

#include <cstdint>

namespace principia {
namespace quantities {
namespace internal_dimensions {

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

using NoDimensions = Dimensions<0, 0, 0, 0, 0, 0, 0, 0>;

}  // namespace internal_dimensions

using internal_dimensions::Dimensions;
using internal_dimensions::NoDimensions;

}  // namespace quantities
}  // namespace principia

#include "quantities/dimensions_body.hpp"
