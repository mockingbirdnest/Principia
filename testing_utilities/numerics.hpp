#pragma once

#include <cstdint>

#include "geometry/r3_element.hpp"
#include "geometry/grassmann.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar);

template<typename T, typename Norm, typename NormType = T>
NormType AbsoluteError(T const& expected, T const& actual,
                       Norm const norm);

// Equivalent to RelativeError(expected, actual, Abs).
double AbsoluteError(double const expected, double const actual);

// Equivalent to RelativeError(expected, actual, Abs<Dimensions>).
template<typename Dimensions>
quantities::Quantity<Dimensions> AbsoluteError(
    quantities::Quantity<Dimensions> const& expected,
    quantities::Quantity<Dimensions> const& actual);

// Uses R3Element.Norm().
template<typename Scalar>
Scalar AbsoluteError(geometry::R3Element<Scalar> const& expected,
                     geometry::R3Element<Scalar> const& actual);

// Uses Multivector.Norm().
template<typename Scalar, typename Frame, unsigned int Rank>
Scalar AbsoluteError(
    geometry::Multivector<Scalar, Frame, Rank> const& expected,
    geometry::Multivector<Scalar, Frame, Rank> const& actual);

template<typename T, typename Norm>
double RelativeError(T const& expected, T const& actual, Norm const norm);

// Equivalent to RelativeError(expected, actual, Abs).
double RelativeError(double const expected, double const actual);

// Equivalent to RelativeError(expected, actual, Abs<Dimensions>).
template<typename Dimensions>
double RelativeError(quantities::Quantity<Dimensions> const& expected,
                     quantities::Quantity<Dimensions> const& actual);

// Uses R3Element.Norm().
template<typename Scalar>
double RelativeError(geometry::R3Element<Scalar> const& expected,
                     geometry::R3Element<Scalar> const& actual);

// Uses Multivector.Norm().
template<typename Scalar, typename Frame, unsigned int Rank>
double RelativeError(geometry::Multivector<Scalar, Frame, Rank> const& expected,
                     geometry::Multivector<Scalar, Frame, Rank> const& actual);

std::int64_t ULPDistance(double const x, double const y);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerics_body.hpp"
