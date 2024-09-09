#pragma once

#include <cstdint>

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace _numerics {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::geometry::_r3_element;
using namespace principia::quantities::_quantities;

template<typename Scalar>
double DoubleValue(Scalar const& scalar);

template<typename T, typename NormType>
NormType AbsoluteError(T const& expected, T const& actual,
                       NormType (T::* norm)() const);

template<typename T, typename NormType, typename NormArg>
NormType AbsoluteError(T const& expected, T const& actual,
                       NormType (*norm)(NormArg const));

// Equivalent to AbsoluteError(expected, actual, &Abs).
double AbsoluteError(double expected, double actual);

// Equivalent to AbsoluteError(expected, actual, &Abs<Dimensions>).
template<typename Dimensions>
Quantity<Dimensions> AbsoluteError(Quantity<Dimensions> const& expected,
                                   Quantity<Dimensions> const& actual);

// Uses `R3Element::Norm`.
template<typename Scalar>
Scalar AbsoluteError(R3Element<Scalar> const& expected,
                     R3Element<Scalar> const& actual);

// Uses Multivector.Norm().
template<typename Scalar, typename Frame, int rank>
Scalar AbsoluteError(
    Multivector<Scalar, Frame, rank> const& expected,
    Multivector<Scalar, Frame, rank> const& actual);

// Uses the underlying multivector.
template<typename Scalar, typename Frame>
Scalar AbsoluteError(
    Point<Multivector<Scalar, Frame, 1>> const& expected,
    Point<Multivector<Scalar, Frame, 1>> const& actual);

// Uses the underlying multivector.
template<typename Scalar>
Scalar AbsoluteError(Point<Scalar> const& expected,
                     Point<Scalar> const& actual);

template<typename T, typename NormType>
double RelativeError(T const& expected, T const& actual,
                     NormType (T::* norm)() const);

template<typename T, typename NormType, typename NormArg>
double RelativeError(T const& expected, T const& actual,
                     NormType (*norm)(NormArg const));

// Equivalent to RelativeError(expected, actual, &Abs).
double RelativeError(double expected, double actual);

// Equivalent to RelativeError(expected, actual, &Abs<Dimensions>).
template<typename Dimensions>
double RelativeError(Quantity<Dimensions> const& expected,
                     Quantity<Dimensions> const& actual);

// Uses `R3Element::Norm`.
template<typename Scalar>
double RelativeError(R3Element<Scalar> const& expected,
                     R3Element<Scalar> const& actual);

// Uses Multivector.Norm().
template<typename Scalar, typename Frame, int rank>
double RelativeError(Multivector<Scalar, Frame, rank> const& expected,
                     Multivector<Scalar, Frame, rank> const& actual);

}  // namespace internal

using internal::AbsoluteError;
using internal::DoubleValue;
using internal::RelativeError;

}  // namespace _numerics
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerics_body.hpp"
