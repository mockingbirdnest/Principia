#pragma once

#include "TestUtilities.hpp"
#include "QuantityComparisons.hpp"
#include "..\Geometry\R3Element.hpp"
#include "..\Geometry\Grassmann.hpp"

namespace principia {
namespace test_utilities {

template<typename Scalar, typename Frame, unsigned int Rank>
void AssertEqual(geometry::Multivector<Scalar, Frame, Rank> left,
                 geometry::Multivector<Scalar, Frame, Rank> right);

template<typename Scalar>
void AssertEqual(geometry::R3Element<Scalar> left,
                 geometry::R3Element<Scalar> right);

}  // namespace test_utilities
}  // namespace principia

#include "GeometryComparisons-body.hpp"
