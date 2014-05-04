#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/R3Element.hpp"
#include "TestUtilities/QuantityComparisons.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace test_utilities {

template<typename Scalar, typename Frame, unsigned int Rank>
void AssertEqual(geometry::Multivector<Scalar, Frame, Rank> left,
                 geometry::Multivector<Scalar, Frame, Rank> right) {
  AssertEqual(left.coordinates, right.coordinates);
}

template<typename Scalar>
void AssertEqual(geometry::R3Element<Scalar> left,
                 geometry::R3Element<Scalar> right) {
  AssertEqual(left.x, right.x);
  AssertEqual(left.y, right.y);
  AssertEqual(left.z, right.z);
}

}  // namespace test_utilities
}  // namespace principia
