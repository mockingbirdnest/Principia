#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/R3Element.hpp"

#include "QuantityComparisons.hpp"
#include "TestUtilities.hpp"

namespace Principia {
namespace TestUtilities {

template<typename Scalar, typename Frame, unsigned int Rank>
void AssertEqual(Geometry::Multivector<Scalar, Frame, Rank> left,
                 Geometry::Multivector<Scalar, Frame, Rank> right) {
  AssertEqual(left.coordinates, right.coordinates);
}

template<typename Scalar>
void AssertEqual(Geometry::R3Element<Scalar> left,
                 Geometry::R3Element<Scalar> right) {
  AssertEqual(left.x, right.x);
  AssertEqual(left.y, right.y);
  AssertEqual(left.z, right.z);
}

}  // namespace TestUtilities
}  // namespace Principia
