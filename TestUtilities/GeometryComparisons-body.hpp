#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/R3Element.hpp"
#include "TestUtilities/QuantityComparisons.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace test_utilities {

template<typename Scalar, typename Frame, unsigned int Rank>
void AssertEqual(geometry::Multivector<Scalar, Frame, Rank> const& left,
                 geometry::Multivector<Scalar, Frame, Rank> const& right,
                 quantities::Dimensionless const& ε) {
  AssertEqual(left.Coordinates(), right.Coordinates(), ε);
}

template<typename Scalar>
void AssertEqual(geometry::R3Element<Scalar> const& left,
                 geometry::R3Element<Scalar> const& right,
                 quantities::Dimensionless const& ε) {
  AssertEqual(left.x, right.x, ε);
  AssertEqual(left.y, right.y, ε);
  AssertEqual(left.z, right.z, ε);
}

}  // namespace test_utilities
}  // namespace principia
