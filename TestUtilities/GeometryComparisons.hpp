#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/R3Element.hpp"
#include "Quantities/Dimensionless.hpp"
#include "TestUtilities/QuantityComparisons.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace test_utilities {

template<typename Scalar, typename Frame, unsigned int Rank>
void AssertEqual(geometry::Multivector<Scalar, Frame, Rank> const& left,
                 geometry::Multivector<Scalar, Frame, Rank> const& right,
                 quantities::Dimensionless const& ε = 0);

template<typename Scalar>
void AssertEqual(geometry::R3Element<Scalar> const& left,
                 geometry::R3Element<Scalar> const& right,
                 quantities::Dimensionless const& ε = 0);

}  // namespace test_utilities
}  // namespace principia

#include "TestUtilities/GeometryComparisons-body.hpp"
