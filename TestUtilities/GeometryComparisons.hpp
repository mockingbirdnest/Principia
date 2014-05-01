#pragma once

#include "TestUtilities.hpp"
#include "QuantityComparisons.hpp"
#include "..\Geometry\R3Element.hpp"
#include "..\Geometry\Grassmann.hpp"

namespace Principia {
namespace TestUtilities {

template<typename Scalar, typename Frame, unsigned int Rank>
void AssertEqual(Geometry::Multivector<Scalar, Frame, Rank> left,
                 Geometry::Multivector<Scalar, Frame, Rank> right);

template<typename Scalar>
void AssertEqual(Geometry::R3Element<Scalar> left,
                 Geometry::R3Element<Scalar> right);


}  // namespace TestUtilities
}  // namespace Principia

#include "GeometryComparisons-body.hpp"
