#pragma once

#include "TestUtilities.hpp"
#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"

namespace Principia {
namespace TestUtilities {

template<typename D>
void AssertEqual(Quantities::Quantity<D> const& left,
                 Quantities::Quantity<D> const& right,
                 Quantities::Dimensionless const& ε) {
  AssertEqualWithin(left, right, ε);
}

template<typename D>
void AssertNotEqual(Quantities::Quantity<D> const& left,
                    Quantities::Quantity<D> const& right,
                    Quantities::Dimensionless const& ε) {
  AssertNotEqualWithin(left, right, ε);
}

}  // namespace TestUtilities
}  // namespace Principia
