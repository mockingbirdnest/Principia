#pragma once

#include "quantities/dimensionless.hpp"
#include "quantities/quantities.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace test_utilities {

template<typename D>
void AssertEqual(quantities::Quantity<D> const& left,
                 quantities::Quantity<D> const& right,
                 quantities::Dimensionless const& ε = 0);

template<typename D>
void AssertNotEqual(quantities::Quantity<D> const& left,
                    quantities::Quantity<D> const& right,
                    quantities::Dimensionless const& ε = 0);

void AssertEqualAbsolute(quantities::Dimensionless const& left,
                         quantities::Dimensionless const& right,
                         quantities::Dimensionless const& ε = 0);

void AssertEqual(quantities::Dimensionless const& left,
                 quantities::Dimensionless const& right,
                 quantities::Dimensionless const& ε = 0);

void AssertNotEqual(quantities::Dimensionless const& left,
                    quantities::Dimensionless const& right,
                    quantities::Dimensionless const& ε = 0);

}  // namespace test_utilities
}  // namespace principia

#include "TestUtilities/QuantityComparisons-body.hpp"
