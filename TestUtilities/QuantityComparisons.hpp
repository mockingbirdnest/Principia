#pragma once

#include "Quantities/Dimensionless.hpp"
#include "Quantities/Quantities.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace test_utilities {

quantities::Dimensionless const tolerance =  1e-14;

template<typename D>
void AssertEqual(quantities::Quantity<D> const& left,
                 quantities::Quantity<D> const& right,
                 quantities::Dimensionless const& ε = tolerance);

template<typename D>
void AssertNotEqual(quantities::Quantity<D> const& left,
                    quantities::Quantity<D> const& right,
                    quantities::Dimensionless const& ε = tolerance);

void AssertEqualAbsolute(quantities::Dimensionless const& left,
                         quantities::Dimensionless const& right,
                         quantities::Dimensionless const& ε = tolerance);

void AssertEqual(quantities::Dimensionless const& left,
                 quantities::Dimensionless const& right,
                 quantities::Dimensionless const& ε = tolerance);

void AssertNotEqual(quantities::Dimensionless const& left,
                    quantities::Dimensionless const& right,
                    quantities::Dimensionless const& ε = tolerance);

}  // namespace test_utilities
}  // namespace principia

#include "QuantityComparisons-body.hpp"
