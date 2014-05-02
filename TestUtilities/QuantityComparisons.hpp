#pragma once

#include "TestUtilities.hpp"
#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"

namespace principia {
namespace test_utilities {

Quantities::Dimensionless const tolerance =  1e-14;

template<typename D>
void AssertEqual(Quantities::Quantity<D> const& left,
                 Quantities::Quantity<D> const& right,
                 Quantities::Dimensionless const& ε = tolerance);

template<typename D>
void AssertNotEqual(Quantities::Quantity<D> const& left,
                    Quantities::Quantity<D> const& right,
                    Quantities::Dimensionless const& ε = tolerance);

void AssertEqualAbsolute(Quantities::Dimensionless const& left,
                         Quantities::Dimensionless const& right,
                         Quantities::Dimensionless const& ε = tolerance);

void AssertEqual(Quantities::Dimensionless const& left,
                 Quantities::Dimensionless const& right,
                 Quantities::Dimensionless const& ε = tolerance);

void AssertNotEqual(Quantities::Dimensionless const& left,
                    Quantities::Dimensionless const& right,
                    Quantities::Dimensionless const& ε = tolerance);

}  // namespace test_utilities
}  // namespace principia

#include "QuantityComparisons-body.hpp"
