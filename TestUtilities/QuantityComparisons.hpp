#pragma once

#include "Quantities/Dimensionless.hpp"
#include "Quantities/Quantities.hpp"

#include "TestUtilities.hpp"

namespace Principia {
namespace TestUtilities {

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

}  // namespace TestUtilities
}  // namespace Principia

#include "QuantityComparisons-body.hpp"
