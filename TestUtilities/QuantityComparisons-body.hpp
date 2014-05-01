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


inline void AssertEqualAbsolute(Quantities::Dimensionless const& left,
                                Quantities::Dimensionless const& right,
                                Quantities::Dimensionless const& ε) {
  std::wstring message = L"Should be equal within " + ToString(ε, 3) +
    L" (absolute): " + ToString(left) + L" and " +
    ToString(right) + L".";
  LogLine(message);
  AssertTrue(Abs(left - right) < ε, message);
  LogLine(L"> Passed!");
}

inline void AssertEqual(Quantities::Dimensionless const& left,
                        Quantities::Dimensionless const& right,
                        Quantities::Dimensionless const& ε) {
  if (left == 0 || right == 0) { AssertEqualAbsolute(left, right, ε); }
  else {AssertEqualWithin(left, right, ε); }
}

inline void AssertNotEqual(Quantities::Dimensionless const& left,
                           Quantities::Dimensionless const& right,
                           Quantities::Dimensionless const& ε) {
  AssertNotEqualWithin(left, right, ε);
}

}  // namespace TestUtilities
}  // namespace Principia
