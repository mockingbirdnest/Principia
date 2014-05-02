#pragma once

#include "TestUtilities.hpp"
#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"

namespace principia {
namespace test_utilities {

template<typename D>
void AssertEqual(quantities::Quantity<D> const& left,
                 quantities::Quantity<D> const& right,
                 quantities::Dimensionless const& ε) {
  AssertEqualWithin(left, right, ε);
}

template<typename D>
void AssertNotEqual(quantities::Quantity<D> const& left,
                    quantities::Quantity<D> const& right,
                    quantities::Dimensionless const& ε) {
  AssertNotEqualWithin(left, right, ε);
}


inline void AssertEqualAbsolute(quantities::Dimensionless const& left,
                                quantities::Dimensionless const& right,
                                quantities::Dimensionless const& ε) {
  std::wstring const message = L"Should be equal within " + ToString(ε, 3) +
                               L" (absolute): " + ToString(left) + L" and " +
                               ToString(right) + L".";
  LogLine(message);
  AssertTrue(Abs(left - right) < ε, message);
  LogLine(L"> Passed!");
}

inline void AssertEqual(quantities::Dimensionless const& left,
                        quantities::Dimensionless const& right,
                        quantities::Dimensionless const& ε) {
  if (left == 0 || right == 0) {
    AssertEqualAbsolute(left, right, ε);
  } else {
    AssertEqualWithin(left, right, ε);
  }
}

inline void AssertNotEqual(quantities::Dimensionless const& left,
                           quantities::Dimensionless const& right,
                           quantities::Dimensionless const& ε) {
  AssertNotEqualWithin(left, right, ε);
}

}  // namespace test_utilities
}  // namespace principia
