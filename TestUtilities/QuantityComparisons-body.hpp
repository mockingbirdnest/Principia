#pragma once

#include <string>

#include "quantities/dimensionless.hpp"
#include "quantities/quantities.hpp"
#include "TestUtilities/TestUtilities.hpp"

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
  std::string const message = "Should be equal within " + ToString(ε, 3) +
                              " (absolute): " + ToString(left) + " and " +
                              ToString(right) + ".";
  LogLine(message);
  AssertTrue(Abs(left - right) <= ε, message);
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
