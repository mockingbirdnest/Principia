#include "Stdafx.hpp"

#include "QuantityComparisons.hpp"

#include "TestUtilities.hpp"
#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"

namespace Principia {
namespace TestUtilities {

void AssertEqualAbsolute(Quantities::Dimensionless const& left,
                         Quantities::Dimensionless const& right,
                         Quantities::Dimensionless const& ε) {
  std::wstring message = L"Should be equal within " + ToString(ε, 3) +
    L" (absolute): " + ToString(left) + L" and " +
    ToString(right) + L".";
  LogLine(message);
  AssertTrue(Abs(left - right) < ε, message);
  LogLine(L"> Passed!");
}

void AssertEqual(Quantities::Dimensionless const& left,
                 Quantities::Dimensionless const& right,
                 Quantities::Dimensionless const& ε) {
  if (left == 0 || right == 0) { AssertEqualAbsolute(left, right, ε); }
  else {AssertEqualWithin(left, right, ε); }
}

void AssertNotEqual(Quantities::Dimensionless const& left,
                    Quantities::Dimensionless const& right,
                    Quantities::Dimensionless const& ε) {
  AssertNotEqualWithin(left, right, ε);
}

}  // namespace TestUtilities
}  // namespace Principia
