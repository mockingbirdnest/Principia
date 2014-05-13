#pragma once

#include <CppUnitTest.h>

#include "quantities/dimensionless.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace test_utilities {
// The Microsoft equivalent only takes a wchar_t*.
void WriteLog(std::wstring const& message);
void NewLine();
// Equivalent to Log(message); Newline();
void LogLine(std::wstring const& message);

// The Microsoft equivalent only takes a wchar_t*.
void AssertFalse(bool const test, std::wstring const& message = L"");
void AssertTrue(bool const test, std::wstring const& message = L"");

// The Microsoft equivalent supports errors only for double.
template<typename ValueType>
void AssertEqualWithin(ValueType const& left,
                       ValueType const& right,
                       quantities::Dimensionless const& ε);
template<typename ValueType>
void AssertEqual(ValueType const& left,
                 ValueType const& right);

template<typename ValueType>
void AssertNotEqualWithin(ValueType const& left,
                          ValueType const& right,
                          quantities::Dimensionless const& ε);
template<typename ValueType>
void AssertNotEqual(ValueType const& left,
                    ValueType const& right);

}  // namespace test_utilities
}  // namespace principia

#include "TestUtilities/TestUtilities-body.hpp"
