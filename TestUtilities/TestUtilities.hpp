#pragma once

#include "CppUnitTest.h"

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"

namespace principia {
namespace test_utilities {
// The Microsoft equivalent only takes a wchar_t*.
void WriteLog(std::wstring const& message);
void NewLine();
// Equivalent to Log(message); Newline();
void LogLine(std::wstring const& message);

// The Microsoft equivalent only takes a wchar_t*.
void AssertTrue(bool const test, std::wstring const& message = L"");

// The Microsoft equivalent supports errors only for double.
template<typename ValueType, typename ErrorType>
void AssertEqualWithin(ValueType const& left,
                       ValueType const& right,
                       ErrorType const& ε);

template<typename ValueType, typename ErrorType>
void AssertNotEqualWithin(ValueType const& left,
                          ValueType const& right,
                          ErrorType const& ε);

}  // namespace test_utilities
}  // namespace principia

#include "TestUtilities-body.hpp"
