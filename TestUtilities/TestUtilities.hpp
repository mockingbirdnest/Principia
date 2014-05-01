#pragma once

#include "CppUnitTest.h"

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"

namespace Principia {
namespace TestUtilities {
// The Microsoft equivalent only takes a wchar_t*.
void WriteLog(std::wstring const& message) {
  Microsoft::VisualStudio::CppUnitTestFramework::Logger::WriteMessage(
      message.c_str());
}
void NewLine() {
  Microsoft::VisualStudio::CppUnitTestFramework::Logger::WriteMessage(L"\n");
}
// Equivalent to Log(message); Newline();
void LogLine(std::wstring const& message) {
  WriteLog(message);
  NewLine();
}
// The Microsoft equivalent only takes a wchar_t*.
void AssertFalse(bool const test, std::wstring const& message = L"") {
  Microsoft::VisualStudio::CppUnitTestFramework::Assert::IsFalse(
      test,
      message.c_str());
}

// The Microsoft equivalent only takes a wchar_t*.
void AssertTrue(bool const test, std::wstring const& message = L"") {
  Microsoft::VisualStudio::CppUnitTestFramework::Assert::IsTrue(
      test,
      message.c_str());
}

// The Microsoft equivalent supports errors only for double.
template<typename ValueType, typename ErrorType>
void AssertEqualWithin(ValueType const& left,
                       ValueType const& right,
                       ErrorType const& ε) {
  std::wstring message = L"Should be equal within " + ToString(ε, 3) +
                         L": " + ToString(left) + L" and " + ToString(right) +
                         L".";
  LogLine(message);
  AssertTrue(left == right || Abs(left / right - 1) < ε, message);
  LogLine(L"> Passed!");
}

template<typename ValueType, typename ErrorType>
void AssertNotEqualWithin(ValueType const& left,
                          ValueType const& right,
                          ErrorType const& ε) {
  std::wstring message = L"Should differ by more than " + ToString(ε, 3) +
                         L": " + ToString(left) + L" and " + ToString(right) +
                         L".";
  LogLine(message);
  AssertTrue(Abs(left / right - 1) > ε, message);
  LogLine(L"> Passed!");
}

}  // namespace TestUtilities
}  // namespace Principia
